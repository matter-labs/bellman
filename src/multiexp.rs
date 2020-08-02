use crate::pairing::{
    CurveAffine,
    CurveProjective,
    Engine
};

use crate::pairing::ff::{
    PrimeField,
    Field,
    PrimeFieldRepr,
    ScalarEngine};

use std::sync::Arc;
use super::source::*;
use std::future::{Future};
use std::task::{Context, Poll};
use std::pin::{Pin};

extern crate futures;

use self::futures::future::{join_all, JoinAll};
use self::futures::executor::block_on;

use super::worker::{Worker, WorkerFuture};

use super::SynthesisError;

use cfg_if;

use hwloc2::{Topology, ObjectType, CpuBindFlags, CpuSet};
/// This genious piece of code works in the following way:
/// - choose `c` - the bit length of the region that one thread works on
/// - make `2^c - 1` buckets and initialize them with `G = infinity` (that's equivalent of zero)
/// - there is no bucket for "zero" cause it's not necessary
/// - go over the pairs `(base, scalar)`
/// - for each scalar calculate `scalar % 2^c` and add the base (without any multiplications!) to the 
/// corresponding bucket
/// - at the end each bucket will have an accumulated value that should be multiplied by the corresponding factor
/// between `1` and `2^c - 1` to get the right value
/// - here comes the first trick - you don't need to do multiplications at all, just add all the buckets together
/// starting from the first one `(a + b + c + ...)` and than add to the first sum another sum of the form
/// `(b + c + d + ...)`, and than the third one `(c + d + ...)`, that will result in the proper prefactor infront of every
/// accumulator, without any multiplication operations at all
/// - that's of course not enough, so spawn the next thread
/// - this thread works with the same bit width `c`, but SKIPS lowers bits completely, so it actually takes values
/// in the form `(scalar >> c) % 2^c`, so works on the next region
/// - spawn more threads until you exhaust all the bit length
/// - you will get roughly `[bitlength / c] + 1` inaccumulators
/// - double the highest accumulator enough times, add to the next one, double the result, add the next accumulator, continue
/// 
/// Demo why it works:
/// ```text
///     a * G + b * H = (a_2 * (2^c)^2 + a_1 * (2^c)^1 + a_0) * G + (b_2 * (2^c)^2 + b_1 * (2^c)^1 + b_0) * H
/// ```
/// - make buckets over `0` labeled coefficients
/// - make buckets over `1` labeled coefficients
/// - make buckets over `2` labeled coefficients
/// - accumulators over each set of buckets will have an implicit factor of `(2^c)^i`, so before summing thme up
/// "higher" accumulators must be doubled `c` times
///
#[cfg(not(feature = "nightly"))]
fn multiexp_inner<Q, D, G, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
    skip: u32,
    c: u32,
    handle_trivial: bool
) -> WorkerFuture< <G as CurveAffine>::Projective, SynthesisError>
    where for<'a> &'a Q: QueryDensity,
          D: Send + Sync + 'static + Clone + AsRef<Q>,
          G: CurveAffine,
          S: SourceBuilder<G>
{
    // Perform this region of the multiexp
    let this = {
        // let bases = bases.clone();
        // let exponents = exponents.clone();
        // let density_map = density_map.clone();

        // This is a Pippenger’s algorithm
        pool.compute(move || {
            // Accumulate the result
            let mut acc = G::Projective::zero();

            // Build a source for the bases
            let mut bases = bases.new();

            // Create buckets to place remainders s mod 2^c,
            // it will be 2^c - 1 buckets (no bucket for zeroes)

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
            let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

            // Sort the bases into buckets
            for (&exp, density) in exponents.iter().zip(density_map.as_ref().iter()) {
                // Go over density and exponents
                if density {
                    if exp == zero {
                        bases.skip(1)?;
                    } else if exp == one {
                        if handle_trivial {
                            bases.add_assign_mixed(&mut acc)?;
                        } else {
                            bases.skip(1)?;
                        }
                    } else {
                        // Place multiplication into the bucket: Separate s * P as 
                        // (s/2^c) * P + (s mod 2^c) P
                        // First multiplication is c bits less, so one can do it,
                        // sum results from different buckets and double it c times,
                        // then add with (s mod 2^c) P parts
                        let mut exp = exp;
                        exp.shr(skip);
                        let exp = exp.as_ref()[0] % (1 << c);

                        if exp != 0 {
                            bases.add_assign_mixed(&mut buckets[(exp - 1) as usize])?;
                        } else {
                            bases.skip(1)?;
                        }
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.into_iter().rev() {
                running_sum.add_assign(&exp);
                acc.add_assign(&running_sum);
            }

            Ok(acc)
        })
    };

    this
}


cfg_if! {
    if #[cfg(feature = "nightly")] {
        #[inline(always)]
        fn multiexp_inner_impl<Q, D, G, S>(
            pool: &Worker,
            bases: S,
            density_map: D,
            exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
            skip: u32,
            c: u32,
            handle_trivial: bool
        ) -> WorkerFuture< <G as CurveAffine>::Projective, SynthesisError>
            where for<'a> &'a Q: QueryDensity,
                D: Send + Sync + 'static + Clone + AsRef<Q>,
                G: CurveAffine,
                S: SourceBuilder<G>
        {
            // multiexp_inner_with_prefetch(pool, bases, density_map, exponents, skip, c, handle_trivial)
            multiexp_inner_with_prefetch_stable(pool, bases, density_map, exponents, skip, c, handle_trivial)
        }
    } else {
        #[inline(always)]
        fn multiexp_inner_impl<Q, D, G, S>(
            pool: &Worker,
            bases: S,
            density_map: D,
            exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
            skip: u32,
            c: u32,
            handle_trivial: bool
        ) -> WorkerFuture< <G as CurveAffine>::Projective, SynthesisError>
            where for<'a> &'a Q: QueryDensity,
                D: Send + Sync + 'static + Clone + AsRef<Q>,
                G: CurveAffine,
                S: SourceBuilder<G>
        {
            // multiexp_inner(pool, bases, density_map, exponents, skip, c, handle_trivial)
            multiexp_inner_with_prefetch_stable(pool, bases, density_map, exponents, skip, c, handle_trivial)
        }
    }  
}

#[cfg(feature = "nightly")]
extern crate prefetch;

#[cfg(feature = "nightly")]
fn multiexp_inner_with_prefetch<Q, D, G, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
    skip: u32,
    c: u32,
    handle_trivial: bool
) -> WorkerFuture< <G as CurveAffine>::Projective, SynthesisError>
    where for<'a> &'a Q: QueryDensity,
          D: Send + Sync + 'static + Clone + AsRef<Q>,
          G: CurveAffine,
          S: SourceBuilder<G>
{
    use prefetch::prefetch::*;
    // Perform this region of the multiexp
    let this = {
        // This is a Pippenger’s algorithm
        pool.compute(move || {
            // Accumulate the result
            let mut acc = G::Projective::zero();

            // Build a source for the bases
            let mut bases = bases.new();

            // Create buckets to place remainders s mod 2^c,
            // it will be 2^c - 1 buckets (no bucket for zeroes)

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
            let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();
            let padding = Arc::new(vec![zero]);

            let mask = 1 << c;

            // Sort the bases into buckets
            for ((&exp, &next_exp), density) in exponents.iter()
                        .zip(exponents.iter().skip(1).chain(padding.iter()))
                        .zip(density_map.as_ref().iter()) {
                // no matter what happens - prefetch next bucket
                if next_exp != zero && next_exp != one {
                    let mut next_exp = next_exp;
                    next_exp.shr(skip);
                    let next_exp = next_exp.as_ref()[0] % mask;
                    if next_exp != 0 {
                        let p: *const <G as CurveAffine>::Projective = &buckets[(next_exp - 1) as usize];
                        prefetch::<Write, High, Data, _>(p);
                    }
                    
                }
                // Go over density and exponents
                if density {
                    if exp == zero {
                        bases.skip(1)?;
                    } else if exp == one {
                        if handle_trivial {
                            bases.add_assign_mixed(&mut acc)?;
                        } else {
                            bases.skip(1)?;
                        }
                    } else {
                        // Place multiplication into the bucket: Separate s * P as 
                        // (s/2^c) * P + (s mod 2^c) P
                        // First multiplication is c bits less, so one can do it,
                        // sum results from different buckets and double it c times,
                        // then add with (s mod 2^c) P parts
                        let mut exp = exp;
                        exp.shr(skip);
                        let exp = exp.as_ref()[0] % mask;

                        if exp != 0 {
                            bases.add_assign_mixed(&mut buckets[(exp - 1) as usize])?;
                        } else {
                            bases.skip(1)?;
                        }
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.into_iter().rev() {
                running_sum.add_assign(&exp);
                acc.add_assign(&running_sum);
            }

            Ok(acc)
        })
    };
    
    this
}

fn multiexp_inner_with_prefetch_stable<Q, D, G, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
    skip: u32,
    c: u32,
    handle_trivial: bool
) -> WorkerFuture< <G as CurveAffine>::Projective, SynthesisError>
    where for<'a> &'a Q: QueryDensity,
          D: Send + Sync + 'static + Clone + AsRef<Q>,
          G: CurveAffine,
          S: SourceBuilder<G>
{
    // Perform this region of the multiexp
    let this = {
        let bases = bases.clone();
        let exponents = exponents.clone();
        let density_map = density_map.clone();

        // This is a Pippenger’s algorithm
        pool.compute(move || {
            // Accumulate the result
            let mut acc = G::Projective::zero();

            // Build a source for the bases
            let mut bases = bases.new();

            // Create buckets to place remainders s mod 2^c,
            // it will be 2^c - 1 buckets (no bucket for zeroes)

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
            let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();
            let padding = Arc::new(vec![zero]);

            let mask = 1 << c;

            // Sort the bases into buckets
            for ((&exp, &next_exp), density) in exponents.iter()
                        .zip(exponents.iter().skip(1).chain(padding.iter()))
                        .zip(density_map.as_ref().iter()) {
                // no matter what happens - prefetch next bucket
                if next_exp != zero && next_exp != one {
                    let mut next_exp = next_exp;
                    next_exp.shr(skip);
                    let next_exp = next_exp.as_ref()[0] % mask;
                    if next_exp != 0 {
                        let p: *const <G as CurveAffine>::Projective = &buckets[(next_exp - 1) as usize];
                        crate::prefetch::prefetch_l3_pointer(p);
                    }
                    
                }
                // Go over density and exponents
                if density {
                    if exp == zero {
                        bases.skip(1)?;
                    } else if exp == one {
                        if handle_trivial {
                            bases.add_assign_mixed(&mut acc)?;
                        } else {
                            bases.skip(1)?;
                        }
                    } else {
                        // Place multiplication into the bucket: Separate s * P as 
                        // (s/2^c) * P + (s mod 2^c) P
                        // First multiplication is c bits less, so one can do it,
                        // sum results from different buckets and double it c times,
                        // then add with (s mod 2^c) P parts
                        let mut exp = exp;
                        exp.shr(skip);
                        let exp = exp.as_ref()[0] % mask;

                        if exp != 0 {
                            bases.add_assign_mixed(&mut buckets[(exp - 1) as usize])?;
                        } else {
                            bases.skip(1)?;
                        }
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.into_iter().rev() {
                running_sum.add_assign(&exp);
                acc.add_assign(&running_sum);
            }

            Ok(acc)
        })
    };

    this
}


/// Perform multi-exponentiation. The caller is responsible for ensuring the
/// query size is the same as the number of exponents.
pub fn future_based_multiexp<G: CurveAffine>(
    pool: &Worker,
    bases: Arc<Vec<G>>,
    exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>
    // bases: &[G],
    // exponents: Arc<Vec<<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr>>
) -> ChunksJoiner< <G as CurveAffine>::Projective >
{
    assert!(exponents.len() <= bases.len());
    let c = if exponents.len() < 32 {
        3u32
    } else {
        let mut width = (f64::from(exponents.len() as u32)).ln().ceil() as u32;
        let mut num_chunks = <G::Scalar as PrimeField>::NUM_BITS / width;
        if <G::Scalar as PrimeField>::NUM_BITS % width != 0 {
            num_chunks += 1;
        }

        if num_chunks < pool.cpus as u32 {
            width = <G::Scalar as PrimeField>::NUM_BITS / (pool.cpus as u32);
            if <G::Scalar as PrimeField>::NUM_BITS % (pool.cpus as u32) != 0 {
                width += 1;
            }
        }
        
        width
    };

    let mut skip = 0;
    let mut futures = Vec::with_capacity((<G::Engine as ScalarEngine>::Fr::NUM_BITS / c + 1) as usize);

    while skip < <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        let chunk_future = if skip == 0 {
            future_based_dense_multiexp_imlp(pool, bases.clone(), exponents.clone(), 0, c, true)
        } else {
            future_based_dense_multiexp_imlp(pool, bases.clone(), exponents.clone(), skip, c, false)
        };

        futures.push(chunk_future);
        skip += c;
    }

    let join = join_all(futures);

    ChunksJoiner {
        join,
        c
    } 
}


fn future_based_dense_multiexp_imlp<G: CurveAffine>(
    pool: &Worker,
    bases: Arc<Vec<G>>,
    exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
    skip: u32,
    c: u32,
    handle_trivial: bool
) -> WorkerFuture< <G as CurveAffine>::Projective, SynthesisError>
{
    // Perform this region of the multiexp
    let this = {
        let bases = bases.clone();
        let exponents = exponents.clone();
        let bases = bases.clone();

        // This is a Pippenger’s algorithm
        pool.compute(move || {
            // Accumulate the result
            let mut acc = G::Projective::zero();

            // Create buckets to place remainders s mod 2^c,
            // it will be 2^c - 1 buckets (no bucket for zeroes)

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
            let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();
            let padding = Arc::new(vec![zero]);

            let mask = 1 << c;

            // Sort the bases into buckets
            for ((&exp, base), &next_exp) in exponents.iter()
                        .zip(bases.iter())
                        .zip(exponents.iter().skip(1).chain(padding.iter())) {
                // no matter what happens - prefetch next bucket
                if next_exp != zero && next_exp != one {
                    let mut next_exp = next_exp;
                    next_exp.shr(skip);
                    let next_exp = next_exp.as_ref()[0] % mask;
                    if next_exp != 0 {
                        let p: *const <G as CurveAffine>::Projective = &buckets[(next_exp - 1) as usize];
                        crate::prefetch::prefetch_l3_pointer(p);
                    }
                    
                }
                // Go over density and exponents
                if exp == zero {
                    continue
                } else if exp == one {
                    if handle_trivial {
                        acc.add_assign_mixed(base);
                    } else {
                        continue
                    }
                } else {
                    // Place multiplication into the bucket: Separate s * P as 
                    // (s/2^c) * P + (s mod 2^c) P
                    // First multiplication is c bits less, so one can do it,
                    // sum results from different buckets and double it c times,
                    // then add with (s mod 2^c) P parts
                    let mut exp = exp;
                    exp.shr(skip);
                    let exp = exp.as_ref()[0] % mask;

                    if exp != 0 {
                        (&mut buckets[(exp - 1) as usize]).add_assign_mixed(base);
                    } else {
                        continue;
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.into_iter().rev() {
                running_sum.add_assign(&exp);
                acc.add_assign(&running_sum);
            }

            Ok(acc)
        })
    };

    this
}

/// Perform multi-exponentiation. The caller is responsible for ensuring the
/// query size is the same as the number of exponents.
pub fn multiexp<Q, D, G, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr>>
) -> ChunksJoiner< <G as CurveAffine>::Projective >
    where for<'a> &'a Q: QueryDensity,
          D: Send + Sync + 'static + Clone + AsRef<Q>,
          G: CurveAffine,
          S: SourceBuilder<G>
{
    let c = if exponents.len() < 32 {
        3u32
    } else {
        (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    if let Some(query_size) = density_map.as_ref().get_query_size() {
        // If the density map has a known query size, it should not be
        // inconsistent with the number of exponents.

        assert!(query_size == exponents.len());
    }

    let mut skip = 0;
    let mut futures = Vec::with_capacity((<G::Engine as ScalarEngine>::Fr::NUM_BITS / c + 1) as usize);

    while skip < <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        let chunk_future = if skip == 0 {
            multiexp_inner_impl(pool, bases.clone(), density_map.clone(), exponents.clone(), 0, c, true)
        } else {
            multiexp_inner_impl(pool, bases.clone(), density_map.clone(), exponents.clone(), skip, c, false)
        };

        futures.push(chunk_future);
        skip += c;
    }

    let join = join_all(futures);

    ChunksJoiner {
        join,
        c
    } 
}

pub struct ChunksJoiner<G: CurveProjective> {
    join: JoinAll< WorkerFuture<G, SynthesisError> >,
    c: u32
}

impl<G: CurveProjective> Future for ChunksJoiner<G> {
    type Output = Result<G, SynthesisError>;

    fn poll(self: Pin<&mut Self>, cx: &mut Context) -> Poll<Self::Output>
    {
        let c = self.as_ref().c;
        let join = unsafe { self.map_unchecked_mut(|s| &mut s.join) };
        match join.poll(cx) {
            Poll::Ready(v) => {
                let v = join_chunks(v, c);
                return Poll::Ready(v);
            },
            Poll::Pending => {
                return Poll::Pending;
            }
        }
    }
}

impl<G: CurveProjective> ChunksJoiner<G> {
    pub fn wait(self) -> <Self as Future>::Output {
        block_on(self)
    }
}

fn join_chunks<G: CurveProjective>
    (chunks: Vec<Result<G, SynthesisError>>, c: u32) -> Result<G, SynthesisError> {
    if chunks.len() == 0 {
        return Ok(G::zero());
    }

    let mut iter = chunks.into_iter().rev();
    let higher = iter.next().expect("is some chunk result");
    let mut higher = higher?;

    for chunk in iter {
        let this = chunk?;
        for _ in 0..c {
            higher.double();
        }

        higher.add_assign(&this);
    }

    Ok(higher)
}


/// Perform multi-exponentiation. The caller is responsible for ensuring that
/// the number of bases is the same as the number of exponents.
#[allow(dead_code)]
pub fn dense_multiexp<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr]
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{
    if exponents.len() != bases.len() {
        return Err(SynthesisError::AssignmentMissing);
    }
    // do some heuristics here
    // we proceed chunks of all points, and all workers do the same work over 
    // some scalar width, so to have expected number of additions into buckets to 1
    // we have to take log2 from the expected chunk(!) length
    let c = if exponents.len() < 32 {
        3u32
    } else {
        let chunk_len = pool.get_chunk_size(exponents.len());
        (f64::from(chunk_len as u32)).ln().ceil() as u32

        // (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    // dense_multiexp_inner_unrolled_with_prefetch(pool, bases, exponents, 0, c, true)
    dense_multiexp_inner(pool, bases, exponents, 0, c, true)
}

// Get thread id from libc
fn get_thread_id() -> libc::pthread_t {
    unsafe { libc::pthread_self() }
}

// Get core nums
fn get_core_num(topo: &Arc<std::sync::Mutex<hwloc2::Topology>>) -> usize{
    let topo_rc = topo.clone();
    let topo_locked = topo_rc.lock().unwrap();
    (*topo_locked)
        .objects_with_type(&ObjectType::Core)
        .unwrap()
        .len()
}

// Load the `CpuSet` for the given core index.
fn cpuset_for_core(topology: &Topology, idx: usize) -> CpuSet {
    let cores = (*topology).objects_with_type(&ObjectType::Core).unwrap();
    match cores.get(idx) {
        Some(val) => val.cpuset().unwrap(),
        None => panic!("No Core found with id {}", idx),
    }
}

// Bind thread to core
fn bind_thread(
    child_topo: &Arc<std::sync::Mutex<hwloc2::Topology>>,
    idx: usize) {
    // Get the current thread id and lock the topology to use.
    let tid = get_thread_id();
    let mut locked_topo = child_topo.lock().unwrap();

    // Thread binding before explicit set.
    let before = locked_topo.get_cpubind_for_thread(tid, CpuBindFlags::CPUBIND_THREAD);

    // load the cpuset for the given core index.
    let mut bind_to = cpuset_for_core(&*locked_topo, idx);

    // Get only one logical processor (in case the core is SMT/hyper-threaded).
    bind_to.singlify();

    // Set the binding.
    locked_topo
        .set_cpubind_for_thread(tid, bind_to, CpuBindFlags::CPUBIND_THREAD)
        .unwrap();

    // Thread binding after explicit set.
    let after = locked_topo.get_cpubind_for_thread(tid, CpuBindFlags::CPUBIND_THREAD);
    println!("Thread {:?}: Before {:?}, After {:?}", tid, before, after);
}

fn dense_multiexp_inner<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr],
    mut skip: u32,
    c: u32,
    handle_trivial: bool
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{   
    use std::sync::{Mutex};
    // Perform this region of the multiexp. We use a different strategy - go over region in parallel,
    // then over another region, etc. No Arc required
    let this = {
        // let mask = (1u64 << c) - 1u64;
        let this_region = Mutex::new(<G as CurveAffine>::Projective::zero());
        let arc = Arc::new(this_region);

        let topo = Arc::new(Mutex::new(Topology::new().unwrap()));

        // Grab the number of cores.
        let num_cores = get_core_num(&topo);
        println!("Found {} cores.", num_cores);

        pool.scope(bases.len(), |scope, chunk| {
            let mut core_idx = 0;
            for (base, exp) in bases.chunks(chunk).zip(exponents.chunks(chunk)) {
                let this_region_rwlock = arc.clone();
                // let handle = 

                let child_topo = topo.clone();

                scope.spawn(move |_| {

                    // binding thread to specific core
                    bind_thread(&child_topo, core_idx % num_cores);

                    let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];
                    // Accumulate the result
                    let mut acc = G::Projective::zero();
                    let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
                    let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

                    for (base, &exp) in base.iter().zip(exp.iter()) {
                        // let index = (exp.as_ref()[0] & mask) as usize;

                        // if index != 0 {
                        //     buckets[index - 1].add_assign_mixed(base);
                        // }

                        // exp.shr(c as u32);

                        if exp != zero {
                            if exp == one {
                                if handle_trivial {
                                    acc.add_assign_mixed(base);
                                }
                            } else {
                                let mut exp = exp;
                                exp.shr(skip);
                                let exp = exp.as_ref()[0] % (1 << c);
                                if exp != 0 {
                                    buckets[(exp - 1) as usize].add_assign_mixed(base);
                                }
                            }
                        }
                    }

                    // buckets are filled with the corresponding accumulated value, now sum
                    let mut running_sum = G::Projective::zero();
                    for exp in buckets.into_iter().rev() {
                        running_sum.add_assign(&exp);
                        acc.add_assign(&running_sum);
                    }

                    let mut guard = match this_region_rwlock.lock() {
                        Ok(guard) => guard,
                        Err(_) => {
                            panic!("poisoned!"); 
                            // poisoned.into_inner()
                        }
                    };

                    (*guard).add_assign(&acc);
                });

                core_idx += 1;
            }
        });

        let this_region = Arc::try_unwrap(arc).unwrap();
        let this_region = this_region.into_inner().unwrap();

        this_region
    };

    skip += c;

    if skip >= <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        // There isn't another region, and this will be the highest region
        return Ok(this);
    } else {
        // next region is actually higher than this one, so double it enough times
        let mut next_region = dense_multiexp_inner(
            pool, bases, exponents, skip, c, false).unwrap();
        for _ in 0..c {
            next_region.double();
        }

        next_region.add_assign(&this);

        return Ok(next_region);
    }
}

#[allow(dead_code)]
pub fn dense_unrolled_multiexp_with_prefetch<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr]
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{
    if exponents.len() != bases.len() {
        return Err(SynthesisError::AssignmentMissing);
    }
    // do some heuristics here
    // we proceed chunks of all points, and all workers do the same work over 
    // some scalar width, so to have expected number of additions into buckets to 1
    // we have to take log2 from the expected chunk(!) length
    let c = if exponents.len() < 32 {
        3u32
    } else {
        let chunk_len = pool.get_chunk_size(exponents.len());
        (f64::from(chunk_len as u32)).ln().ceil() as u32

        // (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    dense_multiexp_inner_unrolled_with_prefetch(pool, bases, exponents, 0, c, true)
}

#[allow(dead_code)]
fn dense_multiexp_inner_unrolled_with_prefetch<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr],
    mut skip: u32,
    c: u32,
    handle_trivial: bool
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{   
    const UNROLL_BY: usize = 8;

    use std::sync::{Mutex};
    // Perform this region of the multiexp. We use a different strategy - go over region in parallel,
    // then over another region, etc. No Arc required
    let this = {
        let mask = (1u64 << c) - 1u64;
        let this_region = Mutex::new(<G as CurveAffine>::Projective::zero());
        let arc = Arc::new(this_region);

        pool.scope(bases.len(), |scope, chunk| {
            for (bases, exp) in bases.chunks(chunk).zip(exponents.chunks(chunk)) {
                let this_region_rwlock = arc.clone();
                // let handle = 
                scope.spawn(move |_| {
                    let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];
                    // Accumulate the result
                    let mut acc = G::Projective::zero();
                    let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
                    let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

                    let unrolled_steps = bases.len() / UNROLL_BY;
                    let remainder = bases.len() % UNROLL_BY;

                    let mut offset = 0;
                    for _ in 0..unrolled_steps {
                        // [0..7]
                        for i in 0..UNROLL_BY {
                            crate::prefetch::prefetch_l3_pointer(&bases[offset+i] as *const _);
                            crate::prefetch::prefetch_l3_pointer(&exp[offset+i] as *const _);
                        }

                        // offset + [0..6]
                        for i in 0..(UNROLL_BY-1) {
                            let this_exp = exp[offset+i];
                            let mut next_exp = exp[offset+i+1];
                            let base = &bases[offset+i];

                            if this_exp != zero {
                                if this_exp == one {
                                    if handle_trivial {
                                        acc.add_assign_mixed(base);
                                    }
                                } else {
                                    let mut this_exp = this_exp;
                                    this_exp.shr(skip);
                                    let this_exp = this_exp.as_ref()[0] & mask;
                                    if this_exp != 0 {
                                        buckets[(this_exp - 1) as usize].add_assign_mixed(base);
                                    }
                                }
                            }

                            {
                                next_exp.shr(skip);
                                let next_exp = next_exp.as_ref()[0] & mask;
                                if next_exp != 0 {
                                    crate::prefetch::prefetch_l3_pointer(&buckets[(next_exp - 1) as usize] as *const _);
                                }
                            }
                        }

                        // offset + 7
                        let this_exp = exp[offset+(UNROLL_BY-1)];
                        let base = &bases[offset+(UNROLL_BY-1)];

                        if this_exp != zero {
                            if this_exp == one {
                                if handle_trivial {
                                    acc.add_assign_mixed(base);
                                }
                            } else {
                                let mut this_exp = this_exp;
                                this_exp.shr(skip);
                                let this_exp = this_exp.as_ref()[0] & mask;
                                if this_exp != 0 {
                                    buckets[(this_exp - 1) as usize].add_assign_mixed(base);
                                }
                            }
                        }

                        // go into next region
                        offset += UNROLL_BY;
                    }

                    for _ in 0..remainder {
                        let this_exp = exp[offset];
                        let base = &bases[offset];

                        if this_exp != zero {
                            if this_exp == one {
                                if handle_trivial {
                                    acc.add_assign_mixed(base);
                                }
                            } else {
                                let mut this_exp = this_exp;
                                this_exp.shr(skip);
                                let this_exp = this_exp.as_ref()[0] & mask;
                                if this_exp != 0 {
                                    buckets[(this_exp - 1) as usize].add_assign_mixed(base);
                                }
                            }
                        }

                        offset += 1;
                    }

                    // buckets are filled with the corresponding accumulated value, now sum
                    let mut running_sum = G::Projective::zero();
                    for exp in buckets.into_iter().rev() {
                        running_sum.add_assign(&exp);
                        acc.add_assign(&running_sum);
                    }

                    let mut guard = match this_region_rwlock.lock() {
                        Ok(guard) => guard,
                        Err(_) => {
                            panic!("poisoned!"); 
                            // poisoned.into_inner()
                        }
                    };

                    (*guard).add_assign(&acc);
                });
        
            }
        });

        let this_region = Arc::try_unwrap(arc).unwrap();
        let this_region = this_region.into_inner().unwrap();

        this_region
    };

    skip += c;

    if skip >= <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        // There isn't another region, and this will be the highest region
        return Ok(this);
    } else {
        // next region is actually higher than this one, so double it enough times
        let mut next_region = dense_multiexp_inner_unrolled_with_prefetch(
            pool, bases, exponents, skip, c, false).unwrap();
        for _ in 0..c {
            next_region.double();
        }

        next_region.add_assign(&this);

        return Ok(next_region);
    }
}


#[allow(dead_code)]
pub fn dense_multiexp_with_manual_unrolling<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr]
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{
    if exponents.len() != bases.len() {
        return Err(SynthesisError::AssignmentMissing);
    }
    // do some heuristics here
    // we proceed chunks of all points, and all workers do the same work over 
    // some scalar width, so to have expected number of additions into buckets to 1
    // we have to take log2 from the expected chunk(!) length
    let c = if exponents.len() < 32 {
        3u32
    } else {
        let chunk_len = pool.get_chunk_size(exponents.len());
        (f64::from(chunk_len as u32)).ln().ceil() as u32

        // (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    dense_multiexp_with_manual_unrolling_impl(pool, bases, exponents, 0, c, true)
    // dense_multiexp_with_manual_unrolling_impl_2(pool, bases, exponents, 0, c, true)
}


#[allow(dead_code)]
fn dense_multiexp_with_manual_unrolling_impl<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr],
    mut skip: u32,
    c: u32,
    handle_trivial: bool
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{   
    const UNROLL_BY: usize = 1024;

    use std::sync::{Mutex};
    // Perform this region of the multiexp. We use a different strategy - go over region in parallel,
    // then over another region, etc. No Arc required
    let this = {
        let mask = (1u64 << c) - 1u64;
        let this_region = Mutex::new(<G as CurveAffine>::Projective::zero());
        let arc = Arc::new(this_region);

        pool.scope(bases.len(), |scope, chunk| {
            for (bases, exp) in bases.chunks(chunk).zip(exponents.chunks(chunk)) {
                let this_region_rwlock = arc.clone();
                // let handle = 
                scope.spawn(move |_| {
                    // make buckets for ALL exponents including 0 and 1
                    let mut buckets = vec![<G as CurveAffine>::Projective::zero(); 1 << c];

                    // let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
                    let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

                    let mut this_chunk_exponents = [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr::default(); UNROLL_BY];
                    let mut next_chunk_exponents = [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr::default(); UNROLL_BY];

                    let mut this_chunk_bases = [G::zero(); UNROLL_BY];
                    let mut next_chunk_bases = [G::zero(); UNROLL_BY];

                    let unrolled_steps = bases.len() / UNROLL_BY;
                    assert!(unrolled_steps >= 2);
                    let remainder = bases.len() % UNROLL_BY;
                    assert_eq!(remainder, 0);

                    // first step is manually unrolled

                    // manually copy to the stack
                    let mut start_idx = 0;
                    let mut end_idx = UNROLL_BY;
                    this_chunk_exponents.copy_from_slice(&exp[start_idx..end_idx]);
                    this_chunk_bases.copy_from_slice(&bases[start_idx..end_idx]);

                    start_idx += UNROLL_BY;
                    end_idx += UNROLL_BY;
                    next_chunk_exponents.copy_from_slice(&exp[start_idx..end_idx]);
                    next_chunk_bases.copy_from_slice(&bases[start_idx..end_idx]);

                    let mut intra_chunk_idx = 0;

                    let mut previous_exponent_index = 0;
                    let mut previous_base = G::zero();

                    let mut this_exponent_index = 0;
                    let mut this_base = G::zero();

                    let this_exp = this_chunk_exponents[intra_chunk_idx];

                    if this_exp == one {
                        if handle_trivial {
                            this_exponent_index = 1;
                        }
                    } else {
                        let mut this_exp = this_exp;
                        this_exp.shr(skip);
                        let this_exp = this_exp.as_ref()[0] & mask;
                        this_exponent_index = this_exp as usize;
                    }

                    this_base = this_chunk_bases[intra_chunk_idx];

                    previous_base = this_base;
                    previous_exponent_index = this_exponent_index;

                    crate::prefetch::prefetch_l2_pointer(&buckets[previous_exponent_index] as *const _);

                    intra_chunk_idx += 1;

                    // now we can roll

                    for _ in 1..(unrolled_steps-1) {
                        while intra_chunk_idx < UNROLL_BY {
                            // add what was processed in a previous step
                            (&mut buckets[previous_exponent_index]).add_assign_mixed(&previous_base);

                            let this_exp = this_chunk_exponents[intra_chunk_idx];

                            if this_exp == one {
                                if handle_trivial {
                                    this_exponent_index = 1;
                                }
                            } else {
                                let mut this_exp = this_exp;
                                this_exp.shr(skip);
                                let this_exp = this_exp.as_ref()[0] & mask;
                                this_exponent_index = this_exp as usize;
                            }

                            this_base = this_chunk_bases[intra_chunk_idx];

                            previous_base = this_base;
                            previous_exponent_index = this_exponent_index;

                            crate::prefetch::prefetch_l2_pointer(&buckets[previous_exponent_index] as *const _);

                            intra_chunk_idx += 1;
                        }

                        // swap and read next chunk

                        this_chunk_bases = next_chunk_bases;
                        this_chunk_exponents = next_chunk_exponents;

                        start_idx += UNROLL_BY;
                        end_idx += UNROLL_BY;
                        next_chunk_exponents.copy_from_slice(&exp[start_idx..end_idx]);
                        next_chunk_bases.copy_from_slice(&bases[start_idx..end_idx]);

                        intra_chunk_idx = 0;
                    }

                    // process the last one
                    {
                        while intra_chunk_idx < UNROLL_BY {
                            // add what was processed in a previous step
                            (&mut buckets[previous_exponent_index]).add_assign_mixed(&previous_base);

                            let this_exp = this_chunk_exponents[intra_chunk_idx];

                            if this_exp == one {
                                if handle_trivial {
                                    this_exponent_index = 1;
                                }
                            } else {
                                let mut this_exp = this_exp;
                                this_exp.shr(skip);
                                let this_exp = this_exp.as_ref()[0] & mask;
                                this_exponent_index = this_exp as usize;
                            }

                            this_base = this_chunk_bases[intra_chunk_idx];

                            previous_base = this_base;
                            previous_exponent_index = this_exponent_index;

                            crate::prefetch::prefetch_l2_pointer(&buckets[previous_exponent_index] as *const _);

                            intra_chunk_idx += 1;
                        }

                        // very last addition
                        (&mut buckets[previous_exponent_index]).add_assign_mixed(&previous_base);
                    }

                    let _: Vec<_> = buckets.drain(..1).collect();

                    let acc: Vec<_> = buckets.drain(..1).collect();
                    let mut acc = acc[0];

                    // buckets are filled with the corresponding accumulated value, now sum
                    let mut running_sum = G::Projective::zero();
                    for exp in buckets.into_iter().rev() {
                        running_sum.add_assign(&exp);
                        acc.add_assign(&running_sum);
                    }

                    let mut guard = match this_region_rwlock.lock() {
                        Ok(guard) => guard,
                        Err(_) => {
                            panic!("poisoned!"); 
                            // poisoned.into_inner()
                        }
                    };

                    (*guard).add_assign(&acc);
                });
        
            }
        });

        let this_region = Arc::try_unwrap(arc).unwrap();
        let this_region = this_region.into_inner().unwrap();

        this_region
    };

    skip += c;

    if skip >= <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        // There isn't another region, and this will be the highest region
        return Ok(this);
    } else {
        // next region is actually higher than this one, so double it enough times
        let mut next_region = dense_multiexp_with_manual_unrolling_impl(
            pool, bases, exponents, skip, c, false).unwrap();
        for _ in 0..c {
            next_region.double();
        }

        next_region.add_assign(&this);

        return Ok(next_region);
    }
}


#[allow(dead_code)]
fn dense_multiexp_with_manual_unrolling_impl_2<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: & [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr],
    mut skip: u32,
    c: u32,
    _handle_trivial: bool
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{   
    // we assume that a single memory fetch is around 10-12 ns, so before any operation
    // we ideally should prefetch a memory unit for a next operation
    const CACHE_BY: usize = 1024;

    use std::sync::{Mutex};
    // Perform this region of the multiexp. We use a different strategy - go over region in parallel,
    // then over another region, etc. No Arc required
    let this = {
        let mask = (1u64 << c) - 1u64;
        let this_region = Mutex::new(<G as CurveAffine>::Projective::zero());
        let arc = Arc::new(this_region);

        pool.scope(bases.len(), |scope, chunk| {
            for (bases, exp) in bases.chunks(chunk).zip(exponents.chunks(chunk)) {
                let this_region_rwlock = arc.clone();
                // let handle = 
                scope.spawn(move |_| {
                    // make buckets for ALL exponents including 0 and 1
                    let mut buckets = vec![<G as CurveAffine>::Projective::zero(); 1 << c];

                    // let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
                    // let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

                    let mut exponents_chunk = [<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr::default(); CACHE_BY];
                    let mut bases_chunk = [G::zero(); CACHE_BY];

                    let unrolled_steps = bases.len() / CACHE_BY;
                    assert!(unrolled_steps >= 2);
                    let remainder = bases.len() % CACHE_BY;
                    assert_eq!(remainder, 0);

                    use std::ptr::NonNull;

                    let mut basket_pointers_to_process = [(NonNull::< <G as CurveAffine>::Projective>::dangling(), G::zero()); CACHE_BY];

                    let basket_pointer = buckets.as_mut_ptr();

                    let mut start_idx = 0;
                    let mut end_idx = CACHE_BY;

                    for _ in 0..(unrolled_steps-1) {
                        exponents_chunk.copy_from_slice(&exp[start_idx..end_idx]);
                        bases_chunk.copy_from_slice(&bases[start_idx..end_idx]);

                        let mut bucket_idx = 0;

                        for (e, b) in exponents_chunk.iter().zip(bases_chunk.iter()) {
                            let mut this_exp = *e;
                            this_exp.shr(skip);
                            let this_exp = (this_exp.as_ref()[0] & mask) as usize;
                            if this_exp != 0 {
                                let ptr = unsafe { NonNull::new_unchecked(basket_pointer.add(this_exp)) };
                                basket_pointers_to_process[bucket_idx] = (ptr, *b);
                                bucket_idx += 1;
                            }
                        }

                        for i in 0..bucket_idx {
                            crate::prefetch::prefetch_l1_pointer(basket_pointers_to_process[i].0.as_ptr() as *const _);
                        }

                        crate::prefetch::prefetch_l2_pointer(&bases[end_idx] as *const _);
                        crate::prefetch::prefetch_l2_pointer(&bases[end_idx+1] as *const _);

                        for i in 0..bucket_idx {
                            let (mut ptr, to_add) = basket_pointers_to_process[i];
                            let point_ref: &mut _ = unsafe { ptr.as_mut()};
                            point_ref.add_assign_mixed(&to_add);
                        }

                        start_idx += CACHE_BY;
                        end_idx += CACHE_BY;
                    }

                    drop(basket_pointer);

                    let _: Vec<_> = buckets.drain(..1).collect();

                    let acc: Vec<_> = buckets.drain(..1).collect();
                    let mut acc = acc[0];

                    // buckets are filled with the corresponding accumulated value, now sum
                    let mut running_sum = G::Projective::zero();
                    for exp in buckets.into_iter().rev() {
                        running_sum.add_assign(&exp);
                        acc.add_assign(&running_sum);
                    }

                    let mut guard = match this_region_rwlock.lock() {
                        Ok(guard) => guard,
                        Err(_) => {
                            panic!("poisoned!"); 
                            // poisoned.into_inner()
                        }
                    };

                    (*guard).add_assign(&acc);
                });
        
            }
        });

        let this_region = Arc::try_unwrap(arc).unwrap();
        let this_region = this_region.into_inner().unwrap();

        this_region
    };

    skip += c;

    if skip >= <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        // There isn't another region, and this will be the highest region
        return Ok(this);
    } else {
        // next region is actually higher than this one, so double it enough times
        let mut next_region = dense_multiexp_with_manual_unrolling_impl_2(
            pool, bases, exponents, skip, c, false).unwrap();
        for _ in 0..c {
            next_region.double();
        }

        next_region.add_assign(&this);

        return Ok(next_region);
    }
}


/// Perform multi-exponentiation. The caller is responsible for ensuring that
/// the number of bases is the same as the number of exponents.
#[allow(dead_code)]
pub fn dense_multiexp_consume<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: Vec<<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr>
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{
    if exponents.len() != bases.len() {
        return Err(SynthesisError::AssignmentMissing);
    }
    let c = if exponents.len() < 32 {
        3u32
    } else {
        (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    dense_multiexp_inner_consume(pool, bases, exponents, c)
}

fn dense_multiexp_inner_consume<G: CurveAffine>(
    pool: &Worker,
    bases: & [G],
    exponents: Vec<<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr>,
    c: u32,
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
{   
    // spawn exactly required number of threads at the time, not more
    // each thread mutates part of the exponents and walks over the same range of bases

    use std::sync::mpsc::{channel};

    let (tx, rx) = channel();

    pool.scope(bases.len(), |scope, chunk| {
        for (base, exp) in bases.chunks(chunk).zip(exponents.chunks(chunk)) {
            let tx = tx.clone();
            scope.spawn(move |_| {
                let mut skip = 0;

                let mut result = G::Projective::zero();

                let mut buckets = vec![<G as CurveAffine>::Projective::zero(); 1 << c];

                let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
                // let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

                let padding = Some(<G::Engine as ScalarEngine>::Fr::zero().into_repr());
                let mask: u64 = (1 << c) - 1;

                loop {
                    let mut next_bucket_index = (exp[0].as_ref()[0] & mask) as usize;
                    let exp_next_constant_iter = exp.iter().skip(1);
                    // let this_exp_to_use = exp.iter();

                    let mut acc = G::Projective::zero();

                    // for ((base, &this_exp_to_use), &next_exp_to_prefetch) in base.iter()
                    //                         .zip(this_exp_to_use)
                    //                         .zip(exp_next_constant_iter.chain(padding.iter()))
                    //     {
                    for (base, &next_exp_to_prefetch) in base.iter()
                            .zip(exp_next_constant_iter.chain(padding.iter()))
                    {
                        let this_bucket_index = next_bucket_index;

                        {
                            // if next_exp_to_prefetch != zero && next_exp_to_prefetch != one {
                            if next_exp_to_prefetch != zero {
                                let mut e = next_exp_to_prefetch;
                                e.shr(skip);
                                next_bucket_index = (next_exp_to_prefetch.as_ref()[0] & mask) as usize;

                                if next_bucket_index > 0 {
                                    // crate::prefetch::prefetch_l3(&buckets[next_bucket_index]);
                                    crate::prefetch::prefetch_l3_pointer(&buckets[next_bucket_index] as *const _);
                                }
                            } else {
                                next_bucket_index = 0;
                            }
                        }

                        if this_bucket_index > 0 {
                            buckets[this_bucket_index].add_assign_mixed(base);
                        }

                        // // now add base to the bucket that we've 
                        // if this_bucket_index > 1 {
                        //     buckets[this_bucket_index].add_assign_mixed(base);
                        // } else {
                        //     acc.add_assign_mixed(base);
                        // }
                    }

                    // buckets are filled with the corresponding accumulated value, now sum
                    let mut running_sum = G::Projective::zero();
                    // now start from the last one and add
                    for exp in buckets.iter().skip(1).rev() {
                        running_sum.add_assign(&exp);
                        acc.add_assign(&running_sum);
                    }

                    for _ in 0..skip {
                        acc.double();
                    }

                    result.add_assign(&acc);

                    skip += c;
                    
                    if skip >= <G::Engine as ScalarEngine>::Fr::NUM_BITS {
                        // next chunk is the last one
                        tx.send(result).unwrap();

                        break;
                    } else {
                        buckets.truncate(0);
                        buckets.resize(1 << c, <G as CurveAffine>::Projective::zero());
                    }
                }
            });
        }
    });

    // do something with rx

    let mut result = <G as CurveAffine>::Projective::zero();

    for value in rx.try_iter() {
        result.add_assign(&value);
    }

    Ok(result)
}


#[test]
fn test_new_multiexp_with_bls12() {
    fn naive_multiexp<G: CurveAffine>(
        bases: Arc<Vec<G>>,
        exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>
    ) -> G::Projective
    {
        assert_eq!(bases.len(), exponents.len());

        let mut acc = G::Projective::zero();

        for (base, exp) in bases.iter().zip(exponents.iter()) {
            acc.add_assign(&base.mul(*exp));
        }

        acc
    }

    use rand::{self, Rand};
    use crate::pairing::bls12_381::Bls12;

    use self::futures::executor::block_on;

    const SAMPLES: usize = 1 << 14;

    let rng = &mut rand::thread_rng();
    let v = Arc::new((0..SAMPLES).map(|_| <Bls12 as ScalarEngine>::Fr::rand(rng).into_repr()).collect::<Vec<_>>());
    let g = Arc::new((0..SAMPLES).map(|_| <Bls12 as Engine>::G1::rand(rng).into_affine()).collect::<Vec<_>>());

    let naive = naive_multiexp(g.clone(), v.clone());

    let pool = Worker::new();

    let fast = block_on(
        multiexp(
            &pool,
            (g, 0),
            FullDensity,
            v
        )
    ).unwrap();

    assert_eq!(naive, fast);
}

#[test]
#[ignore]
fn test_new_multexp_speed_with_bn256() {
    use rand::{self, Rand};
    use crate::pairing::bn256::Bn256;
    use num_cpus;

    let cpus = num_cpus::get();
    const SAMPLES: usize = 1 << 22;

    let rng = &mut rand::thread_rng();
    let v = Arc::new((0..SAMPLES).map(|_| <Bn256 as ScalarEngine>::Fr::rand(rng).into_repr()).collect::<Vec<_>>());
    let g = Arc::new((0..SAMPLES).map(|_| <Bn256 as Engine>::G1::rand(rng).into_affine()).collect::<Vec<_>>());

    let pool = Worker::new();

    use self::futures::executor::block_on;

    let start = std::time::Instant::now();

    let _fast = block_on(
        multiexp(
            &pool,
            (g, 0),
            FullDensity,
            v
        )
    ).unwrap();


    let duration_ns = start.elapsed().as_nanos() as f64;
    println!("Elapsed {} ns for {} samples", duration_ns, SAMPLES);
    let time_per_sample = duration_ns/(SAMPLES as f64);
    println!("Tested on {} samples on {} CPUs with {} ns per multiplication", SAMPLES, cpus, time_per_sample);
}


#[test]
fn test_dense_multiexp_vs_new_multiexp() {
    use rand::{XorShiftRng, SeedableRng, Rand, Rng};
    use crate::pairing::bn256::Bn256;
    use num_cpus;

    // const SAMPLES: usize = 1 << 22;
    const SAMPLES: usize = 1 << 16;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let v = (0..SAMPLES).map(|_| <Bn256 as ScalarEngine>::Fr::rand(rng).into_repr()).collect::<Vec<_>>();
    let g = (0..SAMPLES).map(|_| <Bn256 as Engine>::G1::rand(rng).into_affine()).collect::<Vec<_>>();

    println!("Done generating test points and scalars");

    let pool = Worker::new();

    let start = std::time::Instant::now();

    let dense = dense_multiexp(
        &pool, &g, &v.clone()).unwrap();

    let duration_ns = start.elapsed().as_nanos() as f64;
    println!("{} ns for dense for {} samples", duration_ns, SAMPLES);

    use self::futures::executor::block_on;

    let start = std::time::Instant::now();

    let sparse = block_on(
        multiexp(
            &pool,
            (Arc::new(g), 0),
            FullDensity,
            Arc::new(v)
        )
    ).unwrap();

    let duration_ns = start.elapsed().as_nanos() as f64;
    println!("{} ns for sparse for {} samples", duration_ns, SAMPLES);

    assert_eq!(dense, sparse);
}


#[test]
fn test_bench_sparse_multiexp() {
    use rand::{XorShiftRng, SeedableRng, Rand, Rng};
    use crate::pairing::bn256::Bn256;
    use num_cpus;

    const SAMPLES: usize = 1 << 22;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let v = (0..SAMPLES).map(|_| <Bn256 as ScalarEngine>::Fr::rand(rng).into_repr()).collect::<Vec<_>>();
    let g = (0..SAMPLES).map(|_| <Bn256 as Engine>::G1::rand(rng).into_affine()).collect::<Vec<_>>();

    println!("Done generating test points and scalars");

    let pool = Worker::new();
    let start = std::time::Instant::now();

    let _sparse = multiexp(
        &pool,
        (Arc::new(g), 0),
        FullDensity,
        Arc::new(v)
    ).wait().unwrap();

    let duration_ns = start.elapsed().as_nanos() as f64;
    println!("{} ms for sparse for {} samples", duration_ns/1000.0f64, SAMPLES);
}

#[test]
fn test_bench_dense_consuming_multiexp() {
    use rand::{XorShiftRng, SeedableRng, Rand, Rng};
    use crate::pairing::bn256::Bn256;
    use num_cpus;

    const SAMPLES: usize = 1 << 20;
    let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let v = (0..SAMPLES).map(|_| <Bn256 as ScalarEngine>::Fr::rand(rng).into_repr()).collect::<Vec<_>>();
    let g = (0..SAMPLES).map(|_| <Bn256 as Engine>::G1::rand(rng).into_affine()).collect::<Vec<_>>();

    println!("Done generating test points and scalars");

    let pool = Worker::new();

    let g = Arc::new(g);
    let v = Arc::new(v);

    let start = std::time::Instant::now();

    let _sparse = multiexp(
        &pool,
        (g.clone(), 0),
        FullDensity,
        v.clone()
    ).wait().unwrap();

    println!("{:?} for sparse for {} samples", start.elapsed(), SAMPLES);

    let g = Arc::try_unwrap(g).unwrap();
    let v = Arc::try_unwrap(v).unwrap();

    let start = std::time::Instant::now();

    let _dense = dense_multiexp_consume(
        &pool,
        &g,
        v
    ).unwrap();

    println!("{:?} for dense for {} samples", start.elapsed(), SAMPLES);
}