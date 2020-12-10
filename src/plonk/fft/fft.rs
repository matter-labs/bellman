use crate::pairing::ff::PrimeField;
use crate::worker::*;
use crate::gpu::{self, GPUError, LockedFFTKernel};

use crate::pairing::Engine;


use crate::pairing::bn256::{self, Bn256,Fr};
use crate::domain::{EvaluationDomain, Scalar};
use crate::domain::{self, best_fft_gpu};
use crate::plonk::fft::cooley_tukey_ntt::{bitreverse as ntt_bitreverse, CTPrecomputations, BitReversedOmegas, OmegasInvBitreversed};

use crate::plonk::domains::*;

use log::info;
use crate::plonk::fft::cooley_tukey_ntt;


fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow+1)) <= num {
        pow += 1;
    }

    pow
}


pub(crate) fn bit_rev_best_ct_ntt_2_best_fft_gpu<F: PrimeField, P: CTPrecomputations<F>>(
    a: &mut [F],
    worker: &Worker,
    log_n: u32,
    precomputed_omegas: &P
)
{
    info!("begin bit_rev_best_ct_ntt_2_best_fft_gpu");
    assert_eq!(a.len(), (1<<log_n) as usize);

    let omega_1 = BitReversedOmegas::<Fr>::new_for_domain_size(4).omegas[1];
    let omega_1_r = unsafe{std::mem::transmute::<&Fr, &F>(&omega_1)};
    let omegas_inv_1 = OmegasInvBitreversed::<Fr>::new_for_domain_size(4).omegas[1];
    let omegas_inv_1_r = unsafe{std::mem::transmute::<&Fr, &F>(&omegas_inv_1)};

    let omegas_precom_1 = precomputed_omegas.bit_reversed_omegas()[1];

    let poly_size = a.len();

    let domain = Domain::<F>::new_for_size(poly_size as u64).unwrap();
    let omega = domain.generator;
    let omega_inv = domain.generator.inverse().expect("must exist");

    if *omega_1_r == omegas_precom_1 {
        info!("use best_fft_bn256_gpu with omega");
        best_fft_bn256_gpu(a, worker, &omega, log_n, None);
    }
    else if  *omegas_inv_1_r == omegas_precom_1 {
        info!("use best_fft_bn256_gpu with omega_inv");
        best_fft_bn256_gpu(a, worker, &omega_inv, log_n, None);
    }
    else {
        info!("use native cooley_tukey_ntt::best_ct_ntt");
        cooley_tukey_ntt::best_ct_ntt(a, worker, log_n, Some(worker.cpus), precomputed_omegas);
    }

    let log_n = log_n as usize;
    for k in 0..poly_size {
        let rk = ntt_bitreverse(k, log_n);
        if k < rk {
            a.swap(rk, k);
        }
    }
    info!("end bit_rev_best_ct_ntt_2_best_fft_gpu");
}

// can only be used to omega, cannot be used omegaInv
pub(crate) fn best_ct_ntt_2_best_fft_gpu<F: PrimeField>(
    a: &mut [F],
    worker: &Worker,
    poly_size: usize,
)
{
    let domain = Domain::<Fr>::new_for_size(poly_size as u64).unwrap();
    let omega = domain.generator;
    let omega_r = unsafe{std::mem::transmute::<&Fr, &F>(&omega)};
    let log_n = domain.power_of_two as u32;
    best_fft_bn256_gpu(a, worker, omega_r, log_n, None);

    let log_n = log_n as usize;
    for k in 0..poly_size {
        let rk = ntt_bitreverse(k, log_n);
        if k < rk {
            a.swap(rk, k);
        }
    }

}

pub(crate) fn best_fft_bn256_gpu<F: PrimeField>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    use_cpus_hint: Option<usize>
)
{
    info!("begin best_fft_bn256_gpu");
    let mut fft_kern = Some(LockedFFTKernel::<Bn256>::new(log_n as usize, false));
    let a = unsafe { std::mem::transmute::<&mut [F], &mut [Scalar::<Bn256>]>(a) };
    let omega = unsafe{std::mem::transmute::<&F, &bn256::Fr>(omega)};
    best_fft_gpu(&mut fft_kern, a, &worker, omega,log_n).unwrap();
    drop(fft_kern);
    info!("end best_fft_bn256_gpu");
}

// #### some unkown bugs
// pub(crate) fn best_fft_bn256_gpu<F: PrimeField>(
//     a: &mut [F],
//     worker: &Worker,
//     omega: &F,
//     log_n: u32,
//     use_cpus_hint: Option<usize>
// )
// {
//     if let Some(ref mut kern) = kern {
//
//         let b = unsafe { std::mem::transmute::<&mut [F], &mut [Scalar::<Bn256>]>(a) };
//         let mut c: Vec<Scalar<Bn256>> = b.into_iter().map(|e| (*e)).collect();
//         let omega = unsafe{std::mem::transmute::<&F, &bn256::Fr>(omega)};
//
//         // if kern
//         //     .with(|k: &mut gpu::FFTKernel<Bn256>| gpu_fft(k, b, omega, log_n))
//         //     .is_ok()
//         // {
//         //     return Ok(());
//         // }
//         // else {
//         //    return Err(GPUError::Simple("GPU accelerator does not correctly execute"));
//         // }
//
//         if kern
//             .with(|k: &mut gpu::FFTKernel<Bn256>| gpu_fft(k, &mut c, omega, log_n))
//             .is_ok()
//         {
//             return Ok(());
//         }
//         else {
//             return Err(GPUError::Simple("GPU accelerator does not correctly execute"));
//         }
//
//     }
//     best_fft(a, worker, omega, log_n, use_cpus_hint);
//     return Ok(())
// }


// important
pub(crate) fn best_fft<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, use_cpus_hint: Option<usize>)
{
    let log_cpus = if let Some(hint) = use_cpus_hint {
        assert!(hint <= worker.cpus);
        let hint = if hint > 0 {
            log2_floor(hint)
        } else {
            0
        };

        hint
    } else {
        worker.log_num_cpus()
    };

    if log_cpus == 0 || log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus);
    }
}



pub(crate) fn serial_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32)
{
    #[inline(always)]
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(n / (2*m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = F::one();
            for j in 0..m {
                let mut t = a[(k+j+m) as usize];
                t.mul_assign(&w);
                let mut tmp = a[(k+j) as usize];
                tmp.sub_assign(&t);
                a[(k+j+m) as usize] = tmp;
                a[(k+j) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2*m;
        }
        
        m *= 2;
    }
}

pub(crate) fn parallel_fft<F: PrimeField>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    log_cpus: u32
)
{
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![F::zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    worker.scope(0, |scope, _| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move |_| {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow(&[j as u64]);
                let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                let mut elt = F::one();
                for i in 0..(1 << log_new_n) {
                    for s in 0..num_cpus {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];

                        t.mul_assign(&elt);
                        tmp[i].add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n);
            });
        }
    });

    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                let mask = (1 << log_cpus) - 1;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_cpus];
                    idx += 1;
                }
            });
        }
    });
}