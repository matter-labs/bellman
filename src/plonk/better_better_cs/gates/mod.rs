use super::*;
use super::cs::*;
use crate::smallvec::SmallVec;

use crate::pairing::ff::{Field, PrimeField, PrimeFieldRepr};
use crate::pairing::{Engine, CurveAffine, CurveProjective};
use crate::bit_vec::BitVec;

use crate::{SynthesisError};
use std::marker::PhantomData;

use crate::worker::Worker;
use crate::plonk::domains::*;
use crate::plonk::polynomials::*;

use crate::plonk::cs::variable::*;
use crate::plonk::better_cs::utils::*;
use crate::plonk::fft::cooley_tukey_ntt::*;

pub mod selector_optimized_with_d_next;
pub mod main_gate_with_d_next;
pub mod naive_main_gate;