> [!WARNING]
> This crate has been archived. Development has moved to the [zksync-crypto](https://github.com/matter-labs/zksync-crypto) repository.
> Please use it instead.

# bellman "Community edition"
 
Originally developed for ZCash, it has diverged now and focuses solely on the [PLONK](https://eprint.iacr.org/2019/953) proof system. Uses our "community edition" pairing for Ethereum's BN256 curve. 

## Features

Allows one to design PLONK circuits with custom gates and lookup tables with junction with [franklin-crypto](https://github.com/matter-labs/franklin-crypto) gadget library. At the moment the lookup argument implies using the lookup over the first three state columns (usually refered as A/B/C) and allows to have simultaneously a gate and a lookup applied on the same row of the trace.

## License

Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Code Examples:

- [Edcon2019_material](https://github.com/matter-labs/Edcon2019_material)
- [EDCON Workshop record (youtube): Intro to bellman: Practical zkSNARKs constructing for Ethereum](https://www.youtube.com/watch?v=tUY0YGTpehg&t=74s)

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
