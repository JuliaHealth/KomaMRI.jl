This folder contains legacy redirects for documentation pages that no longer exist.

These old URLs are still referenced in the paper:

> Villacorta-Aylagas, P, Castillo-Passi, C, Kierulf, R.A, Menchón-Lara, R.M, Rodríguez-Galván, J.R, Sierra-Pallares, J.B, Irarrazaval, P, Alberola-López, C. Versatile and Highly Efficient MRI Simulation of Arbitrary Motion in KomaMRI. Magn Reson Med. 2026; 95, no. 3: 1791–1803. doi: 10.1002/mrm.70145.

The redirects are kept so those links continue to work.
They currently redirect to documentation version `v0.14`.

During docs generation, `docs/make.jl` copies files from this folder into `docs/src/public/` (removing the `legacy/` prefix), so the old paths remain reachable in deployed docs.