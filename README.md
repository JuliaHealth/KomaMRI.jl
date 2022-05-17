> :warning: **Koma.jl v0.6 is coming**: A hugely improved version of the package is going public in the next weeks. Please wait for a little bit, `] add Koma` and `using Koma` will work soon. If you are impatient like me, you could always `] add https://github.com/cncastillo/Koma.jl`.

<p align="center">
<img width="300px" src="./src/ui/assets/Logo.svg#gh-light-mode-only"/>
<img width="300px" src="./src/ui/assets/Logo_dark.svg#gh-dark-mode-only"/>
</p>

![Build status](https://github.com/cncastillo/Koma.jl/actions/workflows/ci.yml/badge.svg)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://cncastillo.github.io/Koma.jl/) [![DOI](https://zenodo.org/badge/252201289.svg)](https://zenodo.org/badge/latestdoi/252201289)


Koma.jl (formerly MRIsim.jl), whose name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field. 

This package is meant to simulate general Magnetic Resonance Imaging (MRI) scenarios that could arise in pulse sequence development. 

**Roadmap**:
 - [x] Phantom and Sequence data-types,
 - [x] Spin preccesion in gradient-only blocks (simulation optimization),
 - [x] GPU accelaration using CUDA.jl,
 - [x] RF excitation,
 - [x] GPU accelaration of RF excitation,
 - [x] Scanner data-type: <img src="https://latex.codecogs.com/gif.latex?B_0,\,B_1,\,G_{\max},\,S_{\max}">, etc.,
 - [x] [Pulseq](https://github.com/imr-framework/pypulseq) IO,
 - [x] Signal "Raw Output" dictionary ([ISMRMRD](https://ismrmrd.github.io/)),
 - [x] [MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/) for the reconstruciton,
 - [ ] Documentation,
 - [ ] <img src="https://latex.codecogs.com/gif.latex?T_{2}^{*}"> decay,
 - [ ] [Auxiliary Pulseq functions](https://github.com/imr-framework/pypulseq/tree/master/pypulseq),
 - [ ] Coil sensitivities,
 - [ ] Diffusion models with Laplacian Eigen Functions,
 - [ ] Magnetic susceptibility,
 - [ ] Use [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/dev/apps.html) to build a ditributable core or app.

## Installation
To install just write the following in the Julia REPL:

```julia
] add Koma
```
## First run
Koma.jl comes with a handy GUI that contains a brain phantom with an EPI sequence. To open it use:

```julia
KomaUI()
```
Press the button that says "Simulate!" to do your first simulation :). Then, a notification emerge telling you that the simulation was successful. In this notification, you can either select to (1) see the Raw signal or (2) to procced with the reconstruction.

## How to cite
If you use this package please acknowledge us by citing (currently working in the paper!):

```bibtex
@software{carlos_castillo_passi_2021_5507370,
  author       = {Castillo-Passi, Carlos and Irarrazaval, Pablo},
  title        = {cncastillo/Koma.jl},
  month        = sep,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {v0.3.8},
  doi          = {10.5281/zenodo.5784237},
  url          = {https://doi.org/10.5281/zenodo.5784237}
}
```

---

## Koma GUI

![Koma](others/GUI.png)
