<p align="center">
<img width="300px" src="./src/ui/assets/Logo.svg#gh-light-mode-only"/>
<img width="300px" src="./src/ui/assets/Logo_dark.svg#gh-dark-mode-only"/>
</p>

<div align="center">

##### Main project
| **Documentation**         | **KomaMRI.jl Paper**           | **Build Status**                      |
|:-------------------------:|:------------------------------:|:-------------------------------------:|
| [![][docr-img]][docr-url] | [![][paper-img1]][paper-url1]  | [![][gh-actions-img]][gh-actions-url] |
| [![][docd-img]][docd-url] | [![][arXiv-img1]][arXiv-url1]  | [![][codecov-img]][codecov-url]       |
##### Submodules

| **KomaMRI.jl**            | **KomaMRICore.jl**             | **KomaMRIPlots.jl**                   |
|:-------------------------:|:------------------------------:|:-------------------------------------:|
| [![][gh-actions-img1]][gh-actions-url] | [![][gh-actions-img2]][gh-actions-url] | [![][gh-actions-img3]][gh-actions-url]       |
| [![][codecov-img1]][codecov-url] | [![][codecov-img2]][codecov-url] | [![][codecov-img3]][codecov-url]       |
 ```
📦 KomaMRI.jl (UI)                      
├─ 📦 KomaMRICore.jl (Simulation and IO)
└─ 📦 KomaMRIPlots.jl (Plots)           
```
</div>

[docr-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docr-url]: https://cncastillo.github.io/KomaMRI.jl/stable/

[docd-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docd-url]: https://cncastillo.github.io/KomaMRI.jl/dev/

[gh-actions-img]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img1]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img2]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img3]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-url]: https://github.com/cncastillo/KomaMRI.jl/actions

[codecov-img]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg
[codecov-img1]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=komamri
[codecov-img2]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=core
[codecov-img3]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=plots
[codecov-url]: https://codecov.io/gh/cncastillo/KomaMRI.jl

[arXiv-img1]: https://img.shields.io/badge/arXiv-2107.11000-blue.svg
[arXiv-url1]: https://arxiv.org/abs/2301.02702

[paper-img1]: https://img.shields.io/badge/doi-10.1002/mrm.29635-blue.svg
[paper-url1]: https://doi.org/10.1002/mrm.29635

KomaMRI.jl (formerly MRIsim.jl), whose name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field. 

This package is meant to simulate general Magnetic Resonance Imaging (MRI) scenarios that could arise in pulse sequence development. 

> 🟢 **[29 Jun 2023] [KomaMRI.jl's paper](https://onlinelibrary.wiley.com/doi/10.1002/mrm.29635) was chosen as a July editor's pick in MRM 🥳!**

> 🟢 **[6 Mar 2023] Paper published in MRM 😃!** The open access article is available [here](https://onlinelibrary.wiley.com/doi/10.1002/mrm.29635).

> 🟢 **[8 Dec 2022] KomaMRI v0.7**: Huge code rewrite, this implies improved: performance (now **5x faster**), type stability, extensibility, and more!

> 🟢 **[17 May 2022] ISMRM 2022 presentation**: My abstract presentation is now uploaded [here!](https://www.youtube.com/watch?v=tH_XUnoSJK8)

<details>
<summary> <samp>&#9776; Roadmap</samp></summary>

 v1.0: 
 - [x] Phantom and Sequence data-types,
 - [x] Spin precession in gradient-only blocks (simulation optimization),
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
 - [ ] Cardiac phantoms, and triggers.
 
 Next:
 - [ ] Diffusion models with Laplacian Eigen Functions,
 - [ ] Magnetic susceptibility,
 - [ ] Use [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/dev/apps.html) to build a ditributable core or app.
 
</details>


## Installation
To install just write the following in the Julia REPL:

```repl
] add KomaMRI
```
## First run
KomaMRI.jl comes with a handy GUI that contains a brain phantom with an EPI sequence. To open it use:

```julia
using KomaMRI
KomaUI()
```
Press the button that says "Simulate!" to do your first simulation :). Then, a notification emerge telling you that the simulation was successful. In this notification, you can either select to (1) see the Raw signal or (2) to procced with the reconstruction.

## How to cite
If you use this package please acknowledge us by citing our paper. 

**How to cite:**
> Castillo-Passi, C, Coronado, R, Varela-Mattatall, G, Alberola-López, C, Botnar, R, Irarrazaval, P. KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration. Magn Reson Med. 2023; 1- 14. doi: 10.1002/mrm.29635

**BibTex:**
```bibtex
@article{https://doi.org/10.1002/mrm.29635,
author = {Castillo-Passi, Carlos and Coronado, Ronal and Varela-Mattatall, Gabriel and Alberola-López, Carlos and Botnar, René and Irarrazaval, Pablo},
title = {KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration},
journal = {Magnetic Resonance in Medicine},
keywords = {Bloch equations, GPU, GUI, Julia, open source, simulation},
doi = {https://doi.org/10.1002/mrm.29635},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.29635},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.29635},
}
```
---

## Koma GUI

<img width="100%" src="/others/GUI.svg"/>
