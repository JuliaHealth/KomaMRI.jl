<p align="center">
<img width="300px" src="./assets/logo.svg#gh-light-mode-only" title="Ko-ma (こま) is the Japanese word for spinning-top. They precess due to gravity like spins in a magnetic field."/>
<img width="300px" src="./assets/logo-dark.svg#gh-dark-mode-only" title="Ko-ma (こま) is the Japanese word for spinning-top. They precess due to gravity like spins in a magnetic field."/>
</p>

<p align="center">
<a href="https://github.com/cncastillo/KomaMRI.jl/blob/master/LICENSE"><img src="http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat" /></a>
<a href="https://github.com/invenia/BlueStyle"><img src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" /></a>
<a href="https://pkgs.genieframework.com?packages=KomaMRI"><img src="https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/KomaMRI" /></a>
</p>
 
KomaMRI.jl is a Julia pacakge for highly-efficient ⚡ MRI simulations. KomaMRI was built from the ground up to be: easy-to-use, extensible, cross-platform, and powered by open-source community standards. 

<table>
<tr><td rowspan="2">
<img width="350px" src="./docs/src/assets/ui-simulation.gif"/>
</td>
 <td><b>Features:</b></td>
</tr>
<tr>
<td>
 
- Fast simulations with CPU/GPU parallelization 🏃💨
- Extensible, so anyone can include new features 🆙
- Supports community-standards [🤝](## "Pulseq and ISMRMRD") 
- Interactive visualizations using PlotlyJS.jl 📲
- Cross-platform (Windows, Mac & Linux) 🌐
- Friendly GUI (using web technologies) 😌
- Compatible with modern notebooks [🎈](## "Pluto and Jupyter") 
- Flexible API for advanced users 👨‍💻
</td>
</tr>
<tr>
<td><b>
 
Documentation: <br>[![][docr-img]][docr-url] [![][docd-img]][docd-url]</b></td>
<td><b>
 
Paper: <br>[![][paper-img1]][paper-url1]
</b></td>
</tr>
<tr>
<td colspan="2">

| **Packages**                                | **Julia Compat ≥** | **Stable Version**                           | **Build Status**                       | **Code Coverage**                |
|:------------------------------------------- |:-------------------|:---------------------------------------------|:---------------------------------------|:---------------------------------|
| 📦 [KomaMRI.jl](## "User Interface")        | ![][julia-19]      | [![][komamri-version]][komamri-juliahub]     | [![][gh-actions-img1]][gh-actions-url] | [![][codecov-img1]][codecov-url] |
| └ 📦 [KomaMRIBase.jl](## "Custom Types")    | ![][julia-19]      | [![][komabase-version]][komabase-juliahub]   | [![][gh-actions-img5]][gh-actions-url] | [![][codecov-img5]][codecov-url] |
| └ 📦 [KomaMRICore.jl](## "Simulation")      | ![][julia-19]      | [![][komacore-version]][komacore-juliahub]   | [![][gh-actions-img2]][gh-actions-url] | [![][codecov-img2]][codecov-url] |
| └ 📦 [KomaMRIFiles.jl](## "Input/Output")   | ![][julia-19]      | [![][komafiles-version]][komafiles-juliahub] | [![][gh-actions-img4]][gh-actions-url] | [![][codecov-img4]][codecov-url] |
| └ 📦 [KomaMRIPlots.jl](## "Plots")          | ![][julia-19]      | [![][komaplots-version]][komaplots-juliahub] | [![][gh-actions-img3]][gh-actions-url] | [![][codecov-img3]][codecov-url] |
</td>
</tr>
</table>

[julia-19]: https://img.shields.io/badge/julia-v1.9-9558B2?logo=julia

[komamri-version]: https://juliahub.com/docs/General/KomaMRI/stable/version.svg
[komabase-version]: https://juliahub.com/docs/General/KomaMRIBase/stable/version.svg
[komacore-version]: https://juliahub.com/docs/General/KomaMRICore/stable/version.svg
[komaplots-version]: https://juliahub.com/docs/General/KomaMRIPlots/stable/version.svg
[komafiles-version]: https://juliahub.com/docs/General/KomaMRIFiles/stable/version.svg
[komamri-juliahub]: https://juliahub.com/ui/Packages/General/KomaMRI
[komabase-juliahub]: https://juliahub.com/ui/Packages/General/KomaMRIBase
[komacore-juliahub]: https://juliahub.com/ui/Packages/General/KomaMRICore
[komaplots-juliahub]: https://juliahub.com/ui/Packages/General/KomaMRIPlots
[komafiles-juliahub]: https://juliahub.com/ui/Packages/General/KomaMRIFiles

[docr-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docr-url]: https://juliahealth.github.io/KomaMRI.jl/stable/

[docd-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docd-url]: https://juliahealth.github.io/KomaMRI.jl/dev/

[gh-actions-nightly-img]: https://github.com/cncastillo/KomaMRI.jl/workflows/Nightly/badge.svg
[gh-actions-nightly-url]: https://github.com/cncastillo/KomaMRI.jl/actions/workflows/nightly.yml
[gh-actions-img]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img1]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img2]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img3]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img4]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-img5]: https://github.com/cncastillo/KomaMRI.jl/workflows/CI/badge.svg
[gh-actions-url]: https://github.com/cncastillo/KomaMRI.jl/actions

[codecov-img]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg
[codecov-img1]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=komamri
[codecov-img2]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=core
[codecov-img3]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=plots
[codecov-img4]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=files
[codecov-img5]: https://codecov.io/gh/cncastillo/KomaMRI.jl/branch/master/graph/badge.svg?flag=base
[codecov-url]: https://codecov.io/gh/cncastillo/KomaMRI.jl

[arXiv-img1]: https://img.shields.io/badge/arXiv-2107.11000-blue.svg
[arXiv-url1]: https://arxiv.org/abs/2301.02702

[paper-img1]: https://img.shields.io/badge/doi-10.1002/mrm.29635-blue.svg
[paper-url1]: https://doi.org/10.1002/mrm.29635

## News

> 🟢 **[7 Dec 2023] MRI Together:** Koma was present in [MRI Together](https://mritogether.esmrmb.org/) 😼. The talk is available [here](https://www.youtube.com/watch?v=9mRQH8um4-A). Also, I uploaded the promised [educational example](https://cncastillo.github.io/KomaMRI.jl/dev/educational-1d-simulation/).

> 🟢 **[17 Nov 2023] ISMRM Pulseq's virtual meeting:** Pretty excited of being part of [this meeting](https://github.com/pulseq/ISMRM-Virtual-Meeting--November-15-17-2023). The slides can be seen [here](https://github.com/pulseq/ISMRM-Virtual-Meeting--November-15-17-2023/blob/35a8da7eaa0bf42f2127e1338a440ccd4e3ef53c/slides/day3_KomaMRI_simulator_Quantitative_MRI.pdf).

> 🟢 **[27 Jul 2023] JuliaCon 2023 presentation:** I gave a talk at the MIT 😄! A video of the presentation can be seen [here](https://www.youtube.com/watch?v=WVT9wJegC6Q).

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
 - [ ] [Auxiliary Pulseq functions](https://github.com/imr-framework/pypulseq/tree/master/pypulseq),
 - [ ] Coil sensitivities,
 - [ ] Cardiac phantoms, and triggers.
 - [ ] <img src="https://latex.codecogs.com/gif.latex?T_{2}^{*}"> decay,
 
 Next:
 - [ ] Diffusion models with Laplacian Eigen Functions,
 - [ ] Magnetic susceptibility,
 - [ ] Use [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/dev/apps.html) to build a ditributable core or app.
 
</details>


## Installation
To install just **type** the following in the Julia REPL:

```repl
] add KomaMRI
```
For more information about installation insturctions refer the section [Getting Started](https://cncastillo.github.io/KomaMRI.jl/stable/getting-started/) of the documentation.
## First run
KomaMRI.jl features a convenient GUI with predefined simulation inputs (i.e. `Sequence`, `Phantom`, and `Scanner`). To launch the GUI, use the following command:

```julia
using KomaMRI
KomaUI()
```
Press the button that says "Simulate!" to do your first simulation :). Then, a notification will emerge telling you that the simulation was successful. In this notification, you can either select to (1) see the Raw Data or (2) to procced with the reconstruction.

## How to Contribute
KomaMRI exists thanks to all our contributors:

<a href="https://github.com/cncastillo/KomaMRI.jl/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=cncastillo/KomaMRI.jl" height="40px"/>
</a>

Want to be highlighted here? We welcome contributions from the community! If you're interested in contributing, please read our [Contribution Guidelines](CONTRIBUTING.md) for details on how to get started.


## How to Cite
If you use this package please acknowledge us by citing our paper. 

**Plain Text:**
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

## Tested compatibility

We automatically test KomaMRICore's CPU multi-threading support on the most popular operating systems. Nevertheless, for GPU support the process is more manual (until [#147](https://github.com/cncastillo/KomaMRI.jl/issues/147)). Here is a summary of our automatic CPU tests and local GPU tests for multiple versions of Julia:

| Julia (OS)                            | CPU                                       | GPU (Nvidia) |
|:--------------------------------------|:----------------------------------------------:|:------------------:|
| Julia 1.9 (Windows)  | [![][gh-actions-img]][gh-actions-url]          | ✅                |
| Julia 1.9 (Linux)    | [![][gh-actions-img]][gh-actions-url]          | ✅                |
| Julia 1.9  (Mac)      | [![][gh-actions-img]][gh-actions-url]          | ➖                |
| Julia 1.10 (Windows)           | [![][gh-actions-img]][gh-actions-url]          | ✅                |
| Julia 1.10 (Linux)             | [![][gh-actions-img]][gh-actions-url]          | ✅                |
| Julia 1.10 (Mac)               | [![][gh-actions-img]][gh-actions-url]          | ➖                |
| Julia 1.11 (Windows)          | [![][gh-actions-nightly-img]][gh-actions-nightly-url]  | ❌                 |
| Julia 1.11 (Linux)            | [![][gh-actions-nightly-img]][gh-actions-nightly-url]  | ❌                 |
| Julia 1.11 (Mac)              | [![][gh-actions-nightly-img]][gh-actions-nightly-url]  | ➖                 |

If you see any problem with this information, please let us know in the form of a GitHub issue.
