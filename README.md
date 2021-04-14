# MRIsim.jl

<p align="center">
<img width="250px" src="./docs/src/assets/logo.png"/>
</p>

![Build status](https://github.com/cncastillo/MRIsim.jl/actions/workflows/ci.yml/badge.svg)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://cncastillo.github.io/MRIsim.jl/)

MRIsim.jl is a Julia package to simulate Magnetic Resonance Imaging (MRI) acquisitions. The main focus of this package is to simulate general scenarios that could arise in pulse sequence development.

**TO-DO**:
 - [x] GPU accelaration using CUDA,
 - [ ] RF excitation (under development),
 - [ ] Magnetic susceptibility,
 - [ ] Signal "Raw Output" dictionary ([ISMRMRD](https://ismrmrd.github.io/)),
 - [ ] [MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/) for the reconstruciton,
 - [ ] Scanner object: <img src="https://latex.codecogs.com/gif.latex?B_0,\,B_1,\,G_{\max},\,S_{\max}">, etc. 
 - [ ] Coil sensitivities,
 - [ ] Documentation.


**Documentation**: [cncastillo.github.io/MRIsim.jl](https://cncastillo.github.io/MRIsim.jl/build/index.html).

**Presentation**: [MRIsim - Carlos Castillo.pdf](others/MRIsim-CarlosCastillo.pdf)

---

## SpinLab GUI

![SpinLab](others/GUI.png)

---

## Example 1: Brain phantom for different TEs

![Brain phantom with different TEs](others/TEs.png)

## Example 2: dMRI with multiple diffusion directions  

<img src="others/propagator.gif" width="80%">

## Example 3: Moment-compensated diffusion

<img src="others/MomentCompensation.gif" width="80%">
