# MRIsim.jl

[![Build Status](https://travis-ci.com/cncastillo/MRIsim.jl.svg?branch=master)](https://travis-ci.com/cncastillo/MRIsim.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cncastillo/MRIsim.jl?svg=true)](https://ci.appveyor.com/project/cncastillo/MRIsim-jl)
[![Codecov](https://codecov.io/gh/cncastillo/MRIsim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cncastillo/MRIsim.jl)
[![Coveralls](https://coveralls.io/repos/github/cncastillo/MRIsim.jl/badge.svg?branch=master)](https://coveralls.io/github/cncastillo/MRIsim.jl?branch=master)

MRIsim.jl is a Julia package to simulate Magnetic Resonance Imaging (MRI) acquisitions. The main focus of this package is to simulate general scenarios that could arise in pulse sequence development.

**TO-DO**:

 - [ ] RF excitation (under development),
 - [ ] Magnetic susceptibility,
 - [ ] Signal "Raw Output" dictionary ([ISMRMRD](https://ismrmrd.github.io/)),
 - [ ] [MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/) for the reconstruciton,
 - [ ] Scanner object: <img src="https://latex.codecogs.com/gif.latex?B_0,\,B_1,\,G_{\max},\,S_{\max}">, etc. 
 - [ ] Coil sensitivities,
 - Documentation.


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
