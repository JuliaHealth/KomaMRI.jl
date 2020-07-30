![equation](http://latex.codecogs.com/gif.latex?O_t%3D%5Ctext%20%7B%20Onset%20event%20at%20time%20bin%20%7D%20t)
![equation](http://latex.codecogs.com/gif.latex?s%3D%5Ctext%20%7B%20sensor%20reading%20%7D) 
![equation](http://latex.codecogs.com/gif.latex?P%28s%20%7C%20O_t%20%29%3D%5Ctext%20%7B%20Probability%20of%20a%20sensor%20reading%20value%20when%20sleep%20onset%20is%20observed%20at%20a%20time%20bin%20%7D%20t)

# MRIsim.jl

[![Build Status](https://travis-ci.com/cncastillo/MRIsim.jl.svg?branch=master)](https://travis-ci.com/cncastillo/MRIsim.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cncastillo/MRIsim.jl?svg=true)](https://ci.appveyor.com/project/cncastillo/MRIsim-jl)
[![Codecov](https://codecov.io/gh/cncastillo/MRIsim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cncastillo/MRIsim.jl)
[![Coveralls](https://coveralls.io/repos/github/cncastillo/MRIsim.jl/badge.svg?branch=master)](https://coveralls.io/github/cncastillo/MRIsim.jl?branch=master)

MRIsim.jl is a Julia package to simulate Magnetic Resonance Imaging (MRI) acquisitions. The main focus of this package is to simulate general scenarios that could arise in pulse sequence development.

**TO-DO**:

 - [ ] RF excitation (under development),
 - [ ] Magnetic susceptibility,
 - [ ] Signal "Raw Output" dictionary (ISMRMRD),
 - [ ] MRIReco.jl for the reconstruciton,
 - [ ] Scanner object: $B_0,\,B_1,\,G_{\max},\,S_{\max}$,etc.
 - [ ] Coil sensitivities.


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
