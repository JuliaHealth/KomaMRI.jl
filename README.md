# MRIsim.jl

[![Build Status](https://travis-ci.com/cncastillo/MRIsim.jl.svg?branch=master)](https://travis-ci.com/cncastillo/MRIsim.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/cncastillo/MRIsim.jl?svg=true)](https://ci.appveyor.com/project/cncastillo/MRIsim-jl)
[![Codecov](https://codecov.io/gh/cncastillo/MRIsim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cncastillo/MRIsim.jl)
[![Coveralls](https://coveralls.io/repos/github/cncastillo/MRIsim.jl/badge.svg?branch=master)](https://coveralls.io/github/cncastillo/MRIsim.jl?branch=master)

MRIsim.jl is a Julia package to simulate Magnetic Resonance Imaging (MRI) acquisitions. The main focus of this package is to simulate general scenarios that could arise in pulse sequence development.

**TO-DO**:
 * Coil sensitivities,
 * RF excitation (under development),
 * Magnetic susceptibility,
 * etc.

**Documentation**: [cncastillo.github.io/MRIsim.jl](https://cncastillo.github.io/MRIsim.jl/build/index.html).

**Presentation**: [MRIsim - Carlos Castillo.pdf](others/MRIsim-CarlosCastillo.pdf)

---

## SpinLab GUI
<center>
![SpinLab](others/GUI.png)
</center>
## Example 1: Brain phantom for different TEs
<center>
![Brain phantom with different TEs](others/TEs.png)
</center>
## Example 2: dMRI with multiple diffusion directions  
<center>
<img src="others/propagator.gif" width="80%">
</center>
## Example 3: Moment-compensated diffusion
<center>
<img src="others/MomentCompensation.gif" width="80%">
</center>
