```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: KomaMRI.jl
  text: Fast and Extensible MRI Simulation in Julia
  tagline: Pulseq in, ISMRMRD out. Fast CPU and GPU simulation with a GUI, interactive Plotly visualizations, and dynamic phantoms.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: KomaMRI
  actions:
    - theme: brand
      text: Getting Started
      link: /introduction/1-getting-started
    - theme: alt
      text: View on GitHub
      link: https://github.com/JuliaHealth/KomaMRI.jl
    - theme: alt
      text: API Reference
      link: /reference/1-api

features:
  - icon: 📦
    title: Pulseq and ISMRMRD Community Standards
    details: Build, read, and write Pulseq sequences and export raw data in ISMRMRD. Integrate directly with reconstruction workflows.
    link: /how-to/3-create-your-own-sequence

  - icon: 🚀
    title: Fast, Device Agnostic & Differentiable
    details: Run on CPU and GPU with CUDA, AMDGPU, Metal, and experimental oneAPI support. Differentiable simulations powered by Enzyme and Reactant.
    link: /tutorial/gen-11-ReverseModeRFOptimization

  - icon: 🌊
    title: Dynamic Motion Models
    details: Simulate static and dynamic phantoms. Model motion and complex spin trajectories with reusable HDF5 phantom files.
    link: /explanation/gen-2-motion

  - icon: 🖥️
    title: Interactive GUI
    details: Configure scanners, phantoms, and sequences in KomaUI. Explore results with interactive Plotly visualizations.
    link: /how-to/1-1-use-koma-ui

---
```

## Installation

KomaMRI.jl is a registered Julia package and can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add KomaMRI
```

Or, alternatively, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("KomaMRI")
```

## GPU Support

KomaMRI supports GPU acceleration with CUDA, AMDGPU, Metal, and experimental oneAPI support. To use GPU acceleration, install the corresponding backend package:

:::code-group

```julia [NVIDIA GPUs]
# Install
import Pkg; Pkg.add("CUDA")
# Load
using KomaMRI, CUDA
```

```julia [AMD GPUs]
# Install
import Pkg; Pkg.add("AMDGPU")
# Load
using KomaMRI, AMDGPU
```

```julia [Apple Silicon]
# Install
import Pkg; Pkg.add("Metal")
# Load
using KomaMRI, Metal
```

```julia [Intel GPUs (experimental)]
# Install
import Pkg; Pkg.add("oneAPI")
# Load
using KomaMRI, oneAPI
```

:::

## Citation

If you use KomaMRI.jl in your research, please cite our papers:

**Castillo-Passi, C, Coronado, R, Varela-Mattatall, G, Alberola-López, C, Botnar, R, Irarrazaval, P. KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration. Magn Reson Med. 2023; 90, no. 1: 1-14. doi: 10.1002/mrm.29635**

**Villacorta-Aylagas, P, Castillo-Passi, C, Kierulf, R.A, Menchón-Lara, R.M, Rodríguez-Galván, J.R, Sierra-Pallares, J.B, Irarrazaval, P, Alberola-López, C. Versatile and Highly Efficient MRI Simulation of Arbitrary Motion in KomaMRI. Magn Reson Med. 2026; 95, no. 3: 1791–1803. doi: 10.1002/mrm.70145.**

```bibtex
@article{https://doi.org/10.1002/mrm.29635,
  author = {Castillo-Passi, Carlos and Coronado, Ronal and Varela-Mattatall, Gabriel and Alberola-López, Carlos and Botnar, René and Irarrazaval, Pablo},
  title = {KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration},
  journal = {Magnetic Resonance in Medicine},
  keywords = {Bloch equations, GPU, GUI, Julia, open source, simulation},
  doi = {https://doi.org/10.1002/mrm.29635},
  url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.29635},
  eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.29635}
}

@article{https://doi.org/10.1002/mrm.70145,
  author = {Villacorta-Aylagas, Pablo and Castillo-Passi, Carlos Andrés and Kierulf, Ryan Anders and Menchón-Lara, Rosa María and Rodríguez-Galván, Justino R. and Sierra-Pallares, José Benito and Irarrazaval, Pablo and Alberola-López, Carlos},
  title = {Versatile and Highly Efficient MRI Simulation of Arbitrary Motion in KomaMRI},
  journal = {Magnetic Resonance in Medicine},
  volume = {95},
  number = {3},
  pages = {1791-1803},
  keywords = {motion, MRI simulation, open-source, performance},
  doi = {https://doi.org/10.1002/mrm.70145},
  url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.70145},
  eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.70145},
  year = {2026}
}
```
