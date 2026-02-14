---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: KomaMRI.jl
  text: Fast and Extensible MRI Simulation in Julia
  tagline: Pulseq in, ISMRMRD out. Fast CPU and GPU simulation with a GUI, interactive PlotlyJS visualizations, and dynamic phantoms.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: KomaMRI
  actions:
    - theme: brand
      text: Getting Started
      link: /how-to/1-getting-started
    - theme: alt
      text: View on GitHub
      link: https://github.com/JuliaHealth/KomaMRI.jl
    - theme: alt
      text: API Reference
      link: /reference/1-api

features:
  - icon: 📦
    title: Community Standards
    details: Import Pulseq sequences and export raw data in ISMRMRD. Integrate directly with reconstruction and benchmarking workflows.
    link: /how-to/3-phantom-formats

  - icon: 🚀
    title: Fast, Device Agnostic
    details: Run on CPU and GPU with CUDA, AMDGPU, Metal, and oneAPI support. Differentiable simulations powered by Enzyme and Reactant.
    link: /explanation/6-simulation

  - icon: 🌊
    title: Dynamic Motion Models
    details: Simulate static and dynamic phantoms. Model motion and complex spin trajectories with reusable HDF5 phantom files.
    link: /explanation/3-phantom-format

  - icon: 🖥️
    title: Interactive GUI
    details: Configure scanners, phantoms, and sequences in KomaUI. Explore results with interactive PlotlyJS visualizations.
    link: /how-to/2-1-use-koma-ui

---

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

KomaMRI supports GPU acceleration with CUDA, AMDGPU, Metal, and oneAPI. To use GPU acceleration, install the corresponding backend package:

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

```julia [Intel GPUs]
# Install
import Pkg; Pkg.add("oneAPI")
# Load
using KomaMRI, oneAPI
```

:::

## Citation

If you use KomaMRI.jl in your research, please cite our paper:

**Castillo-Passi, C, Coronado, R, Varela-Mattatall, G, Alberola-López, C, Botnar, R, Irarrazaval, P. KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration. Magn Reson Med. 2023; 1-14. doi: 10.1002/mrm.29635**

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
```
