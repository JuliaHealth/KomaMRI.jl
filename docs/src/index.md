```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: KomaMRI.jl
  text: Fast and Extensible MRI Simulation in Julia
  tagline: Pulseq in, ISMRMRD out. Fast CPU and GPU simulation with a GUI, interactive PlotlyJS visualizations, and dynamic phantoms.
  image:
    src: /assets/logo.svg
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

  - icon: 🖥️
    title: Interactive GUI
    details: Configure scanners, phantoms, and sequences in KomaUI. Explore results with interactive PlotlyJS visualizations.
    link: /how-to/2-1-use-koma-ui

  - icon: 🌊
    title: Dynamic Motion Models
    details: Simulate static and dynamic phantoms. Model motion and complex spin trajectories with reusable HDF5 phantom files.
    link: /explanation/3-phantom-format
---
```

## Introduction

**KomaMRI** is a Julia package for simulating general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma), as they precess due to gravity like spins in a magnetic field.

**KomaMRI** generates **raw data** by solving the **Bloch equations** using the specified **scanner**, **phantom** and **sequence**. It also provides a Graphical User Interface (GUI) that encapsulates the whole imaging pipeline (simulation and reconstruction).

```@raw html
<p align="center"><img class="docs-light-only" width="100%" src="/assets/koma-schema.svg"/></p>
<p align="center"><img class="docs-dark-only"  width="100%" src="/assets/koma-schema-dark.svg"/></p>
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

KomaMRI supports GPU acceleration with CUDA, AMDGPU, Metal, and oneAPI. To use GPU acceleration, install the corresponding backend package:

:::code-group

```julia [NVIDIA GPUs]
using Pkg
Pkg.add("CUDA")
```

```julia [AMD GPUs]
using Pkg
Pkg.add("AMDGPU")
```

```julia [Apple Silicon]
using Pkg
Pkg.add("Metal")
```

```julia [Intel GPUs]
using Pkg
Pkg.add("oneAPI")
```

:::

We organized the documentation following the philosophy presented by [David Laing](https://documentation.divio.com/).

!!! details "How to Cite Koma"
    If you use Koma, please cite our paper:

    **Plain Text:**

    ```
    Castillo-Passi, C, Coronado, R, Varela-Mattatall, G, Alberola-López, C, Botnar, R, Irarrazaval, P. KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration. Magn Reson Med. 2023; 1- 14. doi: 10.1002/mrm.29635
    ```

    **BibTex:**
    
    ```
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

## Features

* Fast simulations by using CPU and GPU parallelization 🏃💨.
* Open Source, so anyone can include additional features 🆙.
* Compatibility with community-standards 🤝 like Pulseq `.seq` and ISMRMRD `.mrd`.
* Compatibility with [Pluto](how-to/2-2-use-koma-notebooks.md#Pluto) and [Jupyter](how-to/2-2-use-koma-notebooks.md#Jupyter) notebooks 🎈
* Interactive visualizations using PlotlyJS.jl 📲
* Cross-platform 🌐 thanks to the use of the Julia programming language.
* Friendly user interface for people with no programming skills 😌.
* Flexible API for advanced users 👨‍💻.

## Potential Use Cases

* The generation of synthetic data to train Machine Learning models.
* To test novel pulse sequences before implementing them directly in a real scanner (with a Pulseq sequence).
* Teaching exercises for **MRI** acquisition or reconstruction.