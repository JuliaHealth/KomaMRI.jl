## Introduction

**KomaMRI.jl** is a Julia package meant to simulate general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field.

KomaMRI generates **raw data** by solving the **Bloch equations** using the specified **scanner**, **phantom** and **sequence**. It also provides a Graphical User Interface (GUI) that encapsulates the whole imaging pipeline (simulation and reconstruction).

```@raw html
<p align="center"><img class="display-light-only" width="100%" src="assets/koma-schema.svg"/></p>
<p align="center"><img class="display-dark-only"  width="100%" src="assets/koma-schema-dark.svg""/></p>
```

KomaMRI can be used by either:
* **Graphical User Interface**: User-friendly interaction. No Julia programming skills are required. Refer to [Graphical User Interface](getting-started.md#Graphical-User-Interface) to check the simplest usage example.

* **Scripts** : Basic knowledge of Julia is required. Refer to [Brain Example](simulation-examples.md#Brain-Example) to check a tutorial. Refer to the [API documentation](api.md) to discover all the functions that the package has to offer.

## Features

Some of the features of **KomaMRI.jl** are:
* Fast simulations by using CPU and GPU parallelization.
* Open Source, so anyone can include additional features.
* Compatibility with community-standards like Pulseq `.seq` and ISMRMRD `.mrd`.
* Cross-platform thanks to the use of the Julia programing language.
* Friendly user interface for people with no programming skills.
* Flexible API for advance users.

## Potential use cases

We see Koma being used in:
* The generation of synthetic data to train Machine Learning models.
* To test novel pulse sequences before implementing them directly in a real scanner (with a Pulseq sequence).
* Teaching exercises for  MRI acquisition or reconstruction.

