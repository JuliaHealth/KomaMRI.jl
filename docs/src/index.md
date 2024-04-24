## Introduction

**KomaMRI** is a Julia package meant to simulate general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field.

**KomaMRI** generates **raw data** by solving the **Bloch equations** using the specified **scanner**, **phantom** and **sequence**. It also provides a Graphical User Interface (GUI) that encapsulates the whole imaging pipeline (simulation and reconstruction).

```@raw html
<p align="center"><img class="display-light-only" width="100%" src="assets/koma-schema.svg"/></p>
<p align="center"><img class="display-dark-only"  width="100%" src="assets/koma-schema-dark.svg""/></p>
```
We organized the documentation following the philosophy presented by [David Laing](https://documentation.divio.com/).

## Features

Some of the features of **KomaMRI** are:
* Fast simulations by using CPU and GPU parallelization 🏃💨.
* Open Source, so anyone can include additional features 🆙.
* Compatibility with community-standards 🤝 like Pulseq `.seq` and ISMRMRD `.mrd`.
* Compatibility with [Pluto](how-to/2-2-use-koma-notebooks.md#Pluto) and [Jupyter](how-to/2-2-use-koma-notebooks.md#Jupyter) notebooks 🎈
* Cross-platform 🌐 thanks to the use of the Julia programming language.
* Friendly user interface for people with no programming skills 😌.
* Flexible API for advanced users 👨‍💻.

## Potential Use Cases

We see Koma being used in:
* The generation of synthetic data to train Machine Learning models.
* To test novel pulse sequences before implementing them directly in a real scanner (with a Pulseq sequence).
* Teaching exercises for **MRI** acquisition or reconstruction.