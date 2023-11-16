## Introduction

**KomaMRI** is a Julia package meant to simulate general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field.

**KomaMRI** generates **raw data** by solving the **Bloch equations** using the specified **scanner**, **phantom** and **sequence**. It also provides a Graphical User Interface (GUI) that encapsulates the whole imaging pipeline (simulation and reconstruction).

```@raw html
<p align="center"><img class="display-light-only" width="100%" src="assets/koma-schema.svg"/></p>
<p align="center"><img class="display-dark-only"  width="100%" src="assets/koma-schema-dark.svg""/></p>
```

**KomaMRI** can be used in different environments:

* **User Interface**: User-friendly interaction. No Julia programming skills are required. Refer to the [User Interface Section](ui-details.md) to dive into the details of how to use the GUI.

* **Scripts** : Basic knowledge of **Julia** is required. Refer to the [Scripts Section](programming-workflow.md) to follow a basic workflow on how to work with **KomaMRI**.

* **Notebooks**: Basic knowledge of **Julia** is required. This serves as an alternative development environment featuring user-friendly interactive tools. For guidance on setting up these environments, refer to the [Notebooks Section](#Notebooks).

If you are new to **KomaMRI**, we recommend starting with the ["Getting Started"](getting-started.md) section to install **Julia**, **KomaMRI**, and perform your first simulation.


## Features

Some of the features of **KomaMRI** are:
* Fast simulations by using CPU and GPU parallelization 🏃💨.
* Open Source, so anyone can include additional features 🆙.
* Compatibility with community-standards 🤝 like Pulseq `.seq` and ISMRMRD `.mrd`.
* Compatibility with [Pluto](notebooks.md#Using-KomaMRI-with-Pluto) and [Jupyter](notebooks.md#Using-KomaMRI-with-Jupyter) notebooks 🎈
* Cross-platform 🌐 thanks to the use of the Julia programming language.
* Friendly user interface for people with no programming skills 😌.
* Flexible API for advanced users 👨‍💻.

## Potential Use Cases

We see Koma being used in:
* The generation of synthetic data to train Machine Learning models.
* To test novel pulse sequences before implementing them directly in a real scanner (with a Pulseq sequence).
* Teaching exercises for  MRI acquisition or reconstruction.
