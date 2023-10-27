## Introduction

**KomaMRI** is a Julia package meant to simulate general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field.

**KomaMRI** generates **raw data** by solving the **Bloch equations** using the specified **scanner**, **phantom** and **sequence**. It also provides a Graphical User Interface (GUI) that encapsulates the whole imaging pipeline (simulation and reconstruction).

```@raw html
<p align="center"><img class="display-light-only" width="100%" src="assets/koma-schema.svg"/></p>
<p align="center"><img class="display-dark-only"  width="100%" src="assets/koma-schema-dark.svg""/></p>
```

**KomaMRI** can be used by either:
* **User Interface**: User-friendly interaction. No Julia programming skills are required. Refer to [Simulation with User Interface](ui-details.md) to dive into the details of how to use the GUI.

* **Scripts** : Basic knowledge of Julia is required. Refer to [Simulation with Scripts](programming-workflow.md) to follow a basic workflow on how to work with **KomaMRI**.

If it is the first time you want to use **KomaMRI**, we suggest to navigate the ["Getting Started"](getting-started.md) section to install **Julia**, **KomaMRI** and perform your first simulation.

If you want to go more in depth on how to use **KomaMRI**, check the [Simulation with User Interface](ui-details.md) and [Simulation with Scripts](programming-workflow.md) in which you can contrast the simulation workflow when using the **User Interface** or **Scripts**.

## Features

Some of the features of **KomaMRI** are:
* Fast simulations by using CPU and GPU parallelization 🏃💨.
* Open Source, so anyone can include additional features 🆙.
* Compatibility with community-standards 🤝 like Pulseq `.seq` and ISMRMRD `.mrd`.
* Cross-platform 🌐 thanks to the use of the Julia programing language.
* Friendly user interface for people with no programming skills 😌.
* Flexible API for advanced users 👨‍💻.

## Potential Use Cases

We see Koma being used in:
* The generation of synthetic data to train Machine Learning models.
* To test novel pulse sequences before implementing them directly in a real scanner (with a Pulseq sequence).
* Teaching exercises for  MRI acquisition or reconstruction.
