# KomaMRI.jl

*Magnetic Resonance Imaging Simulator*

## Introduction

**KomaMRI.jl** is a Julia package meant to simulate general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field.

KomaMRI has three main components:
* **Simulator**: generates the **raw signal** by solving the **Bloch equations** with the **scanner**, **phantom** and **sequence** inputs.
* **Reconstructor**: reconstructs the **image** by using the **MRIReco.jl** package using the **raw signal**.
* **User Interface**: encapsulates the functionalities of the **Simulator** and the **Reconstructor** into an **user-friendly** interface.

```@raw html
<p align="center">
<img width="100%" src="assets/koma-schema.svg"/>
</p>
```

A user can operate with the KomaMRI package by using:
* **User Interface**: user-friendly interaction. No Julia programming skills required. Refer to [Graphical User Interface](getting-started.md#Graphical-User-Interface) to check the simplest usage example.
* **Julia REPL or Scripts** : Basic knowledge of Julia is required. Refer to [Brain Example](simulation-examples.md#Brain-Example) to check the one example. Refer to the [API documentation](api.md) to discover all the possibilities the KomaMRI package offers.

## Features

The main concrete task of the **KomaMRI.jl** julia package is to simulate the **Bloch equations** to get the **raw signal**. The task to reconstruct the **image** is delegated to the **MRIReco.jl** package. The second concrete task of **KomaMRI.jl** is to encapsulate the simulation and reconstruction functionalities into a user interface.

Some of the most outstanding features of **KomaMRI.jl** are:
* Useful in a wide range of applications, since it solves directly the **Block equations**.
* Fast simulation time by exploiting MRI physics and sequence properties, allowing CPU and GPU parallelization.
* Wide compatibility by using standard file formats `.h5`, `.phantom`, `.scanner`, `.seq` and `.mrd`.
* Friendly user interface to people with no programming skill.
* Flexible Julia API for more advance users with programming skills.
* Easy access to MRI concepts in a hands-on way for education an research purposes.
* Appealing to test novel pulse sequences before implementing them directly in a real scanner.
* Useful to generate synthetic data to train Machine Learning models.
* Open Source, so anyone can add extra features and join to the KomaMRI community.
* Cross-platform thanks to the use of the Julia programing language.

Here is a table with feature comparison against other general MRI simulators:

| Name        | GUI | GPU | Open | Cross-platform |
|:---         |:---:|:---:|:----:|:--------------:|
| KomaMRI     | ✅  | ✅ | ✅   | ✅            |
| JEMRIS      | ✅  | ❌ | ✅   | ❌ Windows    |
| MRISIMUL    | ❌  | ✅ | ❌   | ❌            |
| BlockSolver | ✅  | ✅ | ❌   | ✅            |
| MRiLab      | ✅  | ✅ | ✅   | ❌ MacOS      |
| coreMRI     | ✅  | ✅ | ❌   | ✅            |
