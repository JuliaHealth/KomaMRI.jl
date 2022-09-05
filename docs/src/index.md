# KomaMRI

KomaMRI.jl is a Julia package meant to simulate general Magnetic Resonance Imaging (MRI) scenarios. Its name comes from the Japanese word for spinning-top こま (ko-ma) as they precess due to gravity like spins in a magnetic field.

```@raw html
<p align="center">
<img width="50%" src="assets/logo-dark.svg"/>
</p>
```

A user willing to interact with KomaMRI can picture this software as the interaction of the following 3 components:
* **Simulator**: generates the **raw signal** by solving the **bloch equation** with the **scanner**, **phantom** and **sequence** inputs.
* **Reconstructor**: reconstructs the **image** by using the **MRIReco.jl** package from the **raw signal**.
* **User Interface**: encapsulates the functionalities of the **Simulator** and the **Reconstructor** into an **user-friendly** interface.

```@raw html
<p align="center">
<img width="100%" src="assets/koma-schema.svg"/>
</p>
```

A user can use the KomaMRI package by using:
* **User Interface**: user-friendly interaction. No Julia programming skills required. Refer to [UI Example](getting-started.md#UI-Example) to check the simplest example.
* **Julia REPL or Scripts** : command line interface interaction. Basic knowledge of Julia is required. Refer to [1-Spin Example](simulation-examples.md#Spin-Example) to check the simplest example. Refer to the [API documentation](api.md) to discover all the possibilities the KomaMRI package offers.

