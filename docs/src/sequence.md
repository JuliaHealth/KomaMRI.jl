# Useful Information

This section covers some topics that are helpful to get a deeper understanding and more insights about a variety of topics.

## Sequence Structure

This subsection dives into some details about how a sequence is constructed. Let's introduce the following simple sequence figure to extend the ideas from a visual example to a more general sequence definition:

```@raw html
<p align="center">
<img width="90%" src="../assets/sequence-diagram.svg"/>
</p>
```

A **sequence** can be thought as and ordered concatenation of blocks over time. Every block is composed by an **RF** pulse, the ``(x,y,z)`` **gradients**,  and the **acquisition** of the samples. There is also a time **duration** associated to each block. For short, we are going to refer to these components like so:

```math
\begin{matrix*}[l]
i          &: & \text{sequence block ID} \\
RF[i]      &: & \text{RF pulse at the $i$ block} \\
G_j[i]     &: & \text{gradients at the $i$ block}, \: \forall j \in \{x,y,z\} \\
ADC[i]     &: & \text{acquisition at the $i$ block} \\
DUR[i]     &: & \text{duration at the $i$ block}
\end{matrix*}
```

Additionally, there are associated some uniform time resolution parameters or **raster** times for the **RF**, **gradients** and **adquistion**:

```math
\begin{matrix*}[l]
\Delta t_{RF}   &: & \text{raster time for RF pulses}\\
\Delta t_{G}    &: & \text{raster time for gradients}\\
\Delta t_{ADC}  &: & \text{raster time for acquisition}
\end{matrix*}
```
