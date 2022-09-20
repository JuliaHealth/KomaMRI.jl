# Useful Information

This section covers some topics that are helpful to get a deeper understanding and more insights about a variety of topics

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
i          &: & sequence \: block \: ID \\
RF[i]      &: & RF \: pulse \: at \: the \: i \: block \\
G_j[i]     &: &gradients \: at \: the \: i \: block, \: \forall j \in \{x,y,z\} \\
ADC[i]     &: &acquisition \: at \: the \: i \: block \\
DUR[i]     &: &duration \: at \: the \: i \: block
\end{matrix*}
```

Additionally, there are associated some uniform time resolution parameters or **raster** times for the **RF**, **gradients** and **adquistion**:

```math
\begin{matrix*}[l]
\Delta t_{RF}   &: & raster \: time \: for \: RF \: pulses\\
\Delta t_{G}    &: & raster \: time \: for \: gradients\\
\Delta t_{ADC}  &: & raster \: time \: for \: acquisition
\end{matrix*}
```
