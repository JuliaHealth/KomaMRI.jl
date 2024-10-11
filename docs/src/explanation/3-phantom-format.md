# Phantom File Format

## Introduction

While there is already an open and fairly standardised format for MRI sequences 
such as [Pulseq](https://pulseq.github.io/index.html), this is not the case for digital phantoms.
That's why we defined a new ''.phantom'' format, which relies on the [HDF5 standard](https://www.hdfgroup.org/solutions/hdf5/).
HDF5 is specially designed to store large amounts of heterogeneous data and to make it readable 
and writable quickly and easily. In addition, it allows the storage of metadata. 
For all these reasons, it is the ideal file format for storing phantoms. 

## File Format Specification

### Phantom File Tree

```@raw html
<p><img class="docs-light-only" width="80%" src="../../assets/ph-phantom-file-format-light.svg"/></p>
<p><img class="docs-dark-only"  width="80%" src="../../assets/ph-phantom-file-format-dark.svg"/></p>
```

### Action types

```@raw html
<p><img class="docs-light-only" width="80%" src="../../assets/ph-action-types-light.svg"/></p>
<p><img class="docs-dark-only"  width="80%" src="../../assets/ph-action-types-dark.svg"/></p>
```

### TimeSpan types

```@raw html
<p><img class="docs-light-only" width="80%" src="../../assets/ph-timespan-types-light.svg"/></p>
<p><img class="docs-dark-only"  width="80%" src="../../assets/ph-timespan-types-dark.svg"/></p>
```

### SpinSpan types

```@raw html
<p><img class="docs-light-only" width="80%" src="../../assets/ph-spinspan-types-light.svg"/></p>
<p><img class="docs-dark-only"  width="80%" src="../../assets/ph-spinspan-types-dark.svg"/></p>
```