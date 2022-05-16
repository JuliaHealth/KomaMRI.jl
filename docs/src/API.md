# Test
## Sequence
```@autodocs
Modules = [Koma]
Private = false
Pages   = ["Sequence.jl"]
```
## Grad
```@autodocs
Modules = [Koma]
Private = false
Pages   = ["Grad.jl"]
```
There should be a plot here:
```@example 1
using Koma
grad = Grad(30e-3, 1e-3, .5e-3, 3e-3)
```
```@example 1
seq = Sequence([grad])
p = plot_seq(seq; slider=false, height=300) 
using PlotlyJS
savefig(p,"assets/Grad.html"); nothing
```

```@raw html
<object type="text/html" data="../assets/Grad.html" style="width:100%;height:330px;"></object>
```

```@example 1
grad2 = Grad([30,50,40].*1e-3, 1e-3, .5e-3, 1e-3)
```

```@example 1
seq = Sequence([grad;grad2;;])
p = plot_seq(seq; slider=false, height=300) 
savefig(p,"assets/Grad2.html"); nothing
```

```@raw html
<object type="text/html" data="../assets/Grad2.html" style="width:100%;height:330px;"></object>
```

```@example 1
seq = Sequence([grad grad2 -grad -grad2 grad grad2])
p = plot_seq(seq; height=400) 
savefig(p,"assets/Grad3.html"); nothing
```

```@raw html
<object type="text/html" data="../assets/Grad3.html" style="width:100%;height:430px;"></object>
```

```@example 1
p = plot_M0(seq; height=400) 
savefig(p,"assets/Grad4.html"); nothing
```

```@raw html
<object type="text/html" data="../assets/Grad4.html" style="width:100%;height:430px;"></object>
```

```@example 
using Plots
x = 0:.001:1
y = sin.(10*x)
plot(x,y)
```
## RF
```@autodocs
Modules = [Koma]
Private = false
Pages   = ["RF.jl"]
```
## ADC
```@autodocs
Modules = [Koma]
Private = false
Pages   = ["ADC.jl"]
```
## PulseDesigner
```@autodocs
Modules = [Koma]
Private = false
Pages   = ["PulseDesigner.jl"]
```