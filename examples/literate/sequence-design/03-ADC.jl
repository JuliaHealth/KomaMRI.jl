# # ADC Design

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./FILE_NAME.jl)
#md # [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./FILE_NAME.ipynb)

# Let's import the KomaMRI package

using KomaMRI

# ## ADC definition

N, T = 8, 7e-3
adc = ADC(N, T)

g0, r0 = Grad(0, 0), RF(0, 0)
g0_m = reshape([g0; g0; g0], :, 1)
r0_m = reshape([r0], :, 1)
adc_v = [adc]
seq = Sequence(g0_m, r0_m, adc_v)

p1 = plot_seq(seq; slider=false, height=300)
#md savefig(p1, "../../assets/FOLDER_NAME/1-adc.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/1-adc.html" style="width:100%; height:320px;"></object>
#md # ```

# ## ADC phase compensation
