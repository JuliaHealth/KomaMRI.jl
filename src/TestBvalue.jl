using MRIsim
using MRIsim: delay, plot_grads, get_bvalue, dur, DIF_base, Grad_fun, get_qvector
Gmax = 30e-3;
δ = 20e-3;
Δ = 60e-3;
DIF = DIF_base(Gmax,Δ,δ;verbose=true);
plot_grads(DIF)
bvalueR =  get_bvalue(Gmax,Δ,δ) #(2π*γ*G*δ)^2*(Δ-δ/3)
bvalue = get_bvalue(DIF) #My method
p = plot_seq(DIF)
p[1]
## Arbitary waveforms
N = 200
fwf = Grad_fun(x->Gmax,δ,N)
DIF = Sequence([fwf   delay(Δ-δ) -fwf;
                0*fwf delay(Δ-δ) -0*fwf])
plot_grads(DIF)
bvalueR =  get_bvalue(Gmax,Δ,δ)
bvalue =   get_bvalue(DIF)
get_qvector(DIF,type="traj")
