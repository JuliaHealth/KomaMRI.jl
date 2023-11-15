# Create Your Own Sequence

!!! warning
    This section is currently under construction, and some details on how to construct a Sequence may be missing.

This is an example of how to create a **Sequence** struct:

```julia
# Export necessary modules
using KomaMRI

# Create the function that creates a phantom
function sequence_example(FOV::Real, N::Integer)

    # Define initial paramters (TODO: consider when N is even)
    sys = Scanner()
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	Ga ≥ sys.Gmax ? error("Ga=$(Ga*1e3) mT/m exceeds Gmax=$(Gmax*1e3) mT/m, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	ϵ1 = Δτ/(Δτ+ζ)

	# EPI base
	EPI = Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2],  #Gx
	 	[mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2])) #Gy
	EPI.ADC = [mod(i,2)==1 ? ADC(0,Δτ,ζ) : ADC(N,Ta,ζ) for i=0:2*Ny-2]

	# Pre-wind and wind gradients
	ϵ2 = Ta/(Ta+ζ)
    PHASE =   Sequence(reshape(1/2*[Grad(      -Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) # This needs to be calculated differently
	DEPHASE = Sequence(reshape(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) # for even N
	seq = PHASE + EPI + DEPHASE

	# Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")

    # Return the sequence
	return seq
end

# Call the function to create a sequence
FOV, N = 23e-2, 101
seq = sequence_example(FOV, N)

# Plot the sequence in time and its kspace
plot_seq(seq; range=[0 30])
plot_kspace(seq)
```
```@raw html
<object type="text/html" data="../assets/create-your-own-sequence-time.html" style="width:50%; height:420px;"></object><object type="text/html" data="../assets/create-your-own-sequence-kspace.html" style="width:50%; height:420px;"></object>
```