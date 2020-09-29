### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 95164fa0-f7da-11ea-20b3-a5dc6c5e8376
begin
	cd("C:/Users/carlo/Documents/MRIsim.jl")
	import Pkg; Pkg.activate(".")
	using MRIsim
end

# ╔═╡ 5c4f5430-f7dc-11ea-0c16-a92eb2bf2614
using MRIsim: DIF_base, plot_grads, get_bvalue

# ╔═╡ 90e6e480-f7df-11ea-345e-934aa0ac621e
using PlutoUI

# ╔═╡ 428d5360-f7de-11ea-05db-0d216406e737
Gmax = 30e-3;

# ╔═╡ 428f2820-f7de-11ea-2686-ed6820e2ae2c
md"δ = $(@bind δ Slider(1:60))"

# ╔═╡ 429346d0-f7de-11ea-0dc5-b73c3d966ac5
Δ = 60e-3;

# ╔═╡ 429828d0-f7de-11ea-26ff-45b99356ece4
DIF = DIF_base(Gmax,Δ,δ*1e-3;verbose=true);

# ╔═╡ cd40a6d0-f7e1-11ea-39d3-038bbf0e03b3
bvalue = get_bvalue(DIF)*1e-6 #My method

# ╔═╡ 730fb0f0-f7de-11ea-226d-e384481c8f27
plot_grads(DIF)

# ╔═╡ 4a5b4fc0-f7de-11ea-0a76-3b1a0051203b
begin
	idx = ["Gx" "Gy" "Gz"]
	M, N = size(DIF.GR)
	G = [DIF.GR[j,floor(Int,i/2)+1].A for i=0:2*N-1, j=1:M]
	T = [DIF.GR[1,i].T for i=1:N]
	t = [sum(T[1:i]) for i=1:N]
	t = [t[floor(Int,i/2)+1] for i=0:2*N-1]
	t = [0; t[1:end-1]]
	p = [plot() for j=1:M]
	for j=1:size(DIF.GR,1)
		p[j] = plot(t*1e3, G[:,j]*1e3,label=idx[j],line_shape="hv")
		plot!(size=(600,300))
	end
end

# ╔═╡ 4cc5074e-f7df-11ea-28ba-431a64630a6d
p[1]

# ╔═╡ 5dc71c90-f7db-11ea-1c0c-17af91306350


# ╔═╡ c9952f32-f7da-11ea-3a1f-53bb0675fa7f


# ╔═╡ bb97f250-f7da-11ea-3e07-0da664be5e76


# ╔═╡ b3c40280-f7da-11ea-2fa4-93274f77fe4f


# ╔═╡ Cell order:
# ╠═95164fa0-f7da-11ea-20b3-a5dc6c5e8376
# ╠═5c4f5430-f7dc-11ea-0c16-a92eb2bf2614
# ╠═90e6e480-f7df-11ea-345e-934aa0ac621e
# ╠═428d5360-f7de-11ea-05db-0d216406e737
# ╟─428f2820-f7de-11ea-2686-ed6820e2ae2c
# ╠═429346d0-f7de-11ea-0dc5-b73c3d966ac5
# ╠═429828d0-f7de-11ea-26ff-45b99356ece4
# ╠═cd40a6d0-f7e1-11ea-39d3-038bbf0e03b3
# ╠═730fb0f0-f7de-11ea-226d-e384481c8f27
# ╟─4a5b4fc0-f7de-11ea-0a76-3b1a0051203b
# ╠═4cc5074e-f7df-11ea-28ba-431a64630a6d
# ╟─5dc71c90-f7db-11ea-1c0c-17af91306350
# ╟─c9952f32-f7da-11ea-3a1f-53bb0675fa7f
# ╟─bb97f250-f7da-11ea-3e07-0da664be5e76
# ╟─b3c40280-f7da-11ea-2fa4-93274f77fe4f
