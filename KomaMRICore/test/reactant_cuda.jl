using Test
using Enzyme
using KomaMRIBase
using KomaMRICore
using Reactant

Reactant.set_default_backend("cuda")
Reactant.allowscalar(false)

include(joinpath(@__DIR__, "test_files", "ad_utils.jl"))

function reactant_cuda_node_loss_and_gradient(x, params)
    result = Enzyme.gradient(
        Enzyme.ReverseWithPrimal,
        bloch_node_ad_loss,
        x,
        Enzyme.Const(params),
    )
    return result.val, result.derivs[1]
end

@info "Reactant CUDA devices" devices=string.(Reactant.devices())

@testset "RF and density controls through simulate and ADC Reactant Enzyme CUDA accuracy" begin
    for sim_method in (BlochSimple(), Bloch())
        params = bloch_node_ad_parameters(sim_method)
        params_ra = bloch_node_ad_reactant_parameters(params)
        x = Reactant.to_rarray(copy(BLOCH_NODE_AD_X0))
        compiled = Reactant.@compile sync=true reactant_cuda_node_loss_and_gradient(x, params_ra)
        loss, gradient = compiled(x, params_ra)

        @test Reactant.to_number(loss) ≈ bloch_node_ad_loss(BLOCH_NODE_AD_X0, params)
        @test Array(gradient) ≈ bloch_node_ad_fd_gradient(params) rtol=1e-8 atol=1e-10
    end
end
