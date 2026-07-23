struct IntegrationNodeSamplingRule{M<:SimulationMethod,R} <: SamplingRule
    method::M
    rule::R
end

KomaMRIBase.preserved_samples(rule::IntegrationNodeSamplingRule) =
    KomaMRIBase.preserved_samples(rule.rule)

function integration_node_times(method, step_times, waveforms)
    isempty(waveforms.rf.t) && return similar(step_times, 0)
    eval_times = eltype(step_times)[]
    rf_start, rf_stop = first(waveforms.rf.t), last(waveforms.rf.t)
    nodes = integration_nodes(method)
    for i in 1:(length(step_times) - 1)
        Δt = step_times[i + 1] - step_times[i]
        Δt > 0 || continue
        rf_start <= step_times[i] && step_times[i + 1] <= rf_stop || continue
        for node in nodes
            0 < node < 1 || continue
            push!(eval_times, step_times[i] + eltype(step_times)(node) * Δt)
        end
    end
    return eval_times
end

function KomaMRIBase.additional_sampling_times(rule::IntegrationNodeSamplingRule, event_times, context::BlockSamplingContext)
    step_extra_times = KomaMRIBase.additional_sampling_times(rule.rule, event_times, context)
    step_times = isempty(step_extra_times) ? event_times : KomaMRIBase.merge_sampling_times(event_times, step_extra_times)
    eval_node_times = integration_node_times(rule.method, step_times, context.waveforms)
    isempty(eval_node_times) && return step_extra_times
    return KomaMRIBase.merge_sampling_times(step_extra_times, eval_node_times)
end
