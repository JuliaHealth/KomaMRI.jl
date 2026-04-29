_addblock_term!(seq::Sequence, chunk::Sequence) = append!(seq, chunk)
_addblock_term!(seq::Sequence, events::_BlockEventTuple) = addblock!(seq, events)
_addblock_term!(seq::Sequence, events::_BlockEvent...) = addblock!(seq, events...)
_addblock_term!(seq::Sequence, event) = addblock!(seq, event)

_addblock_fresh_term!(seq::Sequence, chunk::Sequence) = _append_owned!(seq, chunk)
_addblock_fresh_term!(seq::Sequence, event) = _addblock_term!(seq, event)

function _addblock_transformed_term!(seq::Sequence, op, events...; x=nothing, y=nothing, z=nothing)
    return _append_owned!(seq, op * _block_sequence(events; x, y, z))
end

function _addblock_check!(seq::Sequence, ctx)
    if get(ctx, :check_timing, false)
        sys = get(ctx, :sys, nothing)
        isnothing(sys) ? check_timing(seq) : check_timing(seq, sys)
    end
    if get(ctx, :check_hw_limits, false)
        sys = get(ctx, :sys, nothing)
        isnothing(sys) ? check_hw_limits(seq) : check_hw_limits(seq, sys)
    end
    return seq
end

function _addblock_terms(ex)
    if ex isa Expr && ex.head == :call && ex.args[1] == :+
        return mapreduce(_addblock_terms, append!, ex.args[2:end]; init=Any[])
    end
    return Any[ex]
end

function _has_addblock_tuple_kwargs(ex)
    ex isa Expr || return false
    ex.head == :(=) && ex.args[1] in (:x, :y, :z) && return true
    ex.head == :tuple && any(arg -> arg isa Expr && arg.head in (:(=), :parameters), ex.args) && return true
    _addblock_semicolon_splat(ex) !== nothing && return true
    return any(_has_addblock_tuple_kwargs, ex.args)
end

function _addblock_semicolon_splat(ex)
    ex isa Expr && ex.head == :block || return nothing
    args = Any[arg for arg in ex.args if !(arg isa LineNumberNode)]
    length(args) == 2 || return nothing
    args[2] isa Expr && args[2].head == :... || return nothing
    return Expr(:tuple, Expr(:parameters, args[2]), args[1])
end

function _addblock_keyword(arg)
    if arg isa Expr && arg.head == :(=)
        key = arg.args[1]
        key in (:x, :y, :z) || error("@addblock only accepts `x=`, `y=`, or `z=` in block tuples.")
        return Expr(:kw, key, arg.args[2])
    elseif arg isa Expr && arg.head == :kw
        key = arg.args[1]
        key in (:x, :y, :z) || error("@addblock only accepts `x`, `y`, or `z` keyword axes in block tuples.")
        return arg
    end
    return nothing
end

_addblock_parameter(arg) = (arg isa Expr && arg.head == :...) ? arg : _addblock_keyword(arg)

function _addblock_tuple(term)
    term isa Expr && term.head == :(=) && term.args[1] in (:x, :y, :z) && return Expr(:tuple, term)
    return term
end

function _addblock_call(seq, term)
    addblock = GlobalRef(@__MODULE__, :addblock!)
    addterm = GlobalRef(@__MODULE__, :_addblock_term!)
    addfresh = GlobalRef(@__MODULE__, :_addblock_fresh_term!)
    term = _addblock_tuple(term)
    term = something(_addblock_semicolon_splat(term), term)
    transformed = _addblock_transformed_call(seq, term)
    isnothing(transformed) || return transformed
    if term isa Expr && term.head == :tuple
        return _addblock_tuple_call(addblock, seq, term)
    end
    if term isa Expr && term.head == :call && term.args[1] in (:+, :-, :*, :/)
        return Expr(:call, addfresh, seq, term)
    end
    return Expr(:call, addterm, seq, term)
end

function _addblock_transformed_call(seq, term)
    term isa Expr && term.head == :call && term.args[1] == :* && length(term.args) == 3 || return nothing
    rhs = _addblock_tuple(term.args[3])
    rhs = something(_addblock_semicolon_splat(rhs), rhs)
    rhs isa Expr && rhs.head == :tuple || return nothing
    return _addblock_tuple_call(GlobalRef(@__MODULE__, :_addblock_transformed_term!), seq, rhs, (term.args[2],))
end

function _addblock_tuple_call(f, seq, term, prefix=())
    events = Any[]
    kws = Any[]
    for arg in term.args
        kw = _addblock_keyword(arg)
        if !isnothing(kw)
            push!(kws, kw)
        elseif arg isa Expr && arg.head == :parameters
            for kwarg in arg.args
                kw = _addblock_parameter(kwarg)
                isnothing(kw) && error("@addblock only accepts `x`, `y`, or `z` keyword axes in block tuples.")
                push!(kws, kw)
            end
        else
            push!(events, arg)
        end
    end
    args = isempty(kws) ? Any[f, seq, prefix..., events...] : Any[f, Expr(:parameters, kws...), seq, prefix..., events...]
    return Expr(:call, args...)
end

function _addblock_assignment(lhs, rhs, ctx)
    lhs_value = gensym(:lhs)
    seq_type = GlobalRef(@__MODULE__, :Sequence)
    if isnothing(ctx)
        seq_branch = Expr(:block, (_addblock_call(lhs_value, term) for term in _addblock_terms(rhs))..., lhs_value)
    else
        incoming = gensym(:incoming)
        seq_branch = Expr(:block,
            Expr(:(=), incoming, Expr(:call, seq_type)),
            :($incoming.DEF = $(GlobalRef(Base, :deepcopy))($lhs_value.DEF)),
            (_addblock_call(incoming, term) for term in _addblock_terms(rhs))...,
            Expr(:call, GlobalRef(@__MODULE__, :_addblock_check!), incoming, ctx),
            Expr(:call, GlobalRef(Base, :append!), lhs_value, incoming),
            lhs_value,
        )
    end
    fallback = _has_addblock_tuple_kwargs(rhs) ?
        Expr(:call, GlobalRef(Base, :error), "@addblock tuple keyword syntax requires a Sequence left-hand side.") :
        :($lhs_value + $rhs)
    return :($lhs = let $lhs_value = $lhs
        if $lhs_value isa $seq_type
            $seq_branch
        else
            $fallback
        end
    end)
end

function _rewrite_addblock(ex, target, ctx)
    if ex isa Expr
        if ex.head == :+= && (target === nothing || ex.args[1] == target)
            return _addblock_assignment(ex.args[1], ex.args[2], ctx)
        end
        return Expr(ex.head, (_rewrite_addblock(arg, target, ctx) for arg in ex.args)...)
    end
    return ex
end

function _split_addblock_args(args)
    isempty(args) && error("@addblock expects an expression.")
    ctx = nothing
    i = 1
    kws = Any[]
    while i < length(args) && args[i] isa Expr && args[i].head == :(=) && args[i].args[1] in (:check_timing, :check_hw_limits, :sys)
        push!(kws, Expr(:kw, args[i].args...))
        i += 1
    end
    if !isempty(kws)
        ctx = Expr(:tuple, Expr(:parameters, kws...))
    end
    i <= length(args) || error("@addblock expects an expression after its options.")
    length(args) == i || error("@addblock expects options followed by one expression.")
    return ctx, args[i]
end

"""
    @addblock expr
    @addblock check_timing=true check_hw_limits=true expr

Rewrite one block expression to append in place.

# Examples
```julia
seq = Sequence()
@addblock seq += (rf, z=gz)
@addblock check_timing=true seq += (adc, x=gx)
```
"""
macro addblock(args...)
    ctx, ex = _split_addblock_args(args)
    if ex isa Expr && ex.head == :(=)
        error("@addblock does not create sequences. Use `seq = Sequence()` or `seq = Sequence(sys)`, then `@addblock seq += ...`.")
    end
    return esc(_rewrite_addblock(ex, nothing, ctx))
end

"""
    @addblocks expr
    @addblocks check_timing=true check_hw_limits=true expr

Rewrite block expressions inside a loop or `begin ... end` block.

# Examples
```julia
@addblocks for ky in 1:Ny
    seq += (rf, z=gz)
    seq += readout(ky)
end
```
"""
macro addblocks(args...)
    ctx, ex = _split_addblock_args(args)
    return esc(_rewrite_addblock(ex, nothing, ctx))
end
