_addblock_term!(seq::Sequence, chunk::Sequence) = append!(seq, chunk)
_addblock_term!(seq::Sequence, events::_BlockEventTuple) = addblock!(seq, events)
_addblock_term!(seq::Sequence, events::_BlockEvent...) = addblock!(seq, events...)
_addblock_term!(seq::Sequence, event) = addblock!(seq, event)

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
    term = _addblock_tuple(term)
    term = something(_addblock_semicolon_splat(term), term)
    if term isa Expr && term.head == :tuple
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
        args = isempty(kws) ? Any[addblock, seq, events...] : Any[addblock, Expr(:parameters, kws...), seq, events...]
        return Expr(:call, args...)
    end
    return Expr(:call, addterm, seq, term)
end

function _addblock_initialization(lhs, rhs)
    seq_type = GlobalRef(@__MODULE__, :Sequence)
    return Expr(:block,
        Expr(:(=), lhs, Expr(:call, seq_type)),
        (_addblock_call(lhs, term) for term in _addblock_terms(rhs))...,
        lhs,
    )
end

function _addblock_assignment(lhs, rhs)
    lhs_value = gensym(:lhs)
    seq_type = GlobalRef(@__MODULE__, :Sequence)
    seq_branch = Expr(:block, (_addblock_call(lhs_value, term) for term in _addblock_terms(rhs))..., lhs_value)
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

function _rewrite_addblock(ex, target)
    if ex isa Expr
        if ex.head == :+= && (target === nothing || ex.args[1] == target)
            return _addblock_assignment(ex.args[1], ex.args[2])
        end
        return Expr(ex.head, (_rewrite_addblock(arg, target) for arg in ex.args)...)
    end
    return ex
end

macro addblock(args...)
    target, ex = if length(args) == 1
        nothing, only(args)
    elseif length(args) == 2
        args[1], args[2]
    else
        error("@addblock expects an expression, or a target sequence and an expression.")
    end
    if target === nothing && ex isa Expr && ex.head == :(=)
        return esc(_addblock_initialization(ex.args[1], ex.args[2]))
    end
    return esc(_rewrite_addblock(ex, target))
end

macro addblocks(ex)
    return esc(_rewrite_addblock(ex, nothing))
end
