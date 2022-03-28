# https://github.com/mbauman/TheNaNTrap.jl/blob/main/src/TheNaNTrap.jl
# Matt Bauman

module TheNaNTrap

using Cassette, ReplMaker

Cassette.@context NaNTrap

nonnan(s, x::Number) = (@assert(!isnan(x), string(s, " returned NaN")); x)
nonnan(s, A::AbstractArray) = (foreach(x->nonnan(s, x), A); A)
nonnan(s, x) = x

Cassette.overdub(::NaNTrap, f, args...; kwargs...) = nonnan(string(f, args), Cassette.recurse(NaNTrap(), f, args...; kwargs...))
function parse_to_expr(s)
    p = Meta.parse(s)
    isa(p, Expr) && p.head in (:using, :import) && return p
    quote $Cassette.@overdub($(NaNTrap()), (()->$p)()) end
end
# This should probably be in an init or something...
initrepl(parse_to_expr,
                prompt_text="NaNTrap> ",
                prompt_color = :blue,
                start_key=')',
                mode_name="no_nan_mode")
end