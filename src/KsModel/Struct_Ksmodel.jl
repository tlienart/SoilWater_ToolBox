# =============================================================
#		module: ksModel
# =============================================================
module ksModel

	Base.@kwdef mutable struct KSMODELτ
		τ₁        :: Vector{Float64}
		τ₂        :: Vector{Float64}
		τ₃        :: Vector{Float64}
		τ₄        :: Vector{Float64}
		τ₅        :: Vector{Float64}
		τ₆        :: Vector{Float64}
		τ₁Mac     :: Vector{Float64}
		τ₂Mac     :: Vector{Float64}
		τ₃Mac     :: Vector{Float64}

		τ₁_Min    :: Vector{Float64}
		τ₂_Min    :: Vector{Float64}
		τ₃_Min    :: Vector{Float64}
		τ₄_Min    :: Vector{Float64}
		τ₅_Min    :: Vector{Float64}
		τ₆_Min    :: Vector{Float64}
		
		τ₁Mac_Min :: Vector{Float64}
		τ₂Mac_Min :: Vector{Float64}
		τ₃Mac_Min :: Vector{Float64}
		
		τ₁_Max    :: Vector{Float64}
		τ₂_Max    :: Vector{Float64}
		τ₃_Max    :: Vector{Float64}
		τ₄_Max    :: Vector{Float64}
		τ₅_Max    :: Vector{Float64}
		τ₆_Max    :: Vector{Float64}
		τ₁Mac_Max :: Vector{Float64}
		τ₂Mac_Max :: Vector{Float64}
		τ₃Mac_Max :: Vector{Float64}

		Nse_τ     :: Vector{Float64}
		Rmse_τ    :: Vector{Float64}
		Wilmot_τ  :: Vector{Float64}
	end # mutable struct KSMODEL

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STRUCT_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STRUCT_KSMODEL(;Nτ_Layer = 2::Int64)
			τ₁       = fill(0.0::Float64, Nτ_Layer)
			τ₂       = fill(0.0::Float64, Nτ_Layer)
			τ₃       = fill(0.0::Float64, Nτ_Layer)
			τ₄       = fill(0.0::Float64, Nτ_Layer)
			τ₅       = fill(0.0::Float64, Nτ_Layer)
			τ₆
			τ₁Mac    = fill(0.0::Float64, Nτ_Layer)
			τ₂Mac    = fill(0.0::Float64, Nτ_Layer)
			τ₃Mac    = fill(0.0::Float64, Nτ_Layer)

			τ₁_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₂_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₃_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₄_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₅_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₆_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₁Mac_Min= fill(0.0::Float64, Nτ_Layer)
			τ₂Mac_Min= fill(0.0::Float64, Nτ_Layer)
			τ₃Mac_Min= fill(0.0::Float64, Nτ_Layer)
			
			τ₁_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₂_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₃_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₄_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₅_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₆_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₁Mac_Max= fill(0.0::Float64, Nτ_Layer)
			τ₂Mac_Max= fill(0.0::Float64, Nτ_Layer)
			τ₃Mac_Max= fill(0.0::Float64, Nτ_Layer)

			Nse_τ  = fill(0.0::Float64, Nτ_Layer)
			Rmse_τ = fill(0.0::Float64, Nτ_Layer)
			Wilmot_τ = fill(0.0::Float64, Nτ_Layer)

			ksmodelτ = KSMODELτ(τ₁, τ₂, τ₃, τ₄, τ₅, τ₆, τ₁Mac, τ₂Mac, τ₃Mac, τ₁_Min, τ₂_Min, τ₃_Min, τ₄_Min, τ₅_Min, τ₆_Min, τ₁Mac_Min, τ₂Mac_Min, τ₃Mac_Min, τ₁_Max, τ₂_Max, τ₃_Max, τ₄_Max, τ₅_Max, τ₆_Max, τ₁Mac_Max,τ₂Mac_Max, τ₃Mac_Max, Nse_τ, Rmse_τ, Wilmot_τ)

		return ksmodelτ 
		end  # function: STRUCT_KSMODEL

end  # module: ksModel
# ............................................................