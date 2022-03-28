# =============================================================
#		module: ksModel
# =============================================================
module ksModel

	Base.@kwdef mutable struct KSMODELτ
		τ₁Max     :: Vector{Float64}
		τ₁ηMin    :: Vector{Float64}
		τ₂        :: Vector{Float64}
		τ₃        :: Vector{Float64}
		τ₅        :: Vector{Float64}
		τ₄        :: Vector{Float64}
		τ₁Mac     :: Vector{Float64}
		τ₂Mac     :: Vector{Float64}
		τ₃Mac     :: Vector{Float64}

		τ₁Max_Min    :: Vector{Float64}
		τ₂_Min    :: Vector{Float64}
		τ₃_Min    :: Vector{Float64}
		τ₁ηMin_Min    :: Vector{Float64}
		τ₅_Min    :: Vector{Float64}
		τ₄_Min    :: Vector{Float64}
		
		τ₁Mac_Min :: Vector{Float64}
		τ₂Mac_Min :: Vector{Float64}
		τ₃Mac_Min :: Vector{Float64}
		
		τ₁Max_Max    :: Vector{Float64}
		τ₂_Max    :: Vector{Float64}
		τ₃_Max    :: Vector{Float64}
		τ₁ηMin_Max :: Vector{Float64}
		τ₅_Max    :: Vector{Float64}
		τ₄_Max    :: Vector{Float64}
		τ₁Mac_Max :: Vector{Float64}
		τ₂Mac_Max :: Vector{Float64}
		τ₃Mac_Max :: Vector{Float64}

		Nse_τ     :: Vector{Float64}
		Rmse_τ    :: Vector{Float64}
		Wilmot_τ  :: Vector{Float64}
		Ccc_τ  :: Vector{Float64}
	end # mutable struct KSMODEL

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STRUCT_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STRUCT_KSMODEL(;Nτ_Layer = 2::Int64)
			τ₁Max    = fill(0.0::Float64, Nτ_Layer)
			τ₁ηMin   = fill(0.0::Float64, Nτ_Layer)
			τ₂       = fill(0.0::Float64, Nτ_Layer)
			τ₃       = fill(0.0::Float64, Nτ_Layer)
			τ₅       = fill(0.0::Float64, Nτ_Layer)
			τ₄       = fill(0.0::Float64, Nτ_Layer)
			τ₁Mac    = fill(0.0::Float64, Nτ_Layer)
			τ₂Mac    = fill(0.0::Float64, Nτ_Layer)
			τ₃Mac    = fill(0.0::Float64, Nτ_Layer)

			τ₁Max_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₂_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₃_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₁ηMin_Min = fill(0.0::Float64, Nτ_Layer)
			τ₅_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₄_Min   = fill(0.0::Float64, Nτ_Layer)
			τ₁Mac_Min= fill(0.0::Float64, Nτ_Layer)
			τ₂Mac_Min= fill(0.0::Float64, Nτ_Layer)
			τ₃Mac_Min= fill(0.0::Float64, Nτ_Layer)
			
			τ₁Max_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₂_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₃_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₁ηMin_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₅_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₄_Max   = fill(0.0::Float64, Nτ_Layer)
			τ₁Mac_Max= fill(0.0::Float64, Nτ_Layer)
			τ₂Mac_Max= fill(0.0::Float64, Nτ_Layer)
			τ₃Mac_Max= fill(0.0::Float64, Nτ_Layer)

			Nse_τ  = fill(0.0::Float64, Nτ_Layer)
			Rmse_τ = fill(0.0::Float64, Nτ_Layer)
			Wilmot_τ = fill(0.0::Float64, Nτ_Layer)
			Ccc_τ = fill(0.0::Float64, Nτ_Layer)

			ksmodelτ = KSMODELτ(τ₁Max, τ₁ηMin, τ₂, τ₃, τ₅, τ₄, τ₁Mac, τ₂Mac, τ₃Mac, τ₁Max_Min,τ₁ηMin_Min, τ₂_Min, τ₃_Min,  τ₅_Min, τ₄_Min, τ₁Mac_Min, τ₂Mac_Min, τ₃Mac_Min, τ₁Max_Max, τ₁ηMin_Max, τ₂_Max, τ₃_Max, τ₅_Max, τ₄_Max, τ₁Mac_Max,τ₂Mac_Max, τ₃Mac_Max, Nse_τ, Rmse_τ, Wilmot_τ, Ccc_τ)

		return ksmodelτ 
		end  # function: STRUCT_KSMODEL

end  # module: ksModel
# ............................................................