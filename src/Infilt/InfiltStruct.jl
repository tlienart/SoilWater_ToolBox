# =============================================================
#		MODULE: infiltStruct
# =============================================================
module infiltStruct
	import ..tool
	export INFILTSTRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Base.@kwdef mutable struct INFILT
			Sorptivity         :: 	Vector{Float64}
			iT_TransSteady_Data :: 	Vector{Int64}
			T_TransSteady_Data  :: 	Vector{Float64}
			Nse_Trans          ::	Vector{Float64}
			Nse_Steady         ::	Vector{Float64}
			Nse			       ::	Vector{Float64}
		end # struct KOSUGI

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTSTRUCT(NiZ)
         FieldName           = Array{Symbol}(undef, 1) # Need to put
         Sorptivity          = zeros(Float64, NiZ)
         iT_TransSteady_Data = zeros(Int64, NiZ)
         T_TransSteady_Data  = zeros(Float64, NiZ)
         Nse_Trans           = zeros(Float64, NiZ)
         Nse_Steady          = zeros(Float64, NiZ)
         Nse                 = zeros(Float64, NiZ)
			
			return infiltOutput = INFILT(Sorptivity, iT_TransSteady_Data, T_TransSteady_Data, Nse_Trans, Nse_Steady, Nse)

		end  # function: INFILTSTRUCT

end # module infiltStruct
# ............................................................