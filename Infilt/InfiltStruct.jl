# =============================================================
#		MODULE: infiltStruct
# =============================================================
module infiltStruct
	import ..option, ..tool

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct INFILT
            Sorptivity         :: 	Vector{Float64}
            iT_TransStead_Data :: 	Vector{Float64}
            T_TransStead_Data  :: 	Vector{Float64}
            Nse_Trans          ::	Vector{Float64}
			Nse_Steady         ::	Vector{Float64}
			Nse			       ::	Vector{Float64}
		
            FieldName          ::	Vector{Symbol} # Need to put
		end # struct KOSUGI

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function INFILTSTRUCT()
        FieldName          = Array{Symbol}(undef, 1) # Need to put
        Sorptivity         = zeros(Float64, N_SoilSelect)
        iT_TransStead_Data = zeros(Float64, N_SoilSelect)
        T_TransStead_Data  = zeros(Float64, N_SoilSelect)
        Nse_Trans          = zeros(Float64, N_SoilSelect)
        Nse_Steady         = zeros(Float64, N_SoilSelect)
		Nse                = zeros(Float64, N_SoilSelect)
		
		infilt = INFILT(Sorptivity, iT_TransStead_Data, T_TransStead_Data, Nse_Trans, Nse_Steady, Nse, FieldName)

		return infilt = tool.readWrite.FIELDNAME_2_STRUCT(INFILT, infilt) # Saving the FieldNames

	end  # function: INFILTSTRUCT
	
end # module infiltStruct
# ............................................................