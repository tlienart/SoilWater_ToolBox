# =============================================================
#		MODULE: hydroStruct
# =============================================================
module psdStruct

	import ...option, ...param, ...tool

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mutable struct IMP
        ξ1             :: Vector{Float64}
        ∑Psd_2_ξ2_β1   :: Vector{Float64}
        ∑Psd_2_ξ2_β2   :: Vector{Float64}
        Subclay        :: Vector{Float64}
        Psd_2_θr_α1    :: Vector{Float64}
        Psd_2_θr_α2    :: Vector{Float64}
        ∑Psd_2_ξ2_Size :: Vector{Int64}
        Nse            :: Vector{Float64}

		FieldName		::Vector{Symbol} # Need to put
	end # struct IMP


	mutable struct CHANG
        ξ1        :: Vector{Float64}
        Nse       :: Vector{Float64}
        FieldName :: Vector{Symbol} # Need to put
	end # struct CHANG

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	function PSDSTRUCT(N_SoilSelect)
		FieldName = Array{Symbol}(undef, 1) # Need to put

		if option.psd.Model == "IMP"
            ξ1             = zeros(Float64, N_SoilSelect)
            ∑Psd_2_ξ2_β1   = zeros(Float64, N_SoilSelect)
            ∑Psd_2_ξ2_β2   = zeros(Float64, N_SoilSelect)
            Subclay        = zeros(Float64, N_SoilSelect)
            Psd_2_θr_α1    = zeros(Float64, N_SoilSelect)
            Psd_2_θr_α2    = zeros(Float64, N_SoilSelect)
            ∑Psd_2_ξ2_Size = zeros(Int, N_SoilSelect)
            Nse            = zeros(Float64, N_SoilSelect)

			# Initializing
			for iSoil=1:N_SoilSelect
				ξ1[iSoil] 				= param.psd.imp.ξ1
				∑Psd_2_ξ2_β1[iSoil] 	= param.psd.imp.∑Psd_2_ξ2_β1
				∑Psd_2_ξ2_β2[iSoil] 	= param.psd.imp.∑Psd_2_ξ2_β2
				Subclay[iSoil] 			= param.psd.imp.Subclay
				∑Psd_2_ξ2_Size[iSoil] 	= param.psd.imp.∑Psd_2_ξ2_Size
			end
			
			psdparam = IMP(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, Psd_2_θr_α1, Psd_2_θr_α2, ∑Psd_2_ξ2_Size, Nse, FieldName)

			return psdparam = tool.readWrite.FIELDNAME_2_STRUCT(IMP, psdparam) # Saving the FieldNames
	
		elseif option.psd.Model == "Chang2019Model"
			ξ1		= Array{Float64}(undef, (N_SoilSelect))

			for iSoil=1:N_SoilSelect
				ξ1[iSoil] 				= param.psd.chan.ξ1
			end

			psdparam = CHANG(ξ1, FieldName)

			return psdparam = tool.readWrite.FIELDNAME_2_STRUCT(CHANG, Nse, psdparam) # Saving the FieldNames
		end # option.hydro.HydroModel
	end #  function HYDROSTRUCT
	
end # module psdStruct
# ............................................................