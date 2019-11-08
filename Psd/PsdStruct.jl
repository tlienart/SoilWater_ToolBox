# =============================================================
#		MODULE: hydroStruct
# =============================================================
module psdStruct

	import ...option, ..param

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mutable struct IMP
		ξ1 				:: Vector{Float64}
		∑Psd_2_ξ2_β1 	:: Vector{Float64}
		∑Psd_2_ξ2_β2 	:: Vector{Float64}
		Subclay 		:: Vector{Float64}
		Psd_2_θr_α1 	:: Vector{Float64}
		Psd_2_θr_α2 	:: Vector{Float64}
		∑Psd_2_ξ2_Size 	:: Vector{Int64}

	end # struct IMP

	mutable struct CHANG
		ξ1 				:: Vector{Float64}
	end # struct CHANG

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	function PSDSTRUCT(N_SoilSelect)
		if option.psd.Model == "IMP"
			ξ1 				= Array{Float64}(undef, (N_SoilSelect))
			∑Psd_2_ξ2_β1 	= Array{Float64}(undef, (N_SoilSelect))
			∑Psd_2_ξ2_β2 	= Array{Float64}(undef, (N_SoilSelect))
			Subclay 		= Array{Float64}(undef, (N_SoilSelect))
			Psd_2_θr_α1 	= Array{Float64}(undef, (N_SoilSelect))
			Psd_2_θr_α2 	= Array{Float64}(undef, (N_SoilSelect))
			∑Psd_2_ξ2_Size 	= Array{Int64}(undef, (N_SoilSelect))

			# Initializing
			for iSoil=1:N_SoilSelect
				ξ1[iSoil] 				= param.psd.imp.ξ1
				∑Psd_2_ξ2_β1[iSoil] 	= param.psd.imp.∑Psd_2_ξ2_β1
				∑Psd_2_ξ2_β2[iSoil] 	= param.psd.imp.∑Psd_2_ξ2_β2
				Subclay[iSoil] 			= param.psd.imp.Subclay
				∑Psd_2_ξ2_Size[iSoil] 	= param.psd.imp.∑Psd_2_ξ2_Size
			end
			
			return psdparam = IMP(ξ1, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, Subclay, Psd_2_θr_α1, Psd_2_θr_α2, ∑Psd_2_ξ2_Size )
	
		elseif option.psd.Model == "Chang2019Model"
			ξ1		= Array{Float64}(undef, (N_SoilSelect))

			for iSoil=1:N_SoilSelect
				ξ1[iSoil] 				= param.psd.chan.ξ1
			end

			return psdparam = CHANG(ξ1)
		end # option.HydroModel
	end #  function HYDROSTRUCT
	
	
end # module psdStruct
# ............................................................