# =============================================================
#		MODULE: hydroStruct
# =============================================================
module psdStruct
	import ..tool

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Base.@kwdef mutable struct IMP
		ξ1             :: Vector{Float64}
		ξ2             :: Vector{Float64}
		∑Psd_2_ξ2_β1   :: Vector{Float64}
		∑Psd_2_ξ2_β2   :: Vector{Float64}
		∑Psd_2_ξ2_Size :: Vector{Int64}
		Subclay        :: Vector{Float64}
		Psd_2_θr_α1    :: Vector{Float64}
		Psd_2_θr_α2    :: Vector{Float64}
		θr_Psd         :: Vector{Float64}
		Err_θr_Psd     :: Vector{Float64}
		Nse            :: Vector{Float64}

		# FieldName      :: Vector{Symbol} # Need to put
	end # struct IMP


	Base.@kwdef mutable struct CHANG
      ξ1          :: Vector{Float64}
      Psd_2_θr_α1 :: Vector{Float64}
      Psd_2_θr_α2 :: Vector{Float64}
      θr_Psd      :: Vector{Float64}
		Err_θr_Psd  :: Vector{Float64}
		Nse         :: Vector{Float64}
	end # struct CHANG

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	function PSDSTRUCT(NiZ, option, param)
		FieldName = Array{Symbol}(undef, 1) # Need to put

        θr_Psd      = zeros(Float64, NiZ)
        Psd_2_θr_α1 = zeros(Float64, NiZ)
        Psd_2_θr_α2 = zeros(Float64, NiZ)
        Err_θr_Psd  = zeros(Float64, NiZ)
        Nse         = zeros(Float64, NiZ)

		if option.psd.Model⍰ == "IMP"
			ξ1             = zeros(Float64, NiZ)
			ξ2             = zeros(Float64, NiZ)
			∑Psd_2_ξ2_β1   = zeros(Float64, NiZ)
			∑Psd_2_ξ2_β2   = zeros(Float64, NiZ)
			Subclay        = zeros(Float64, NiZ)
			∑Psd_2_ξ2_Size = zeros(Int, NiZ)


		# Initializing
			for iZ=1:NiZ
				ξ1[iZ]             = param.psd.imp.ξ1
				ξ2[iZ]             = 0.0
				∑Psd_2_ξ2_β1[iZ]   = param.psd.imp.∑Psd_2_ξ2_β1
				∑Psd_2_ξ2_β2[iZ]   = param.psd.imp.∑Psd_2_ξ2_β2
				Subclay[iZ]        = param.psd.imp.Subclay
				∑Psd_2_ξ2_Size[iZ] = param.psd.imp.∑Psd_2_ξ2_Size
			end
			
			return paramPsd = IMP(ξ1, ξ2, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, ∑Psd_2_ξ2_Size, Subclay, Psd_2_θr_α1, Psd_2_θr_α2, Err_θr_Psd, Nse, θr_Psd)

			# return paramPsd = tool.readWrite.FIELDNAME_2_STRUCT_VECT(IMP, paramPsd) # Saving the FieldNames
	
		elseif option.psd.Model⍰ == "Chang2019Model"
			ξ1		= zeros(Float64, NiZ)
			θr_Psd  = zeros(Float64, NiZ)

			for iZ=1:NiZ
				ξ1[iZ] = param.psd.chang.ξ1
			end

			return paramPsd = CHANG(ξ1, Psd_2_θr_α1, Psd_2_θr_α2, θr_Psd, Err_θr_Psd, Nse)

			# return paramPsd = tool.readWrite.FIELDNAME_2_STRUCT_VECT(CHANG, paramPsd) # Saving the FieldNames
		end # option.hydro.HydroModel⍰
	end #  function HYDROSTRUCT
	
end # module psdStruct
# ............................................................