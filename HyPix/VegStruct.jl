# =============================================================
#		module: veg
# =============================================================
module vegStruct
	import ..tool
	export VEGSTRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : veg
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mutable struct VEG	
      Zroot                :: Float64
      ΔRdf_Top             :: Float64
      Zroot_Top            :: Float64
      Lai                  :: Float64
      Lai_Min              :: Float64
      Lai_Max              :: Float64
      CropCoeficient       :: Float64
      CropCoeficient_Min   :: Float64
      CropCoeficient_Max   :: Float64
      ExtinctCoefRadiation :: Float64
      Ψfeddes1             :: Float64
      Ψfeddes2             :: Float64
      Ψfeddes3             :: Float64
      Ψfeddes4             :: Float64
      RootWaterUptakeComp  :: Float64
      Zevapo               :: Float64
      Sint_Lai             :: Float64
      Sint_Sat             :: Float64

		FieldName 				::	Vector{Symbol} # Need to pu
	end # struct VEG

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : VEG
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function VEGSTRUCT()
		FieldName            = Array{Symbol}(undef, 1) # Need to put
		
      Zroot                = 0.0
      ΔRdf_Top             = 0.0
      Zroot_Top            = 0.0
      Lai                  = 0.0
      Lai_Min              = 0.0
      Lai_Max              = 0.0
      CropCoeficient       = 0.0
      CropCoeficient_Min   = 0.0
      CropCoeficient_Max   = 0.0
      ExtinctCoefRadiation = 0.0
      Ψfeddes1             = 0.0
      Ψfeddes2             = 0.0
      Ψfeddes3             = 0.0
      Ψfeddes4             = 0.0
      RootWaterUptakeComp  = 0.0
      Zevapo               = 0.0
      Sint_Lai             = 0.0
      Sint_Sat             = 0.0

		veg = VEG(Zroot, ΔRdf_Top, Zroot_Top, Lai, Lai_Min, Lai_Max, CropCoeficient, CropCoeficient_Min, CropCoeficient_Max, ExtinctCoefRadiation, Ψfeddes1, Ψfeddes2, Ψfeddes3, Ψfeddes4, RootWaterUptakeComp, Zevapo, Sint_Lai, Sint_Sat, FieldName)
		
		return veg = tool.readWrite.FIELDNAME_2_STRUCT_VECT(VEG, veg) # Saving the FieldNames
	end # function VEG
	
end  # module vegStruct
# ............................................................		