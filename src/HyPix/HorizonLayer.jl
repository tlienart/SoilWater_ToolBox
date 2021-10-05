# =============================================================
#		momodule horizonLayer
# =============================================================
module horizonLayer

	import ..hydroStruct
	export HYDROHORIZON_2_HYDRO

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :   HORIZON_2_LAYER
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ::Int64, option)

		hydro = hydroStruct.HYDROSTRUCT(option.hyPix, NiZ)

		# Field names of the structure
			FieldName_Array = propertynames(hydroHorizon)
		
		# looping through every fieldnames of the structure
			for FieldName in FieldName_Array
				Value_Array = getfield(hydroHorizon, FieldName)
				
				Vector = fill(0.0::Float64, NiZ)
				for iZ = 1:NiZ
					Vector[iZ] = getfield(hydroHorizon, FieldName)[Layer[iZ]]
				end
				setfield!(hydro, Symbol(FieldName), Vector)
			end

	return hydro
	end # function HORIZON_2_LAYER

end  # mmodule horizonLayer
# ............................................................