# =============================================================
#		MODULE: psdInitiate
# =============================================================
module psdInitialize
	import ..psdStruct, ..psdThetar, ..param
	export PSD_INITIALIZE
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PSD_INITIALIZE(N_Psd, N_SoilSelect, ∑Psd, Φ_Psd, hydro, hydroPsd)			
		# Compute new N_Psd to take into account when ∑Psd 
		# Correction for N_Psd such that to determine the real maximum Rpart size
			for iSoil=1:N_SoilSelect
				N_Psd_New = 1
				for iPsd = 1:N_Psd[iSoil]
					if ∑Psd[iSoil, iPsd] >= 0.99999
						N_Psd_New = iPsd
						break
					end
				end #  iPsd = 1:N_Psd[iSoil]

				N_Psd[iSoil] = N_Psd_New
			end # looping over soils

			# MAximum number of Psd data
			N_Psd_Max = maximum(N_Psd[1:N_SoilSelect])


		# Compute PSD from ∑PSD
			Psd = zeros(Float64, (N_SoilSelect, N_Psd_Max))
			for iSoil=1:N_SoilSelect
				Psd[iSoil, 1:N_Psd[iSoil]] = ∑PSD_2_PSD(∑Psd[iSoil,1:N_Psd[iSoil]], N_Psd[iSoil])
			end # iSoil=1:N_SoilSelect

		# DERIVING THE STRUCTURE PARAMETERS
			paramPsd = psdStruct.PSDSTRUCT(N_SoilSelect)

		# COMPUTING θr FROM PSD DATA
			paramPsd = psdThetar.PSD_2_θr(N_SoilSelect, ∑Psd, hydro, paramPsd)
			
		# COMPUTING θs FROM TOTAL POROSITY
			θs_Psd = Array{Float64}(undef, N_SoilSelect)
			for iSoil=1:N_SoilSelect
				θs_Psd[iSoil] = param.hydro.Coeff_Φ_2_θs * Φ_Psd[iSoil]
				hydroPsd.θs[iSoil] = θs_Psd[iSoil]
				hydroPsd.θr[iSoil] = paramPsd.θr_Psd[iSoil]
			end 
		
		return N_Psd, N_Psd_Max, Psd, θs_Psd, hydroPsd, paramPsd
	end  # function: PSD_INITIALIZE


	# =========================================
	#          ∑PSD -> PSD
	# =========================================
		function ∑PSD_2_PSD(∑Psd, N_Psd)
			Psd = zeros(Float64, N_Psd)

			Psd[1] = ∑Psd[1]
			for iRpart =2:N_Psd
				Psd[iRpart] = ∑Psd[iRpart] - ∑Psd[iRpart-1]
			end
			return Psd
		end # function ∑PSD_2_PSD
	
end  # module psdInitiate
# ............................................................
