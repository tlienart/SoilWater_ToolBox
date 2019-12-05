# =============================================================
#		MODULE: HYDRO INITIALIZE
# =============================================================
module infiltInitialize

	import ..hydroStruct, ...psdThetar, ..param, ..timeTransSteady
	export INFILT_INITIALIZE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILT_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILT_INITIALIZE(Tinfilt, N_SoilSelect, N_Infilt, infiltParam, ∑Psd, ∑Infilt_Obs)

			# DERIVING THE STRUCTURE PARAMETERS
				hydroInfilt = hydroStruct.HYDROSTRUCT(N_SoilSelect)

			# CORRECTION FOR θr & θs
				for iSoil=1:N_SoilSelect
					# θr computation
						hydroInfilt.θr[iSoil] = psdThetar.PSD_2_θr_FUNC(iSoil, ∑Psd)

						hydroInfilt.θr[iSoil] =  min(hydroInfilt.θr[iSoil], infiltParam.θ_Ini[iSoil])

					# θs computation
						hydroInfilt.θs[iSoil] = param.hydro.Coeff_Φ_2_θs *  infiltParam.Φ[iSoil]	
				end  # for iSoil=1:N_SoilSelect

			# TIME FLUX CORRECTION
				N_Infilt_Max = maximum(N_Infilt[1:N_SoilSelect])

				Tinfilt_Flux = Array{Float64}(undef, (N_SoilSelect, N_Infilt_Max))
				for iSoil=1:N_SoilSelect
					Tinfilt_Flux[iSoil,1] = ( (Tinfilt[iSoil,1] ^0.5 + (0.0)^0.5) / 2.0 ) ^ 2.0

					for iInfilt=2:N_Infilt[iSoil]
						Tinfilt_Flux[iSoil,iInfilt] = ( (Tinfilt[iSoil,iInfilt] ^ 0.5 + Tinfilt[iSoil,iInfilt-1] ^ 0.5) / 2.0 ) ^ 2.0
					end	
				end #  iSoil=1:N_SoilSelect


				iT_TransStead_Data, T_TransStead_Data = timeTransSteady. ∑INFIlT_2_TIMETRANSSTEADY(N_SoilSelect, N_Infilt, Tinfilt_Flux, ∑Infilt_Obs) 

			# Initializing Infilt		
				Infilt = Array{Float64}(undef, (N_SoilSelect, N_Infilt_Max))

			return Infilt, Tinfilt_Flux, hydroInfilt 
		end  # function: INFILT_INITIALIZE


end # module hydroInitialize
# ............................................................