# =============================================================
#		MODULE: HYDRO INITIALIZE
# =============================================================
module infiltInitialize
	import ..hydroStruct, ..psdThetar, ..timeTransSteady, ..infiltStruct
	export INFILT_INITIALIZE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILT_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILT_INITIALIZE(∑Infilt_Obs, ∑Psd, hydroInfilt, infiltParam, N_Infilt, N_iZ, Tinfilt)

			# CORRECTION FOR θr & θs
				for iZ=1:N_iZ
					# θr computation
						# if option.run.IntergranularMixingPsd
							hydroInfilt.θr[iZ] = psdThetar.PSD_2_θr_FUNC(∑Psd, hydroInfilt, iZ, param)
						# else
						# 	hydroInfilt.θr[iZ] = param.hydro.θr
						# end # option.run.IntergranularMixingPsd

					# θr < θ_Ini
						infiltParam.θ_Ini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θ_Ini[iZ])

					# θs computation
						hydroInfilt.θs[iZ] = param.hydro.Coeff_Φ_2_θs * hydroInfilt.Φ[iZ]

					# Checking for θ_Ini
						infiltParam.θ_Ini[iZ] = min(infiltParam.θ_Ini[iZ], hydroInfilt.θs[iZ] - 0.0015) 
				end  # for iZ=1:N_iZ

			# TIME FLUX CORRECTION
				N_Infilt_Max = maximum(N_Infilt[1:N_iZ])

				T = Array{Float64}(undef, (N_iZ, N_Infilt_Max))
				for iZ=1:N_iZ
					# T[iZ,1] = ( (Tinfilt[iZ,1] ^0.5 + (0.0)^0.5) / 2.0 ) ^ 2.0
					T[iZ,1] = Tinfilt[iZ,1]

					for iInfilt=2:N_Infilt[iZ]
						T[iZ,iInfilt] = ( (Tinfilt[iZ,iInfilt] ^ 0.5 + Tinfilt[iZ,iInfilt-1] ^ 0.5) / 2.0 ) ^ 2.0
					end	
				end #  iZ=1:N_iZ

			# STRUCTURE OF INFILTRATION PARAMETERS
				infiltOutput = infiltStruct.INFILTSTRUCT(N_iZ)

			# DETERMENING WHEN STEADY STATE OCCURES
				infiltOutput = timeTransSteady.∑INFIlT_2_TIMETRANSSTEADY(T, N_iZ, N_Infilt, infiltOutput, ∑Infilt_Obs) 

			# Initializing Infilt		
				∑Infilt_3D = Array{Float64}(undef, (N_iZ, N_Infilt_Max))
				∑Infilt_1D = Array{Float64}(undef, (N_iZ, N_Infilt_Max))

			return T, infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D

		end  # function: INFILT_INITIALIZE


end # module hydroInitialize
# ............................................................