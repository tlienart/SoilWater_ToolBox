module timeTransSteady
	import ..stats, ...param
	export ∑INFIlT_2_TIMETRANSSTEADY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∑INFIlT_2_TIMETRANSSTEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		function  ∑INFIlT_2_TIMETRANSSTEADY(N_SoilSelect, N_Infilt, T, ∑Infilt_Obs; N_LastInfiltPoint_Select=3) 
		
            T_TransSteady_Data  = Array{Float64}(undef, (N_SoilSelect))
            iT_TransSteady_Data = Array{Int64}(undef, (N_SoilSelect))
		
			# FOR EVERy SOIL
			for iSoil=1:N_SoilSelect

				println("Max $iSoil, $(N_Infilt[iSoil])")

                Flag_TransSteady     = false
                ∑Infilt_Model        = Array{Float64}(undef, (N_Infilt[iSoil]))
				ΔErr                 = Array{Float64}(undef, (N_Infilt[iSoil]))
				
				# Selecting the number of last points, there need to be at least 3 points
				N_LastInfiltPoint = min(max(N_Infilt[iSoil] - 3,3), N_LastInfiltPoint_Select)

				Intercept, Slope = stats.LINEAR_REGRESSION(T[N_LastInfiltPoint:N_Infilt[iSoil]], ∑Infilt_Obs[N_LastInfiltPoint:N_Infilt[iSoil]])
				
				# Determining starting from the last points
				for i in N_Infilt[iSoil] - N_LastInfiltPoint:-1:1
				
					if i-1 >= 1
						# Determening if a linear equation is valid
						# Intercept, Slope = stats.LINEAR_REGRESSION(T[i:N_Infilt[iSoil]], ∑Infilt_Obs[i:N_Infilt[iSoil]])
		
						∑Infilt_Model[i-1] = T[i-1] * Slope + Intercept

						ΔErr[i-1] = abs( ∑Infilt_Obs[iSoil,i-1] - ∑Infilt_Model[i-1]) / (T[iSoil,i] - T[iSoil,i-1])
				
						if ΔErr[i-1] >= param.infilt.ΔErrMax_TransSteady && !Flag_TransSteady # To catch only the very beginning
							Flag_TransSteady = true
							iT_TransSteady_Data[iSoil] = i-1
							break
						end # ΔErr[i-1] >= param.infilt.ΔErrMax_TransSteady && !Flag_TransSteady 
					else
						iT_TransSteady_Data[iSoil] = 1
						break	
					end # 	if i-1 >= 1

				end # for i in N_Infilt[iSoil] - N_LastInfiltPoint:-1:1
			
				T_TransSteady_Data[iSoil] = T[iSoil,iT_TransSteady_Data[iSoil]]
		
				println("$iSoil iTransSteady= , $(iT_TransSteady_Data[iSoil]), $(T[iSoil,iT_TransSteady_Data[iSoil]])")
		
			end # for iSoil=1:N_SoilSelect

			return iT_TransSteady_Data, T_TransSteady_Data

		end # function: INFIlTOBS_2_iTIME_TRANS_STEADy

	
end  # macro timeTransSteady

