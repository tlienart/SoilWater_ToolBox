module climate
	import Dates: value, DateTime
	export CLIMATE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :   CLIMATE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CLIMATE(clim, optionHypix)

			fill(0.0::Float64, clim.N_Climate)

			∑Pr_Climate      = fill(0.0::Float64, clim.N_Climate)
			∑Pet_Climate	  = fill(0.0::Float64, clim.N_Climate)
			∑T_Climate       = fill(0.0::Float64, clim.N_Climate)
			Temp             = fill(0.0::Float64, clim.N_Climate)

			 # Taking into acount that ΔT is the difference of time of T[iT]-T[iT-1]
				∑Pr_Climate[1]  = 0.0::Float64
				∑Pet_Climate[1] = 0.0::Float64
				∑T_Climate[1]   = 0.0::Float64
				Temp[1]         = 0.0::Float64
		 
			 for iT = 2:clim.N_Climate
				#Computing cumulative time 
				∑T_Climate[iT] = value(clim.Date[iT] - clim.Date[1]) / 1000

				# Cumulative
				if !(optionHypix.RainfallInterception) && optionHypix.TopBoundary⍰ ≠ "Ψ"
               ∑Pr_Climate[iT]  = ∑Pr_Climate[iT-1] + clim.Pr[iT]  # Cumulate Pr
               ∑Pet_Climate[iT] = ∑Pet_Climate[iT-1] + clim.Pet[iT] # Cumulative Potential evaporation
				end
				
				# Tempoerature no change 
					Temp[iT] = clim.Temp[iT]
			end # for

			N_∑T_Climate = Int(∑T_Climate[end])

		return ∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp
		end # function CLIMATE

end # module climate
# ...........................................................................................