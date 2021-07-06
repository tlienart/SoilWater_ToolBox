# =============================================================
#		module: thetaObs
# =============================================================
module thetaObs
	import Dates: value, DateTime
	export thetaObs
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΘOBS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΘOBS(obsTheta, clim, discret, Z)

			# CHECKING DATA CONSISTENCY
				if obsTheta.Date[1] < clim.Date[2]
					error("\n Hypix error: Starting date of obsTheta  $(obsTheta.Date[1]) < starting date of climate data $(clim.Date[1])")
				end # Error checking

				# Checking the celLs
				if obsTheta.Z[obsTheta.Ndepth] > Z[discret.N_iZ]
					error("\n Hypix error: depth of measured θ deeper than the max depth of discretisation: obsTheta.Z[obsTheta.Ndepth] > discret.Z[discret.N_iZ]") 
				end

			# COMPUTING CUMULATIVE TIME
				Start_Date = clim.Date[1] 
				for iT=1:obsTheta.N_iT
					obsTheta.∑T[iT] = value(obsTheta.Date[iT] - Start_Date) / 1000
				end  # for it=1:obsTheta.N_iT

				# TRANSFORM THE DEPTH OF MEASURED Θ -> CELL DEPTH
				for iDepth = 1:obsTheta.Ndepth
					for iZ = 1:discret.N_iZ
						if iZ == 1
							if 0.0 ≤ obsTheta.Z[iDepth] ≤ Z[1]
								obsTheta.ithetaObs[iDepth] = 1
								break  
							end  # if: discret.Z_CellUp
						elseif iZ ≠ 1
							if Z[iZ-1] ≤ obsTheta.Z[iDepth] ≤ Z[iZ]
								obsTheta.ithetaObs[iDepth] = iZ
								break  
							end  # if: discret.Z_CellUp
						end # if iZ == 1
					end # iZ = 2:discret.N_iZ						
				end  # for iDepth = 1:obsTheta.Ndepth
		
		return obsTheta
		end  # function: θOBS

end # module: thetaObs