# =============================================================
#		module: zobs
# =============================================================
module zobs

	import Dates: value, DateTime
	export θOBS
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ZOBS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ZOBS(calibr, clim, discret, Z)

			# CHECKING DATA CONSISTENCY
				if calibr.Date[1] < clim.Date[2]
					error("\n Hypix error: Starting date of calibr  $(calibr.Date[1]) < starting date of climate data $(clim.Date[1])")
				end # Error checking

				# Checking the celLs
				if calibr.Z[calibr.Ndepth] > Z[discret.N_iZ]
					error("\n Hypix error: depth of measured θ deeper than the max depth of discretisation: calibr.Z[calibr.Ndepth] > discret.Z[discret.N_iZ]") 
				end

			# COMPUTING CUMULATIVE TIME OF
				Start_Date = clim.Date[1] 
				for iT=1:calibr.N_iT
					calibr.∑T[iT] = value(calibr.Date[iT] - Start_Date) / 1000
				end  # for it=1:calibr.N_iT

				# TRANSFORM THE DEPTH OF MEASURED Θ -> CELL DEPTH
				for iDepth = 1:calibr.Ndepth
					for iZ = 1:discret.N_iZ
						if iZ == 1
							if 0.0 ≤ calibr.Z[iDepth] ≤ Z[1]
								calibr.iZobs[iDepth] = 1
								break  
							end  # if: discret.Z_CellUp
						elseif iZ ≠ 1
							if Z[iZ-1] ≤ calibr.Z[iDepth] ≤ Z[iZ]
								calibr.iZobs[iDepth] = iZ
								break  
							end  # if: discret.Z_CellUp
						end # if iZ == 1
					end # iZ = 2:discret.N_iZ						
				end  # for iDepth = 1:calibr.Ndepth
		
			return calibr
		end  # function: θOBS

end # module: Zobs