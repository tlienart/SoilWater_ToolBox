# =============================================================
#		MODULE: veg
# =============================================================
module rootWaterUptake
	export ROOT_WATER_UPTAKE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ROOT_WATER_UPTAKE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 function ROOT_WATER_UPTAKE(CropCoeficient::Float64, iT::Int64, N_iRoot::Int64, optionHypix, veg, ΔPet_Transp::Float64, ΔRootDensity::Vector{Float64}, ΔSink::Matrix{Float64}, Ψ::Matrix{Float64})

		if optionHypix.RootWaterUptakeComp
			for iZ = 1:N_iRoot
				RootCompensation = rootWaterUptake.rootCompensation.ROOT_COMPENSATION(iT, iZ, N_iRoot, veg, ΔRootDensity, Ψ)

				ΔSink[iT,iZ] = CropCoeficient * ΔPet_Transp * rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(iT, iZ, veg, Ψ) * RootCompensation
			end # for
		else # optionHypix.RootWaterUptakeComp

			for iZ = 1:N_iRoot
				ΔSink[iT,iZ] = CropCoeficient * ΔPet_Transp * ΔRootDensity[iZ] * rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(iT, iZ, veg, Ψ)
			end  # for
		end # if:

	return ΔSink
	end  # function ROOT_WATER_UPTAKE
	

	# =============================================================
	#		MODULE: rootDistribution
	# =============================================================
	module rootDistribution
		import Optim
		export ROOT_DENSITY, N_IROOT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : COMPUTTING N OF ROOTING DEPTH
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function N_IROOT(NiZ::Int64, veg, Z::Vector{Float64})::Int64
				iZ = 1
				while iZ ≤ NiZ && veg.Zroot ≥ Z[iZ]
					iZ +=1
				end
			return iZ - 1
			end  # function: N_IROOT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : rootdistribution
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		"""Optimizing the value of RootDensityParam such that for the top soil of depth Zroot_Top we have ΔRdf_Top of roots and that it complies with the maximum depth"""
		function ROOT_DENSITY(discret, N_iRoot::Int64, veg, Z::Vector{Float64})
		
			# Deriving the root density parameters
				function OF_ROOTDENSITY(veg, N_iRoot, RootDensityParam, Z)
					RootDensityParam2 = (1.0 - 10.0 ^ -RootDensityParam)
					# RootDensityParam2 ^ 0.0 = 1
					return Of = abs((1.0 - RootDensityParam2 ^ veg.Zroot_Top) / (1.0 - RootDensityParam2 ^ veg.Zroot) - veg.ΔRdf_Top)
				end
				Optimization =  Optim.optimize(RootDensityParam ->  OF_ROOTDENSITY(veg, N_iRoot, RootDensityParam, Z), 0.0, 20.0, Optim.GoldenSection())
				RootDensityParam = (Optim.minimizer(Optimization)[1])
				
			# Computing the RootDensityFunction
				ΔRootDensity = fill(0.0::Float64, N_iRoot)
				RootDensityParam₂ = (1.0 - 10.0 ^ -RootDensityParam)
				for iZ = 1:N_iRoot-1
				 	ΔRootDensity[iZ] = (RootDensityParam₂ ^ discret.Z_CellUp[iZ] - RootDensityParam₂ ^ Z[iZ]) / (1.0 - RootDensityParam₂ ^ veg.Zroot)
				end
				ΔRootDensity[N_iRoot] = (RootDensityParam₂ ^ discret.Z_CellUp[N_iRoot] - RootDensityParam₂ ^ veg.Zroot) / (1.0 - RootDensityParam₂ ^ veg.Zroot)

		return ΔRootDensity		
		end  # function rootdistribution
		
	end # module rootdistribution
	# ............................................................



	# =============================================================
	#		MODULE: rootWaterUptake
	# =============================================================
	module stressReduction 
		export WATER_STRESS_FUNCTION

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : WATER_STRESS_FUNCTION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function WATER_STRESS_FUNCTION(iT::Int64, iZ::Int64, veg, Ψ::Matrix{Float64})
				if veg.Ψfeddes3 < Ψ[iT-1,iZ] < veg.Ψfeddes4
					return RootWaterUptake = (Ψ[iT-1,iZ] - veg.Ψfeddes4) / (veg.Ψfeddes3 - veg.Ψfeddes4)
				
				elseif Ψ[iT-1,iZ] ≤ veg.Ψfeddes1  || Ψ[iT-1,iZ] ≥ veg.Ψfeddes4
					return RootWaterUptake = 0.0
					
				elseif veg.Ψfeddes2 ≤ Ψ[iT-1,iZ] ≤ veg.Ψfeddes3
					return RootWaterUptake = 1.0
					
				elseif veg.Ψfeddes1 < Ψ[iT-1,iZ] < veg.Ψfeddes2
					return RootWaterUptake = (Ψ[iT-1,iZ] - veg.Ψfeddes1) / (veg.Ψfeddes2 - veg.Ψfeddes1)
				end # Ψ >= veg.Ψfeddes3 && Ψ <= veg.Ψfeddes4 	
			end  # function WATER_STRESS_FUNCTION
			
	end  # module rootWaterUptake
	# ............................................................



	# =============================================================
	#		MODULE: rootcompensation
	# =============================================================
	module rootCompensation
		import ..rootWaterUptake
		export ROOT_COMPENSATION

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ROOTCOMPENSATION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ROOT_COMPENSATION(iT::Int64, iZ::Int64, N_iRoot::Int64, veg, ΔRootDensity::Vector{Float64}, Ψ::Matrix{Float64})
				
				#Compute the denominator which is used to normalize
				∑RootCompensation = 0.0
				for iZ = 1:N_iRoot
					∑RootCompensation += rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(iT, iZ, veg, Ψ) * ΔRootDensity[iZ] ^ veg.RootWaterUptakeComp
				end
				
				# Compute the Compensation
				if ∑RootCompensation > 0.0-0 
					# return RootCompensation = (rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(iT, iZ, veg, Ψ) * ΔRootDensity[iZ] ^ (veg.RootWaterUptakeComp - 1.0)) / ∑RootCompensation

					return RootCompensation = (rootWaterUptake.stressReduction.WATER_STRESS_FUNCTION(iT, iZ, veg, Ψ) * ΔRootDensity[iZ] ^ veg.RootWaterUptakeComp) / ∑RootCompensation
				else
					return RootCompensation = 1.0
				end
			end  # function ROOT-COMPENSATION

	end  # module rootcompensation
	# ............................................................
	
end  # module veg
# ............................................................