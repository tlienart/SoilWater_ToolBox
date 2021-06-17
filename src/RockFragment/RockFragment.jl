# =============================================================
#		module: rockFragment
# =============================================================
module rockFragment
	import Polynomials

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  ρᵦ_2_Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			if option.rockFragment.RockInjectedIncluded == :InjectRock
				return Φ = rockFragment.injectRock.ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)

			elseif option.rockFragment.RockInjectedIncluded == :Included
				return Φ = rockFragment.included.ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			end
		end  # function: function ρᵦ_2_Φ
	

	# =============================================================
	#		module: injectRock
	# =============================================================
	module injectRock
		export CORECTION_Φ!, CORECTION_θΨ!, ρᵦ_2_Φ, CORECTION_Φ!

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : CORECTION_θΨ!
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function CORECTION_θΨ!(N_iZ, N_θΨobs, RockFragment, θ_θΨobs)
				for iZ = 1:N_iZ 
					for iθ=1:N_θΨobs[iZ]
						θ_θΨobs[iZ,iθ] =  θ_θΨobs[iZ,iθ] * (1.0 - RockFragment[iZ])
					end # for iθ=1:N_θΨobs[iZ]
				end #  for iZ = 1:N_iZ	
			return  θ_θΨobs
			end  # function: STONECORRECTION_NONEWETABLE

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION :  CORECTION_Φ!
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function CORECTION_Φ!(N_iZ, option, RockFragment, Φ)		
				if option.run.RockCorection
					for iZ = 1:N_iZ 
						Φ[iZ] = Φ[iZ] * (1.0 - RockFragment[iZ])
					end
				end
			return  Φ
			end  # function: STONECORRECTION_NONEWETABLE


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  ρᵦ_2_Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			Φ = fill(0.0::Float64, N_iZ)

			for iZ=1:N_iZ
				if option.run.RockCorection
					Φ[iZ] = (1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]) * (1.0 - RockFragment[iZ])
				else
					Φ[iZ] = 1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]
				end
			end # for
		return Φ
		end  # function: Φ			
	end  # module: injectRock
	# ............................................................


	# =============================================================
	#		module: included
	# =============================================================
	module included
		export ρᵦ_2_Φ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION :  ρᵦ_2_Φ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ρᵦ_2_Φ(N_iZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
				Φ = fill(0.0::Float64, N_iZ)
				for iZ=1:N_iZ
					if option.run.RockCorection					
						Φ[iZ] = 1.0 - (RockFragment[iZ] * ρᵦ_Soil[iZ] / ρₚ_Rock[iZ]) - ((1.0 - RockFragment[iZ]) * ρᵦ_Soil[iZ] / ρₚ_Fine[iZ])
					else
						Φ[iZ] = 1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]
					end
				end # for
			return Φ
			end  # function: Φ
		
	end  # module: included
	# ............................................................
	
end  # module: rockFragment
# ............................................................