# =============================================================
#		module: rockFragment
# =============================================================
module rockFragment

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  ρᵦ_2_Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			if option.rockFragment.RockInjectedIncluded == :Injected
				return Φ = rockFragment.injectRock.ρᵦ_2_Φ(N_SoilSelect, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)

			elseif  option.rockFragment.RockInjectedIncluded == :Included
				return Φ = rockFragment.injectRock.ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			end
		end  # function: function ρᵦ_2_Φ
	

	# =============================================================
	#		module: injectRock
	# =============================================================
	module injectRock
		export CORECTION_Φ!, CORECTION_θΨ!, ρᵦ_2_Φ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : STONECORRECTION_NONEWETABLE
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function CORECTION_θΨ!(N_SoilSelect, N_θΨ, RockFragment, θ_θΨ)
				for iZ = 1:N_SoilSelect 
					for iθ=1:N_θΨ[iZ]
						θ_θΨ[iZ,iθ] =  θ_θΨ[iZ,iθ] * (1.0 - RockFragment[iZ])
					end # for iθ=1:N_θΨ[iZ]
				end #  for iZ = 1:N_SoilSelect
				
			return  θ_θΨ
			end  # function: STONECORRECTION_NONEWETABLE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  ρᵦ_2_Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρᵦ_2_Φ(N_SoilSelect, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			Φ = fill(0.0::Float64, N_SoilSelect)

			for iZ=1:N_SoilSelect
				Φ[iZ] = (1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]) * (1.0 - RockFragment[iZ])
			end # for
		return Φ
		end  # function: Φ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  CORECTION_Φ!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CORECTION_Φ!(N_SoilSelect, RockFragment, Φ)	
			for iZ=1:N_SoilSelect
				Φ[iZ] = Φ[iZ] * (1.0 - RockFragment[iZ])
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
			function ρᵦ_2_Φ(N_SoilSelect, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
				Φ = fill(0.0::Float64, N_SoilSelect)
				for iZ=1:N_SoilSelect
					if option.run.RockCorection					
						Φ[iZ] = 1.0 - (RockFragment[iZ] * ρᵦ_Soil[iZ] / ρₚ_Rock[iZ]) - ((1.0 - RockFragment[iZ]) * ρᵦ_Soil[iZ] / ρₚ_Fine[iZ])
					else
						Φ[iZ] = (1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ])
					end
				end # for
			return Φ
			end  # function: Φ
		
	end  # module: included
	# ............................................................
	
end  # module: rockFragment
# ............................................................