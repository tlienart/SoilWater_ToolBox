# =============================================================
#		module: rockFragment
# =============================================================
module rockFragment
	export CORECTION_θΨ!

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
	
end  # module: rockFragment
# ............................................................