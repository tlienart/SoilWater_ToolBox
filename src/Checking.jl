# =============================================================
#		module: checking
# =============================================================
module checking

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : CHECKING
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function CHECKING(option, optionₘ, optim)

        # ------------ CHECKING HydroLabθΨ---------------------

		  if "Ks" ∈ optim.ParamOpt && !(option.data.Kθ) 
			error("*** If Ks ∈ optim.ParamOpt ⇒option.data.θΨ ***")

			elseif option.run.HydroLabθΨ⍰ ≠ :No && !option.data.θΨ
				error("*** If option.run.HydroLabθΨ⍰ ⇒option.data.θΨ ***")

			elseif option.run.RockCorection && option.rockFragment.RockInjectedIncluded⍰ ==:InjectRock && !( option.data.BulkDensity && option.data.θΨ)
				error("*** If option.run.RockCorrection && option.rockFragment.RockInjectedIncluded⍰ ==:InjectRock ⇒ option.data.BulkDensity OR option.data.θΨ ***")

			elseif optionₘ.HydroModel⍰==:Kosugi && "θsMacMat" ∈ optim.ParamOpt
				error("*** optionₘ.HydroModel⍰==:Kosugi && optionₘ.HydroModel⍰==:Bimodal THAN optionₘ.HydroModel⍰ == :Φ ***")
					
			elseif optionₘ.θrOpt⍰≠:Opt && "θr" ∈ optim.ParamOpt
				error("*** optionₘ.θrOpt≠:Opt && θr ∈ ParamOpt THAN do not optimize θr ***")
			              
			elseif optionₘ.σ_2_Ψm⍰ ==:UniqueRelationship && "Ψm" ∈ optim.ParamOpt
				error("*** optionₘ.σ_2_Ψm⍰ ==:UniqueRelationship THAN Ψm does not need to be optimised ***")
			
			elseif optionₘ.HydroModel⍰ == :Kosugi && (optionₘ.σ_2_Ψm⍰ ==:Constrained && "Ψm" ∉ optim.ParamOpt) 
				error("*** optionₘ.σ_2_Ψm⍰ ==:Constrained THAN Ψm needs to be optimised ***")

			elseif  (optionₘ.θrOpt⍰==:σ_2_θr) && ("θr" ∈ optim.ParamOpt)
				error("*** optionₘ.θrOpt⍰==:σ_2_θr THAN θr does not need to be optimized ***") 

			elseif (optionₘ.θrOpt⍰==:σ_2_θr) && ("σ" ∉ optim.ParamOpt)
				error("*** optionₘ.θrOpt⍰==:σ_2_θr THAN σ needs to be optimized ***")
				
			elseif  (optionₘ.θrOpt⍰==:ParamPsd) && ("θr"∈ optim.ParamOpt) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt⍰==:ParamPsd && θr does not need to be optimized ***")

			elseif  (optionₘ.θrOpt⍰==:ParamPsd) && ("θr"∈ optim.ParamOpt) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt⍰==:ParamPsd && θr does not need to be optimized ***")

			elseif  (optionₘ.θrOpt⍰==:ParamPsd) && ("θr"∉ optim.ParamOpt) && !(option.run.IntergranularMixingPsd) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt⍰==:ParamPsd THAN option.run.IntergranularMixingPsd=true ***")

        	elseif option.data.SimulationKosugiθΨK && "Ks" ∉ optim.ParamOpt
            	error("***  Ks  ∉ optim.ParamOpt && option.smap.SimulationKosugiθΨK THAN Ks ∈ optim.ParamOpt***")

			# ------------ CHECKING Infiltration model--------------------
			elseif option.run.Infilt && !(option.data.Infilt)
				error("***  option.run.Infilt ⇒option.data.Infilt ***")

			# ------------ CHECKING Particle Size Distribution model--------------------
			elseif option.run.IntergranularMixingPsd && !(option.data.Psd)
				error("***  option.run.IntergranularMixingPsd ⇒option.data.Psd ***")

			elseif option.run.IntergranularMixingPsd && "Ks" ∈ optim.ParamOpt
					error("*** option.run.IntergranularMixingPsd ⇒ Ks ∉ optim.ParamOpt ***")
			
				# ------------ CHECKING Smap--------------------
			elseif option.run.Smap && option.data.Pedological⍰≠:Smap
				error("*** If option.run.Smap ⇒option.data.Pedological==Smap ***")

			elseif option.rockFragment.CorectStoneRockWetability && !(option.run.Smap)
				@warn("*** If option.data.RockWetability ⇒ option.run.Smap ***")
	

			# ------------ CHECKING Smap_2_Hypix---------------------

			elseif option.run.Smap2Hypix && (option.data.Pedological⍰ ≠ :Smap)
				error("*** option.run.Smap2Hypix ⇒option.data.Pedological⍰ == :Smap ***")
		
			end # Check error

    	return nothing
    	end  # function: CHECKING
   
end  # module: checking
# ............................................................