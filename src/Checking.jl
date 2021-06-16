# =============================================================
#		module: checking
# =============================================================
module checking

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : CHECKING
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function CHECKING(option, optionₘ, optim)
        # CHECKING FOR UNCONSISTENCY WITH OPTIONS

			if option.run.HydroLabθΨ ≠ :No && !option.data.θΨ
				error("*** If option.run.HydroLabθΨ => option.data.θΨ ***")

			elseif option.run.RockCorection && option.rockFragment.RockInjectedIncluded==:Injected && !( option.data.BulkDensity && option.data.θΨ)
				error("*** If option.run.RockCorrection && option.rockFragment.RockInjectedIncluded==:Injected =>  option.data.BulkDensity OR option.data.θΨ ***")

			elseif option.run.ρᵦ_2_Φ && !(option.data.BulkDensity)
				error("*** option.run.ρᵦ_2_Φ => option.data.BulkDensity ***")

			elseif optionₘ.HydroModel==:Kosugi && "θsMacMat" ∈ optim.ParamOpt
				error("*** optionₘ.HydroModel==:Kosugi && optionₘ.HydroModel==:Bimodal THAN optionₘ.HydroModel == :Φ ***")
					
			elseif optionₘ.θrOpt≠:Opt && "θr" ∈ optim.ParamOpt
				error("*** optionₘ.θrOpt≠:Opt && θr ∈ ParamOpt THAN do not optimize θr ***")
			              
			elseif optionₘ.σ_2_Ψm ==:UniqueRelationship && "Ψm" ∈ optim.ParamOpt
				error("*** optionₘ.σ_2_Ψm ==:UniqueRelationship THAN Ψm does not need to be optimised ***")
			
			elseif optionₘ.HydroModel == :Kosugi && (optionₘ.σ_2_Ψm ==:Constrained && "Ψm" ∉ optim.ParamOpt) 
				error("*** optionₘ.σ_2_Ψm ==:Constrained THAN Ψm needs to be optimised ***")

			elseif optionₘ.KsOpt == :Data && "Ks" ∈ optim.ParamOpt
				error("*** optionₘ.KsOpt == :Data THAN Ks does not need to be optimized ***")

			elseif  (optionₘ.θrOpt==:σ_2_θr) && ("θr" ∈ optim.ParamOpt)
				error("*** optionₘ.θrOpt==:σ_2_θr THAN θr does not need to be optimized ***") 

			elseif (optionₘ.θrOpt==:σ_2_θr) && ("σ" ∉ optim.ParamOpt)
				error("*** optionₘ.θrOpt==:σ_2_θr THAN σ needs to be optimized ***")
				
			elseif  (optionₘ.θrOpt==:ParamPsd) && ("θr"∈ optim.ParamOpt) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt==:ParamPsd && θr does not need to be optimized ***")

			elseif  (optionₘ.θrOpt==:ParamPsd) && ("θr"∈ optim.ParamOpt) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt==:ParamPsd && θr does not need to be optimized ***")

			elseif  (optionₘ.θrOpt==:ParamPsd) && ("θr"∉ optim.ParamOpt) && !(option.run.IntergranularMixingPsd) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt==:ParamPsd THAN option.run.IntergranularMixingPsd=true ***")
			
			elseif ("Ks" ∈ optim.ParamOpt) && (optionₘ.KunsatΨ==false)
				error("***  (Ks ∈ optim.ParamOpt) && (KunsatΨ==false) THAN option.KunsatΨ=true ***")

        	elseif option.smap.UsePointKosugiBimodal && optionₘ.KunsatΨ && "Ks" ∉ optim.ParamOpt
            	error("***  Ks  ∉ optim.ParamOpt && option.smap.UsePointKosugiBimodal && optionₘ.KunsatΨ THAN Ks ∈ optim.ParamOpt***")
			
			elseif  optionₘ.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) && optionₘ.KunsatΨ && option.dataFrom.Smap
				error("optionₘ.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) THEN optionₘ.KunsatΨ=false OR UsePointKosugiBimodal = true")

			elseif  optionₘ.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) &&  "Ks" ∈ optim.ParamOpt && option.Smap
				error("optionₘ.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) THEN Ks ∉ optim.ParamOpt")
			end # Check error

    	return nothing
    	end  # function: CHECKING
   
end  # module: checking
# ............................................................