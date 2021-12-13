# =============================================================
#		module: memory

# =============================================================
module memory
   export MEMORY_MULTISTEPOPTIMISATION, MEMORY
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MEMORY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MEMORY(clim, N_∑T_Climate::Int64, NiZ::Int64, obsTheta, param)

		# N_Memory = ceil(Int, N_∑T_Climate / param.hyPix.ΔT_Min) + Int(N_∑T_Climate % param.hyPix.ΔT_Min + 1)

      N_Memory = ceil(Int, N_∑T_Climate / param.hyPix.ΔT_Min) + 10
		
      ΔEvaporation = fill(0.0::Float64, N_Memory)
      ΔHpond       = fill(0.0::Float64, N_Memory)
      ΔPet         = fill(0.0::Float64, N_Memory)
      ΔPr          = fill(0.0::Float64, N_Memory)
      ΔT           = fill(0.0::Float64, N_Memory)
      ∑Pet         = fill(0.0::Float64, N_Memory)
      ∑Pr          = fill(0.0::Float64, N_Memory)
      ∑T           = fill(0.0::Float64, N_Memory)

      ΔSink = fill(0.0::Float64, N_Memory, NiZ)
      Ψ     = fill(0.0::Float64, N_Memory, NiZ)
      θ     = fill(0.0::Float64, N_Memory, NiZ)
		
      Q     = fill(0.0::Float64, N_Memory, NiZ+1)
		
      Residual = fill(0.0::Float64, NiZ)
      ΔLnΨmax  = fill(0.0::Float64, NiZ)
      Ψ_Max    = fill(0.0::Float64, NiZ)
      Ψ_Min    = fill(0.0::Float64, NiZ)
      Ψbest    = fill(0.0::Float64, NiZ)
      ∂K∂Ψ     = fill(0.0::Float64, NiZ)
      ∂R∂Ψ     = fill(0.0::Float64, NiZ)
      ∂R∂Ψ△    = fill(0.0::Float64, NiZ)
      ∂R∂Ψ▽    = fill(0.0::Float64, NiZ)
      
      Nit_Reduced                  = param.hyPix.iOpt_End - param.hyPix.iOpt_Start + 1

      iNonConverge_iOpt          = fill(0  ::Int64, Nit_Reduced)

      Laiᵀ= fill(0.0::Float64, clim.N_Climate)
		CropCoeficientᵀ = fill(0.0::Float64, clim.N_Climate)

      θSim = fill(0.0::Float64, obsTheta.Nit, NiZ)
		
		return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, iNonConverge_iOpt, Laiᵀ, Q, Residual, ΔEvaporation, ΔHpond, ΔLnΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest
	end  # function: MEMORY


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : MEMORY_STEOPT
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function MEMORY_MULTISTEPOPTIMISATION(param)
      Nit_Reduced = param.hyPix.iOpt_End - param.hyPix.iOpt_Start + 1

      Efficiency                 = fill(0.0::Float64, Nit_Reduced)
      Global_WaterBalance        = fill(0.0::Float64, Nit_Reduced)
      Global_WaterBalance_NormPr = fill(0.0::Float64, Nit_Reduced)
      NseBest                    = fill(0.0::Float64, Nit_Reduced)
      CccBest                    = fill(0.0::Float64, Nit_Reduced)
      WilmotBest                 = fill(0.0::Float64, Nit_Reduced)
      SwcRoots                   = fill(0.0::Float64, Nit_Reduced)
      WofBest                    = fill(0.0::Float64, Nit_Reduced)
      ΔRunTimeHypix              = fill(0.0::Float64, Nit_Reduced)
      ΔT_Average                 = fill(0.0::Float64, Nit_Reduced)
      ∑ΔQ_Bot                    = fill(0.0::Float64, Nit_Reduced)
      ∑∑ΔSink                    = fill(0.0::Float64, Nit_Reduced)
      
   return ∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, NseBest, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average
   end  # function: MEMORY_STEOPT
   # ------------------------------------------------------------------

end  # module: memory 

# ............................................................