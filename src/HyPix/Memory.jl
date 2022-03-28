# =============================================================
#		module: memory

# =============================================================
module memory
   export MEMORY_MULTISTEPOPTIMISATION, MEMORY
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MEMORY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MEMORY(clim, N_∑T_Climate::Int64, NiZ::Int64, obsTheta, paramHypix)

      N_Memory = ceil(Int, N_∑T_Climate / paramHypix.ΔT_Min) + 10
		
      ΔEvaporation = fill(0.0::Float64, N_Memory)
      Hpond        = similar(ΔEvaporation)
      ΔPet         = similar(ΔEvaporation)
      ΔPr          = similar(ΔEvaporation)
      ΔRunoff      = similar(ΔEvaporation)
      ΔT           = similar(ΔEvaporation)
      ∑Pet         = similar(ΔEvaporation)
      ∑Pr          = similar(ΔEvaporation)
      ∑T           = similar(ΔEvaporation)

      ΔSink = fill(0.0::Float64, N_Memory, NiZ)
      Ψ     = similar(ΔSink)
      θ     = similar(ΔSink)
		
      Q     = fill(0.0::Float64, N_Memory, NiZ+1)
		
      Residual = fill(0.0::Float64, NiZ)
      ΔLnΨmax  = similar(Residual)
      Ψ_Max    = similar(Residual)
      Ψ_Min    = similar(Residual)
      Ψbest    = similar(Residual)
      ∂K∂Ψ     = similar(Residual)
      ∂R∂Ψ     = similar(Residual)
      ∂R∂Ψ△    = similar(Residual)
      ∂R∂Ψ▽    = similar(Residual)
      
      Nit_Reduced                  = paramHypix.opt.iOptMultiStep_End - paramHypix.opt.iOptMultiStep_Start + 1

      iNonConverge_iOpt          = fill(0::Int64  ::Int64, Nit_Reduced)

      Laiᵀ= fill(0.0::Float64, clim.N_Climate)
		CropCoeficientᵀ = similar(Laiᵀ)

      θSim = fill(0.0::Float64, obsTheta.Nit, NiZ)
		
	return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, Hpond, iNonConverge_iOpt, Laiᵀ, Q, Residual, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPr, ΔRunoff, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest
	end  # function: MEMORY


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : MEMORY_STEOPT
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function MEMORY_MULTISTEPOPTIMISATION(paramHypix)
      Nit_Reduced = paramHypix.opt.iOptMultiStep_End - paramHypix.opt.iOptMultiStep_Start + 1

      Efficiency                 = fill(0.0::Float64, Nit_Reduced)
      Global_WaterBalance        = similar(Efficiency)
      Global_WaterBalance_NormPr = similar(Efficiency)
      NseBest                    = similar(Efficiency)
      CccBest                    = similar(Efficiency)
      WilmotBest                 = similar(Efficiency)
      SwcRoots                   = similar(Efficiency)
      WofBest                    = similar(Efficiency)
      ΔRunTimeHypix              = similar(Efficiency)
      ΔT_Average                 = similar(Efficiency)
      ∑ΔQ_Bot                    = similar(Efficiency)
      ∑∑ΔSink                    = similar(Efficiency)
      
   return ∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, NseBest, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average
   end  # function: MEMORY_STEOPT
   # ------------------------------------------------------------------

end  # module: memory 

# ............................................................