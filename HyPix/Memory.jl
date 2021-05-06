# =============================================================
#		module: memory

# =============================================================
module memory
	import ..option, ..param
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MEMORY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MEMORY(calibr, clim, N_∑T_Climate::Int64, N_iZ::Int64)

		N_Memory = ceil(Int, N_∑T_Climate / param.hypix.ΔT_Min) + Int(N_∑T_Climate % param.hypix.ΔT_Min + 1)
		
      ΔEvaporation = fill(0.0::Float64, N_Memory)
      ΔHpond       = fill(0.0::Float64, N_Memory)
      ΔPet         = fill(0.0::Float64, N_Memory)
      ΔPr          = fill(0.0::Float64, N_Memory)
      ΔT           = fill(0.0::Float64, N_Memory)
      ∑Pet         = fill(0.0::Float64, N_Memory)
      ∑Pr          = fill(0.0::Float64, N_Memory)
      ∑T           = fill(0.0::Float64, N_Memory)

      ΔSink = fill(0.0::Float64, N_Memory, N_iZ)
      Ψ     = fill(0.0::Float64, N_Memory, N_iZ)
      θ     = fill(0.0::Float64, N_Memory, N_iZ)
		
      Q     = fill(0.0::Float64, N_Memory, N_iZ+1)
		
      Residual = fill(0.0::Float64, N_iZ)
      ΔΨmax  = fill(0.0::Float64, N_iZ)
      Ψ_Max    = fill(0.0::Float64, N_iZ)
      Ψ_Min    = fill(0.0::Float64, N_iZ)
      Ψbest    = fill(0.0::Float64, N_iZ)
      ∂K∂Ψ     = fill(0.0::Float64, N_iZ)
      ∂R∂Ψ     = fill(0.0::Float64, N_iZ)
      ∂R∂Ψ△    = fill(0.0::Float64, N_iZ)
      ∂R∂Ψ▽    = fill(0.0::Float64, N_iZ)
      
      N_∑T_Plot                  = param.hypix.iSim_End - param.hypix.iSim_Start + 1

      iNonConverge_iSim          = fill(0  ::Int64, N_∑T_Plot)
      
      Efficiency                 = fill(0.0::Float64, N_∑T_Plot)
      Global_WaterBalance        = fill(0.0::Float64, N_∑T_Plot)
      Global_WaterBalance_NormPr = fill(0.0::Float64, N_∑T_Plot)
      RmseBest                   = fill(0.0::Float64, N_∑T_Plot)
      SwcRoots                   = fill(0.0::Float64, N_∑T_Plot)
      WofBest                    = fill(0.0::Float64, N_∑T_Plot)
      ΔRunTimeHypix              = fill(0.0::Float64, N_∑T_Plot)
      ΔT_Average                 = fill(0.0::Float64, N_∑T_Plot)
      ∑ΔQ_Bot                    = fill(0.0::Float64, N_∑T_Plot)
      ∑∑ΔSink                    = fill(0.0::Float64, N_∑T_Plot)

      Laiᵀ= fill(0.0::Float64, clim.N_Climate)
		CropCoeficientᵀ = fill(0.0::Float64, clim.N_Climate)

      θSim = Array{Float64}(undef, calibr.N_iT, N_iZ)

		# INITIALIZE 
         # ∑T[1]           = 0.0
         # ∑Pr[1]          = 0.0
         # ΔPr[1]          = 0.0
         # ΔT[1]           = 0.0
         # ΔPet[1]         = 0.0
         # ΔEvaporation[1] = 0.0
         # ∑Pet[1]         = 0.0
         # ΔHpond[1]       = 0.0
         # θ[iT,1:N_iZ]    .= 0.0
		
		return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑∑ΔSink, ∑Pet, ∑Pr, ∑T, ∑ΔQ_Bot, CropCoeficientᵀ, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iSim, Laiᵀ, N_Memory, Q, Residual, RmseBest, SwcRoots, WofBest, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔRunTimeHypix, ΔSink, ΔT, ΔT_Average, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest
	end  # function: MEMORY

end  # module: memory 

# ............................................................