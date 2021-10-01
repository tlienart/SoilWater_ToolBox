# =============================================================
#		module: memory

# =============================================================
module memory
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MEMORY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function MEMORY(clim, iOpt_Count::Int64, N_∑T_Climate::Int64, N_iZ::Int64, obsTheta, param)

		# N_Memory = ceil(Int, N_∑T_Climate / param.hyPix.ΔT_Min) + Int(N_∑T_Climate % param.hyPix.ΔT_Min + 1)

      N_Memory = ceil(Int, N_∑T_Climate / param.hyPix.ΔT_Min)
		
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
      
      N_∑Treduced                  = param.hyPix.iOpt_End - param.hyPix.iOpt_Start + 1

      iNonConverge_iOpt          = fill(0  ::Int64, N_∑Treduced)
      


      Laiᵀ= fill(0.0::Float64, clim.N_Climate)
		CropCoeficientᵀ = fill(0.0::Float64, clim.N_Climate)

      θSim = fill(0.0::Float64, obsTheta.N_iT, N_iZ)
		
		return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, iNonConverge_iOpt, Laiᵀ, Q, Residual, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest
	end  # function: MEMORY

end  # module: memory 

# ............................................................