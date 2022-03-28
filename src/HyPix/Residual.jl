# =============================================================
#		MODULE: residual
# =============================================================
module residual
	import ..flux, ..kunsat, ..ponding, ..wrc
	import ForwardDiff: derivative
	export ∂R∂Ψ_FORWARDDIFF, ∂R∂Ψ△_FORWARDDIFF, ∂R∂Ψ▽_FORWARDDIFF, ∂RESIDUAL∂Ψ, RESIDUAL, RESIDUAL_DIFF

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 		FUNCTION : RESIDUAL_DIFF DERIVATIVE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL(discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Q::Matrix{Float64}, Residual::Vector{Float64}, Hpond::Vector{Float64}, ΔPr::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})
			
			if iZ==1
				Q[iT,1] = flux.Q!(optionHypix, discret, hydro, 1, iT, NiZ, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT,1], Ψ[iT,1])
			end

			Q[iT,iZ+1] = flux.Q!(optionHypix, discret, hydro, iZ+1, iT, NiZ, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, min(iZ+1, NiZ)], Ψ[iT,iZ])

			θ[iT,iZ] = wrc.Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)

			Residual[iZ] = discret.ΔZ[iZ] * ((θ[iT,iZ] - θ[iT-1,iZ]) - hydro.So[iZ] * (Ψ[iT,iZ]- Ψ[iT-1,iZ]) * (θ[iT,iZ] / hydro.θs[iZ])) - ΔT[iT] * (Q[iT,iZ] - Q[iT,iZ+1]) + ΔSink[iT,iZ]

		return Q, Residual, θ
		end  # function: RESIDUAL
	#----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂RESIDUAL∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂RESIDUAL∂Ψ(∂K∂Ψ::Vector{Float64}, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

			Sw = hydro.So[iZ] / hydro.θs[iZ]

			∂K∂Ψ[iZ] = kunsat.∂K∂ΨMODEL(optionHypix, Ψ[iT,iZ], iZ, hydro)

			∂θ∂Ψ₁ = wrc.∂θ∂Ψ(optionHypix, Ψ[iT,iZ], iZ, hydro)

			∂Q∂Ψ₁ = flux.∂q∂Ψ.∂Q∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Ψ)

			∂Q∂Ψ△₁ = flux.∂q∂Ψ.∂Q∂Ψ△(∂K∂Ψ, discret, hydro, iT, iZ,  NiZ, optionHypix, paramHypix, Ψ)

			∂Q▽∂Ψ₁ = flux.∂q∂Ψ.∂Q▽∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ+1,  NiZ, optionHypix, paramHypix, Ψ)

			∂Q▽∂Ψ▽₁ = flux.∂q∂Ψ.∂Q▽∂Ψ▽(∂K∂Ψ, discret, hydro, iT, iZ,  NiZ, optionHypix, paramHypix, Ψ)
		
			if iZ ≥ 2
				∂R∂Ψ△ = - ΔT[iT] * ∂Q∂Ψ△₁
			else
				∂R∂Ψ△ = 0.0::Float64
			end

			∂R∂Ψ = discret.ΔZ[iZ] * (∂θ∂Ψ₁ * (1.0 - Sw * (Ψ[iT,iZ] - Ψ[iT-1,iZ]) ) - Sw * θ[iT,iZ]) - ΔT[iT] * (∂Q∂Ψ₁ - ∂Q▽∂Ψ₁)
		
			if iZ ≤ NiZ-1
				∂R∂Ψ▽ = ΔT[iT] * ∂Q▽∂Ψ▽₁
			else
				∂R∂Ψ▽ = 0.0::Float64
			end

		return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽
		end 
	#------------------------------------------------------------------- 

	# =================================================================================
	# 		AUTOMATIC DIFFERENTIATION
	# =================================================================================


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RESIDUAL_DIFF DERIVATIVE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL_DIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)

			# Q[iT,iZ] format for ForwardDiff
				Q₁ = flux.Q!(optionHypix, discret, hydro, iZ, iT, NiZ, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ_, Ψ▲)
			# Q[iT,iZ+1] format for ForwardDiff
				Q₂ = flux.Q!(optionHypix,  discret, hydro, iZ+1, iT, NiZ, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▼, Ψ_)		
			# θ[iT,iZ] format for ForwardDiff
				θ₂ = wrc.Ψ_2_θDual(optionHypix,Ψ_, iZ, hydro)

		return discret.ΔZ[iZ] * ((θ₂ - θ[iT-1,iZ]) - hydro.So[iZ] * (Ψ_ - Ψ₀) * (θ[iT,iZ] / hydro.θs[iZ])) - ΔT[iT] * (Q₁ - Q₂) + ΔSink[iT,iZ] 
		end  # function: RESIDUAL_DIFF
	#----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂R∂Ψ_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂R∂Ψ_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)	

			ψ = Ψ_

			∂R∂Ψ_Func(ψ) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, max(ψ,0.0::Float64), Ψ▼, Ψ_Max)[1]

			∂R∂Ψ_Derivative_1 = ψ -> derivative(∂R∂Ψ_Func, ψ)	

		return ∂R∂Ψ_Derivative_1(ψ)
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix,Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)

			if iZ ≤ NiZ-1
				ψ▼ = Ψ▼

				∂R∂Ψ_Func(ψ▼) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, max(ψ▼,0.0::Float64), Ψ_Max)[1]

				∂R∂Ψ_Derivative_1 = ψ▼ -> derivative(∂R∂Ψ_Func, ψ▼)			
				
				return ∂R∂Ψ_Derivative_1(ψ▼)
			else
				return 0.0::Float64
			end
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)

			if iZ ≥ 2
				ψ▲ = Ψ▲

				∂R∂Ψ_Func(ψ▲) =  RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, max(ψ▲,0.0), Ψ₀, Ψbest_,Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)[1]
				
				∂R∂Ψ_Derivative_1 = ψ▲ -> derivative(∂R∂Ψ_Func, ψ▲)			
				
				return ∂R∂Ψ_Derivative_1(ψ▲)
			else
				return 0.0::Float64
			end
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION2 : ∂∂R∂Ψ_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂∂R∂Ψ_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)	
			ψ = Ψ_

			∂R∂Ψ_Func(ψ) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, max(ψ,0.0::Float64), Ψ▼, Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ -> derivative(∂R∂Ψ_Func, ψ)	

			∂R∂Ψ_Derivative_2 = ψ -> derivative(∂R∂Ψ_Derivative_1, ψ)	
			
			∂R∂Ψ = ∂R∂Ψ_Derivative_1(ψ)

			∂∂R∂Ψ = ∂R∂Ψ_Derivative_2(ψ)

			return ∂R∂Ψ, ∂∂R∂Ψ
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION2 : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)
			ψ▼ = Ψ▼

			∂R∂Ψ_Func(ψ▼) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, max(ψ▼,0.0::Float64), Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ▼ -> derivative(∂R∂Ψ_Func, ψ▼)
			
			∂R∂Ψ_Derivative_2 = ψ▼ -> derivative(∂R∂Ψ_Derivative_1 , ψ▼)	

		return ∂R∂Ψ_Derivative_2(ψ▼)
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL
	#-------------------------------------------------------------------

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION2 : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)
			ψ▲ = Ψ▲

			∂R∂Ψ_Func(ψ▲) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, NiZ, optionHypix, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, max(ψ▲,0.0::Float64), Ψ₀, Ψbest_,Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ▲ -> derivative(∂R∂Ψ_Func, ψ▲)
			
			∂R∂Ψ_Derivative_2 = ψ▲ -> derivative(∂R∂Ψ_Derivative_1, ψ▲)		
			
			# ∂R∂Ψ△1 = ∂R∂Ψ_Derivative_2(ψ▲)

		return ∂R∂Ψ_Derivative_2(ψ▲)
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL
	#-------------------------------------------------------------------

end  # module: residual
# ............................................................