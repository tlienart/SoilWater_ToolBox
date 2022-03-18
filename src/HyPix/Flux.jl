module flux
	import ..kunsat:Ψ_2_KUNSAT 
	export Q!, K_AVER!


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : K_AVER!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function K_AVER!(optionHypix, paramHypix, discret, hydro, iZ::Int64, NiZ::Int64, ψ_, ψ▲; Pₘₑₐₙ=1)
			if iZ == 1 # <>=<>=<>=<>=<>
				if optionHypix.TopBoundary⍰  == "Flux" # <> = <> = <> = <> = <>
					return 0.0::Float64

				elseif optionHypix.TopBoundary⍰  == "Ψ" # <> = <> = <> = <> = <>
					return max(Ψ_2_KUNSAT(optionHypix, ψ_, 1, hydro), 1.0E-14)

				else 	# <> = <> = <> = <> = <>
					error("K_AVER! optionHypix.TopBoundary⍰ not found")

				end

			elseif 2 ≤ iZ ≤ NiZ # <>=<>=<>=<>=<>
				return max((discret.ΔZ_W[iZ] * Ψ_2_KUNSAT(optionHypix, ψ_, iZ, hydro) ^ Pₘₑₐₙ + (1.0 - discret.ΔZ_W[iZ]) * Ψ_2_KUNSAT(optionHypix, ψ▲, iZ-1, hydro) ^ Pₘₑₐₙ) ^ inv(Pₘₑₐₙ), 1.0E-14)
				# return max(discret.ΔZ_W[iZ] * Ψ_2_KUNSAT(optionHypix, ψ_, iZ, hydro) + (1.0 - discret.ΔZ_W[iZ]) * Ψ_2_KUNSAT(optionHypix, ψ▲, iZ-1, hydro), 1.0E-14)

			else
				return max(Ψ_2_KUNSAT(optionHypix, ψ_, NiZ, hydro), 1.0E-14)
			end
		end  # function: K_AVER!
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Q
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Q!(optionHypix, discret, hydro, iZ::Int64, iT::Int64, NiZ::Int64, paramHypix, Hpond, ΔPr, ΔSink, ΔT, θ, ψ_, ψ▲)
			if iZ == 1  # <>=<>=<>=<>=<>
				if optionHypix.TopBoundary⍰  == "Flux" 
					return max(ΔPr[iT] + Hpond[iT-1] - Hpond[iT], 0.0::Float64) / ΔT[iT]

				elseif optionHypix.TopBoundary⍰ == "Ψ" 
					K_Aver = K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, ψ_, ψ▲)

					return K_Aver * (((ψ_ - paramHypix.Ψ_Top) / discret.ΔZ_⬓[1]) + paramHypix.Cosα)

				else
					error("Q! optionHypix.TopBoundary⍰ not found")
				end

			elseif 2 ≤ iZ ≤ NiZ # <>=<>=<>=<>=<>
				K_Aver = K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, ψ_, ψ▲)
					return K_Aver * ( ((ψ_ - ψ▲) / discret.ΔZ_Aver[iZ]) + paramHypix.Cosα)

			else
				if optionHypix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
					K_Aver = K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, ψ_, ψ▲)
						return K_Aver * paramHypix.Cosα

				elseif optionHypix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
					K_Aver = K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, ψ_, ψ▲)
						return K_Aver * ( ((paramHypix.Ψ_Botom - ψ_) / discret.ΔZ_⬓[NiZ]) + paramHypix.Cosα)

				# elseif optionHypix.BottomBoundary⍰ == "Q" # <>=<>=<>=<>=<>
				# 	return paramHypix.Q_Botom

				else
					error("Q! optionHypix.BottomBoundary⍰ not found")

				end
			end # Case

		end  # function: Q!
	#-----------------------------------------------------------------



	# =============================================================
	#		module: ∂Q∂ψ
	# 		only in use if ∂R∂Ψ_NumericalAuto" = false
	# =============================================================
	module ∂q∂Ψ
		import ..flux
		export ∂Q∂Ψ, ∂Q∂Ψ△, ∂Q▽∂Ψ, ∂Q▽∂Ψ▽

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ(∂K∂Ψ::Vector{Float64}, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Ψ::Matrix{Float64})

				if iZ == 1 # <>=<>=<>=<>=<>
					if optionHypix.TopBoundary⍰ == "Flux" # =<>=<>=<>=<>=<>
						return 0.0::Float64

					elseif optionHypix.TopBoundary⍰ == "Ψ" # =<>=<>=<>=<>=<>
						K_Aver = flux.K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ])

						return ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - paramHypix.Ψ_Top) / discret.ΔZ_⬓[1] + paramHypix.Cosα) + K_Aver / discret.ΔZ_⬓[1]
					else
						error("optionHypix.TopBoundary⍰ not found: ∂Q∂Ψ")

					end

				else # elseif 2 ≤ iZ ≤ NiZ 	<>=<>=<>=<>=<>
					K_Aver = flux.K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return discret.ΔZ_W[iZ] * ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) + K_Aver / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ△
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ△(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Ψ)
				if iZ == 1 						# <>=<>=<>=<>=<>
					return 0.0::Float64

				else #elseif 2 ≤ iZ ≤ NiZ 	# <>=<>=<>=<>=<>
					K_Aver = flux.K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] *  ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) - K_Aver / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ△
		#-----------------------------------------------------------------


		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Ψ)
				if iZ ≤ NiZ 	# <>=<>=<>=<>=<>
					K_Aver▽ = flux.K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) - K_Aver▽ / discret.ΔZ_Aver[iZ]	
				
				else # <>=<>=<>=<>=<>
					if optionHypix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
						return ∂K∂Ψ[NiZ] * paramHypix.Cosα
		
					elseif optionHypix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
						K_Aver▽ = flux.K_AVER!(optionHypix, paramHypix, discret, hydro, iZ, NiZ, Ψ[iT,NiZ], Ψ[iT,NiZ])

						return ∂K∂Ψ[NiZ] * ((paramHypix.Ψ_Botom - Ψ[iT,NiZ]) / discret.ΔZ_⬓[NiZ] + paramHypix.Cosα) - K_Aver▽ /  discret.ΔZ_⬓[NiZ]

					# elseif optionHypix.BottomBoundary⍰ == "Q" # <>=<>=<>=<>=<>
					# 	return  0.0::Float64

					else
						error(" ∂Q▽∂Ψ optionHypix.BottomBoundary⍰ not found")
					end	
				end # if iZ
			end  # function: ∂Q▽∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ▽
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ▽(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, optionHypix, paramHypix, Ψ)
				if iZ ≤ NiZ-1 	# <>=<>=<>=<>=<>

					K_Aver▽ = flux.K_AVER!(optionHypix, paramHypix, discret, hydro, iZ+1, NiZ, Ψ[iT,iZ+1], Ψ[iT,iZ])

					return discret.ΔZ_W[iZ+1] * ∂K∂Ψ[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + paramHypix.Cosα) + K_Aver▽ / discret.ΔZ_Aver[iZ+1]
				
				else #elseif iZ == NiZ <>=<>=<>=<>=<>
					return 0.0::Float64
				end
			end  # function: ∂Q▽∂Ψ▽

	end  # module ∂q∂ψ
	# ............................................................


end # MODULE flux