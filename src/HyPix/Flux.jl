module flux
	import ..kunsat:Ψ_2_KUNSAT 
	export Q!, K_AVER!


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : GRAVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function GRAVITY(discret, iZ, option, param, ψ_, ψ▲)

			∂ψ∂Z = abs(ψ_ - ψ▲) / discret.ΔZ_Aver[iZ]

			if (iZ ≥ 2) && (ψ_ ≤ param.hyPix.ψϱ) && (ψ▲ ≤ param.hyPix.ψϱ) && (∂ψ∂Z > -1.0) && option.hyPix.GravityCorection
				return 0.0
			else
				return  1.0
			end
		end  # function: GRAVITY
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : K_AVER!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function K_AVER!(option, param, discret, hydro, iZ::Int64, NiZ::Int64, ψ_, ψ▲)
			if iZ == 1 # <>=<>=<>=<>=<>
				if option.hyPix.TopBoundary⍰  == "Flux" # <> = <> = <> = <> = <>
					return 0.0

				elseif option.hyPix.TopBoundary⍰  == "Ψ" # <> = <> = <> = <> = <>
					return Ψ_2_KUNSAT(option.hyPix, ψ_, 1, hydro)

				else 	# <> = <> = <> = <> = <>
					error("K_AVER! option.hyPix.TopBoundary⍰ not found")

				end

			elseif 2 ≤ iZ ≤ NiZ # <>=<>=<>=<>=<>
				return discret.ΔZ_W[iZ] * Ψ_2_KUNSAT(option.hyPix, ψ_, iZ, hydro) + (1.0 - discret.ΔZ_W[iZ]) * Ψ_2_KUNSAT(option.hyPix, ψ▲, iZ-1, hydro)

			else
				return Ψ_2_KUNSAT(option.hyPix, ψ_, NiZ, hydro)

			end
		end  # function: K_AVER!
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Q
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Q!(option, discret, hydro, iZ::Int64, iT::Int64, NiZ::Int64, param, ΔHpond, ΔPr, ΔT, θ, ψ_, ψ▲)
			if iZ == 1  # <>=<>=<>=<>=<>
				if option.hyPix.TopBoundary⍰  == "Flux" 
					return (ΔPr[iT] + ΔHpond[iT-1] - ΔHpond[iT]) / ΔT[iT]

				elseif option.hyPix.TopBoundary⍰ == "Ψ" 
					K_Aver = K_AVER!(option, param, discret, hydro, iZ, NiZ, ψ_, ψ▲)

					return K_Aver * (((ψ_ - param.hyPix.Ψ_Top) / discret.ΔZ_⬓[1]) + param.hyPix.Cosα)

				else
					error("Q! option.hyPix.TopBoundary⍰ not found")
				end

			elseif 2 ≤ iZ ≤ NiZ # <>=<>=<>=<>=<>
				K_Aver = K_AVER!(option, param, discret, hydro, iZ, NiZ, ψ_, ψ▲)

				return K_Aver * ( ((ψ_ - ψ▲) / discret.ΔZ_Aver[iZ]) + param.hyPix.Cosα * GRAVITY(discret, iZ, option, param, ψ_, ψ▲))

			else
				if option.hyPix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
					K_Aver = K_AVER!(option, param, discret, hydro, iZ, NiZ, ψ_, ψ▲)

					return K_Aver * param.hyPix.Cosα

				elseif option.hyPix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
					K_Aver = K_AVER!(option, param, discret, hydro, iZ, NiZ, ψ_, ψ▲)

					return K_Aver * ( ((param.hyPix.Ψ_Botom - ψ_) / discret.ΔZ_⬓[NiZ]) + param.hyPix.Cosα)

				elseif option.hyPix.BottomBoundary⍰ == "Q" # <>=<>=<>=<>=<>
					if param.hyPix.Q_Botom ≥ 0.0
						return min(param.hyPix.Q_Botom,  discret.ΔZ_Aver[NiZ] * (θ[iT-1,NiZ] - hydro.θr[NiZ]) / ΔT[iT])
					else
						return max(param.hyPix.Q_Botom, -  discret.ΔZ_Aver[NiZ] * (hydro.θs[NiZ] - θ[iT-1,iZ]) / ΔT[iT])
					end

				else
					error("Q! option.hyPix.BottomBoundary⍰ not found")

				end
			end # Case

		end  # function: Q!
	#-----------------------------------------------------------------



	# =============================================================
	#		module: ∂Q∂ψ
	# 		only in use if ∂R∂Ψ_Numerical = false
	# =============================================================
	module ∂q∂Ψ
		import ..flux
		export ∂Q∂Ψ, ∂Q∂Ψ△, ∂Q▽∂Ψ, ∂Q▽∂Ψ▽

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ(∂K∂Ψ::Vector{Float64}, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, option, param, Ψ::Matrix{Float64})

				if iZ == 1 # <>=<>=<>=<>=<>
					if option.hyPix.TopBoundary⍰ == "Flux" # =<>=<>=<>=<>=<>
						return 0.0

					elseif option.hyPix.TopBoundary⍰ == "Ψ" # =<>=<>=<>=<>=<>
						K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ])

						return ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - param.hyPix.Ψ_Top) / discret.ΔZ_⬓[1] + param.hyPix.Cosα) + K_Aver / discret.ΔZ_⬓[1]
					else
						error("option.hyPix.TopBoundary⍰ not found: ∂Q∂Ψ")

					end

				else # elseif 2 ≤ iZ ≤ NiZ 	<>=<>=<>=<>=<>
					K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return discret.ΔZ_W[iZ] * ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + param.hyPix.Cosα * flux.GRAVITY(discret, iZ, option, param, Ψ[iT,iZ], Ψ[iT,iZ-1])) + K_Aver / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ△
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ△(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, option, param, Ψ)
				if iZ == 1 						# <>=<>=<>=<>=<>
					return 0.0

				else #elseif 2 ≤ iZ ≤ NiZ 	# <>=<>=<>=<>=<>
					K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] *  ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + param.hyPix.Cosα * flux.GRAVITY(discret, iZ, option,  param, Ψ[iT,iZ], Ψ[iT,iZ-1])) - K_Aver / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ△
		#-----------------------------------------------------------------


		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, option, param, Ψ)
				if iZ ≤ NiZ 	# <>=<>=<>=<>=<>
					K_Aver▽ = flux.K_AVER!(option, param, discret, hydro, iZ, NiZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + param.hyPix.Cosα * flux.GRAVITY(discret, iZ, option, param, Ψ[iT,iZ], Ψ[iT,iZ-1])) - K_Aver▽ / discret.ΔZ_Aver[iZ]	
				
				else # <>=<>=<>=<>=<>
					if option.hyPix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
						return ∂K∂Ψ[NiZ] * param.hyPix.Cosα
		
					elseif option.hyPix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
						K_Aver▽ = flux.K_AVER!(option, param, discret, hydro, iZ, NiZ, Ψ[iT,NiZ], Ψ[iT,NiZ])

						return ∂K∂Ψ[NiZ] * ((param.hyPix.Ψ_Botom - Ψ[iT,NiZ]) / discret.ΔZ_⬓[NiZ] + param.hyPix.Cosα) - K_Aver▽ /  discret.ΔZ_⬓[NiZ]

					elseif option.hyPix.BottomBoundary⍰ == "Q" # <>=<>=<>=<>=<>
						return  0.0

					else
						error(" ∂Q▽∂Ψ option.hyPix.BottomBoundary⍰ not found")
					end	
				end # if iZ
			end  # function: ∂Q▽∂Ψ
		#-----------------------------------------------------------------



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ▽
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ▽(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, NiZ::Int64, option, param, Ψ)
				if iZ ≤ NiZ-1 	# <>=<>=<>=<>=<>

					K_Aver▽ = flux.K_AVER!(option, param, discret, hydro, iZ+1, NiZ, Ψ[iT,iZ+1], Ψ[iT,iZ])

					if isnan(discret.ΔZ_W[iZ+1] * ∂K∂Ψ[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + param.hyPix.Cosα * flux.GRAVITY(discret, iZ, option, param, Ψ[iT,iZ], Ψ[iT,max(iZ-1,1)])) + K_Aver▽ / discret.ΔZ_Aver[iZ+1])

						 println("Flux===", Ψ[iT,iZ]," , " ,Ψ[iT,max(iZ-1,1)])
					end

					return discret.ΔZ_W[iZ+1] * ∂K∂Ψ[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + param.hyPix.Cosα * flux.GRAVITY(discret, iZ, option, param, Ψ[iT,iZ], Ψ[iT,max(iZ-1,1)])) + K_Aver▽ / discret.ΔZ_Aver[iZ+1]
				
				else #elseif iZ == NiZ <>=<>=<>=<>=<>
					return 0.0

				end
			end  # function: ∂Q▽∂Ψ▽

	end  # module ∂q∂ψ
	# ............................................................


end # MODULE flux