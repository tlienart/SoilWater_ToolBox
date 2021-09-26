module flux
	import ..kunsat:Ψ_2_KUNSAT 
	export Q!, K_AVER!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : K_AVER!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function K_AVER!(option, param, discret, hydro, iZ::Int64, N_iZ::Int64, ψ_, ψ▲)
			if iZ == 1 # <>=<>=<>=<>=<>
				error("K_AVER! ≠ 1")
				return K_Aver = NaN

			elseif 2 ≤ iZ ≤ N_iZ # <>=<>=<>=<>=<>
				return K_Aver = discret.ΔZ_W[iZ] * Ψ_2_KUNSAT(option.hyPix,  ψ_, iZ, hydro) + (1.0 - discret.ΔZ_W[iZ]) * Ψ_2_KUNSAT(option.hyPix, ψ▲, iZ-1, hydro)

			elseif iZ == N_iZ + 1 # <>=<>=<>=<>=<>
				if option.hyPix.BottomBoundary⍰ == "Free"
					return K_Aver = Ψ_2_KUNSAT(option.hyPix, ψ▲, N_iZ, hydro)

				elseif option.hyPix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<># <>=<>=<>=<>=<>
					return K_Aver =  ( Ψ_2_KUNSAT(option.hyPix, param.hyPix.Ψ_Botom, N_iZ, hydro) + Ψ_2_KUNSAT(option.hyPix, ψ_, N_iZ, hydro) ) / 2.0

				end
			end
		end  # function: K_AVER!
	#-------------------------------------------------------------------



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Q
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Q!(option, discret, hydro, iZ::Int64, iT::Int64, N_iZ::Int64, param, ΔHpond, ΔPr, ΔT, ψ_, ψ▲)

			if iZ == 1  # <>=<>=<>=<>=<>
				return Q = (ΔPr[iT] + ΔHpond[iT-1] - ΔHpond[iT]) / ΔT[iT]

			elseif 2 ≤ iZ ≤ N_iZ # <>=<>=<>=<>=<>
				K_Aver = K_AVER!(option, param, discret, hydro, iZ, N_iZ, ψ_, ψ▲)

				return Q = K_Aver * ( ((ψ_ - ψ▲) / discret.ΔZ_Aver[iZ]) + param.hyPix.Cosα )

			elseif iZ == N_iZ + 1 # <>=<>=<>=<>=<>
				K_Aver = K_AVER!(option, param, discret, hydro, N_iZ+1, N_iZ, ψ_, ψ▲)

				if option.hyPix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
					return Q = K_Aver * param.hyPix.Cosα

				elseif option.hyPix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
					Q = K_Aver * ((( ψ_ - param.hyPix.Ψ_Botom) / discret.ΔZ[N_iZ]) + param.hyPix.Cosα)

					return Q = K_Aver * ( ((param.hyPix.Ψ_Botom - discret.ΔZ_Aver[iZ]/2.0 - ψ_) / discret.ΔZ_Aver[iZ]) + param.hyPix.Cosα )
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
			function ∂Q∂Ψ(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, option, param, Ψ)

				if iZ == 1 # <>=<>=<>=<>=<>
					return ∂Q∂Ψ = 0.0

				else # elseif 2 ≤ iZ ≤ N_iZ 	<>=<>=<>=<>=<>
					K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, N_iZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return ∂Q∂Ψ = discret.ΔZ_W[iZ] * ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + param.hyPix.Cosα) + K_Aver / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ△
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ△(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, option, param, Ψ)
				if iZ == 1 						# <>=<>=<>=<>=<>
					return ∂Q∂Ψ△ = NaN
				else #elseif 2 ≤ iZ ≤ N_iZ 	# <>=<>=<>=<>=<>
					K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, N_iZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return ∂Q∂Ψ△ = (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] *  ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + param.hyPix.Cosα) - K_Aver / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ△
		#-----------------------------------------------------------------

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, option, param, Ψ)
				if iZ ≤ N_iZ 	# <>=<>=<>=<>=<>
					K_Aver▽ = flux.K_AVER!(option, param, discret, hydro, iZ, N_iZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

					return ∂Q▽∂Ψ = (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] *  ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + param.hyPix.Cosα) - K_Aver▽ / discret.ΔZ_Aver[iZ]	
				
				else # <>=<>=<>=<>=<>
					if option.hyPix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
						return ∂Q▽∂Ψ = ∂K∂Ψ[N_iZ] * param.hyPix.Cosα
		
					elseif option.hyPix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
						K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, N_iZ, Ψ[iT,N_iZ], Ψ[iT,N_iZ-1])

						return ∂Q▽∂Ψ = ∂K∂Ψ[N_iZ] * ((param.hyPix.Ψ_Botom - Ψ[iT,N_iZ]) / discret.ΔZ[N_iZ] + param.hyPix.Cosα) - K_Aver / discret.ΔZ[N_iZ]	
					end	
				end # if iZ
			end  # function: ∂Q▽∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ▽
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ▽(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, option, param, Ψ)
				if iZ ≤ N_iZ-1 	# <>=<>=<>=<>=<>

					K_Aver▽ = flux.K_AVER!(option, param, discret, hydro, iZ+1, N_iZ, Ψ[iT,iZ+1], Ψ[iT,iZ])

					return ∂Q▽∂Ψ▽ = discret.ΔZ_W[iZ+1] * ∂K∂Ψ[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + param.hyPix.Cosα) + K_Aver▽ / discret.ΔZ_Aver[iZ+1]
				
				else #elseif iZ == N_iZ  			# <>=<>=<>=<>=<>
					if option.hyPix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
						return ∂Q▽∂Ψ▽ = NaN

					else option.hyPix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
						K_Aver = flux.K_AVER!(option, param, discret, hydro, iZ, N_iZ, Ψ[iT,iZ], Ψ[iT,iZ-1])

						return ∂Q▽∂Ψ▽ = K_Aver / discret.ΔZ_Aver[N_iZ+1]
					end
				end # if iZ
			end  # function: ∂Q▽∂Ψ▽

		end  # module ∂q∂ψ
	# ............................................................

end # MODULE flux