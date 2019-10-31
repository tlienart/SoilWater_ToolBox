module psdFunc
	using ..option

	export PSD_MODEL
	import BlackBoxOptim
		
	# =========================================
   	#       PSD MODELS
	# ========================================
		function PSD_MODEL(Nrpart, Psd, Rpart, Subclay, ∑Psd, θr_Psd, θs, ξ1, ξ2)
			
			if option.psd.model == "IMP"
				# Correction for the small PSD
				Psd, ∑Psd = imp.SUBCLAY_CORRECTION(∑Psd, Subclay, Nrpart) 
				θ_Rpart = imp.RPART_2_θ(θs, θr_Psd, Psd[1:Nrpart], Rpart[1:Nrpart], Nrpart, ξ1, ξ2) 	# Computing θ from Psd
				Ψ_Rpart = imp.RPART_2_ΨRPART(Rpart, Nrpart) 											# Computing ψ from Psd

			elseif option.psd.model == "Chang2019Model"
				Psd, ∑Psd = imp.SUBCLAY_CORRECTION(∑Psd, Subclay, Nrpart)  						#? Not sure if this should be put here?
				θ_Rpart = chang.RPART_2_θ(θs, Psd[1:Nrpart], Rpart[1:Nrpart], Nrpart, ξ1) 		# Computing θ from Psd
				Ψ_Rpart = chang.RPART_2_ΨRPART(Rpart, Nrpart)
				 									# Computing ψ from Psd
			end # option.psd.Chang2019Model
			
			return θ_Rpart, Ψ_Rpart
		end # function PSD_MODEL
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>	



	# =============================================================
	#		MODULE: imp
	# =============================================================
	module imp
		import ...cst, ...param
		export ∑PSD_2_ξ2, ∑PSD_2_PSD, SUBCLAY_CORRECTION, INTERGRANULARMIXING

		# =========================================
		#      Rpart -> Ψ_Rpart
		# =========================================
			function RPART_2_ΨRPART(Rpart, Nrpart) 
				Ψ_Rpart = zeros(Float64, Nrpart)
				
				# It is to be noted that
				Rpart_Max = Rpart[Nrpart]
				Rpart_Min = Rpart[1]
				
				return Ψ_Rpart =  param.psd.Ψ_Max .* ( ( (cst.Y  ./ Rpart[1:Nrpart]) .- (cst.Y ./ Rpart_Max) ) ./ ((cst.Y  ./ Rpart_Min) - (cst.Y  ./ Rpart_Max)) ) .^ param.psd.λ 
			end # function RPART_2_ΨRPART


	# =========================================
	#      INTERGRANULARMIXING MODELS
	# =========================================
		function INTERGRANULARMIXING(Rpart, ξ1, ξ2)
			return IntergranularMixing = min.(max.(ξ1 .* exp(.-(Rpart .^ .-ξ2)), 0.0), param.psd.ξ_Max)
		end # function INTERGRANULARMIXING


	# =========================================
	#      UNIVERSAL INTERGRANULARMIXING MODEL
	# =========================================
		function ∑PSD_2_ξ2(∑Psd; ∑Psd_2_ξ2_β1=param.psd.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.psd.∑Psd_2_ξ2_β2)
			return ξ2 = min(∑Psd_2_ξ2_β1 * exp(∑Psd_2_ξ2_β2 * ∑Psd), param.psd.ξ2_Max)   # ξ2 = ∑Psd_2_ξ2_β1 + (∑Psd_2_ξ2_β2 * ∑Psd) 
		end # function ∑PSD_2_ξ2

		function MODEL_ξ1(P_ξ1=param.psd.P_ξ1)
			return ξ1 = P_ξ1
		end # function MODEL_ξ1


	# =========================================
	#       Rpart -> θ
	# =========================================
		function RPART_2_θ(θs, θr_Psd, Psd, Rpart, Nrpart, ξ1, ξ2)
			θ_Rpart = zeros(Float64, Nrpart)

			# Computing the divisor
				∑θRpart = Psd[1] / (Rpart[1] ^ INTERGRANULARMIXING(Rpart[1], ξ1, ξ2))
				@simd for iRpart=2:Nrpart
					∑θRpart +=  Psd[iRpart] / (Rpart[iRpart] ^ INTERGRANULARMIXING(Rpart[iRpart], ξ1, ξ2))
				end
		
			# Computing the dividend
				θ_Rpart[1] =  Psd[1] / (Rpart[1] ^ INTERGRANULARMIXING(Rpart[1], ξ1, ξ2))
				@simd for iRpart=2:Nrpart
					θ_Rpart[iRpart] = θ_Rpart[iRpart-1] + Psd[iRpart] / (Rpart[iRpart] ^ INTERGRANULARMIXING(Rpart[iRpart], ξ1, ξ2))
				end

			# Computing θ_Rpart
				@simd for iRpart=1:Nrpart
					θ_Rpart[iRpart] =  (θs - θr_Psd) * (θ_Rpart[iRpart] / ∑θRpart) + θr_Psd
				end

			return θ_Rpart
		end # function RPART_2_θ


	# =========================================
	#          ∑PSD -> PSD
	# =========================================
		function ∑PSD_2_PSD(∑Psd, Nrpart)
			Psd = zeros(Float64, Nrpart)

			Psd[1] = ∑Psd[1]
			@simd  for iRpart =2:Nrpart
				Psd[iRpart] = ∑Psd[iRpart] - ∑Psd[iRpart-1]
			end
			return Psd
		end # function ∑PSD_2_PSD


	# =========================================
	#          Subclay -> ∑PSD, PSD
	# =========================================
		function SUBCLAY_CORRECTION(∑Psd, Subclay, Nrpart)
			# Correction for the small PSD
			# Subclay = 1.0 # no subclay correction applied
			∑Psd[1] = ∑Psd[1] * Subclay
			Psd = ∑PSD_2_PSD(∑Psd[1:Nrpart], Nrpart)
			return Psd, ∑Psd
		end # Subclay


	# =========================================
	#          PSD -> ∑PSD
	# =========================================
		function PSD_2_∑PSD(Psd, Nrpart)
			∑Psd = zeros(Float64, Nrpart)
			∑Psd[1] = Psd[1]
			@simd  for iRpart = 2:Nrpart
				∑Psd[iRpart] = ∑Psd[iRpart-1] + Psd[iRpart]
			end
			return ∑Psd
		end # function PSD_2_∑PSD

		# =========================================
		#      Rpart -> r_pore FOR NEXT!!!!     r_pore = Rpart.*(((1.0./(θs.*ρp)).* ∑psd.*(Rpart).^(-ξ)).^(3-ξ))
		# =========================================
		# =========================================
		#      r_pore -> Ψ_Rpart FOR NEXT!!!!   Ψ_Rpart = Y ./ r_pore      #Young Laplace
		# =========================================

	end  # module: imp
	# ............................................................


	# =============================================================
	#		MODULE: Chang et al., 2019
	# =============================================================
	module chang
		import ...cst, ...param

		# ==============================================
		#      Rpart -> Ψ_Rpart  from Chang et al., 2019
		# ==============================================
			function RPART_2_ΨRPART_Chang(Rpart, Nrpart) 
				Ψ_Rpart = zeros(Float64, Nrpart)
				return Ψ_Rpart =  cst.Y ./ (0.3 .* Rpart[1:Nrpart]) 
			end # function RPART_2_ΨRPART_Chang


		# ==============================================
		#        Rpart -> θ  from Chang et al., 2019
		# ==============================================
			function RPART_2_θ_Chang(θs, Psd, Rpart, Nrpart, β)
				θ_Rpart = zeros(Float64, Nrpart)
				δ = zeros(Float64, Nrpart)
				∑ = zeros(Float64, Nrpart)

				Clay = Psd[1]
				∑[1] = Clay ^ β
				@simd  for iRpart = 2:Nrpart
					δ[iRpart] = Psd[iRpart] / sum(Psd[2:Nrpart])
					∑[iRpart] = ∑[iRpart-1] + Psd[iRpart] - ((Clay ^ β) - Clay) * δ[iRpart]
				end

				θ_Rpart[1] = θs * ∑[1] 
				@simd  for iRpart = 2:Nrpart
					θ_Rpart[iRpart] = θs * ∑[iRpart] 
				end
				return θ_Rpart
			end # function RPART_2_θ_Chang
	
	end  # module chang
	# ............................................................
	
end # module psdFunc