# =============================================================
#		MODULE WRC
# =============================================================
module wrc
	export ∂Se∂Ψ, ∂θ∂Ψ, ∂Ψ∂Se, ∂Ψ∂θ, Se_2_θ, θ_2_Se, θ_2_ΨDual, Ψ_2_SeDual, Ψ_2_θDual

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_Se
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θ_2_Se(θ₂, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr = hydroParam.θr[iZ])
			Se = (θ₂ - θr) / (θs - θr)
			return max( min(Se, 1.0), 0.0)
		end # function θ_2_Se
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_θ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Se_2_θ(Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ])
			θ₂ = Se * (θs - θr) + θr
			return max(min(θ₂, θs), θr)
		end # function Se_2_θ`
	#-----------------------------------------------------------------

		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_θDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_θDual(optionₘ, Ψ₁, iZ, hydroParam)

			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg. Ψ_2_θDual(optionₘ,Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return wrc.vgJules.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for Ψ_2_θDual is not yet available")
			end
		end # function Ψ_2_θDual
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_SeDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydroParam)

			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.Ψ_2_SeDual(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return wrc.vgJules.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for Ψ_2_SeDual is not yet available")
			end # function Ψ_2_θDual
		end
	#-----------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_ΨDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θ_2_ΨDual(optionₘ, θ₁, iZ, hydroParam)
			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.θ_2_ΨDual(optionₘ, θ₁, iZ, hydroParam)
			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.θ_2_Ψ(θ₁, iZ, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.θ_2_Ψ(θ₁, iZ, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.θ_2_Ψ(θ₁, iZ, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for θ_2_ΨDual is not yet available")
			end # function θ_2_ΨDual
		end
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_ΨDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Se_2_ΨDual(optionₘ, Se₁, iZ, hydroParam)
			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.Se_2_ΨDual(optionₘ, Se₁, iZ, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for Se_2_ΨDual is not yet available")
			end # function Se_2_ΨDual
		end
	#-----------------------------------------------------------------

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂θ∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			
			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂θ∂Ψ is not yet available")	
			end
		end # function ∂θ∂Ψ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂Ψ∂θ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂Ψ∂θ(optionₘ, θ₁, iZ::Int64, hydroParam)	
			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂Ψ∂θ(optionₘ, θ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂Ψ∂θ is not yet available")	
			end
		end # function ∂Ψ∂θ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂Se∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂Se∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			
			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂Se∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂Se∂Ψ is not yet available")	
			end
		end # function ∂Se∂Ψ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂Se∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂Ψ∂Se(optionₘ, Se₁, iZ::Int64, hydroParam)	
			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂Ψ∂Se(optionₘ, Se₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂Ψ∂Se is not yet available")
			end
		end # function ∂Se∂Ψ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : GREEN-AMPT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function GREEN_AMPT(optionₘ, iZ::Int64, hydroParam)
			if optionₘ.HydroModel⍰ == "BrooksCorey"
				return  wrc.bc.GREEN_AMPT(optionₘ, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.GREEN_AMPT(optionₘ, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for GREEN_AMPT is not yet available")
			end
		end # function: GREEN_AMPT



	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	
	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import SpecialFunctions: erfc, erfcinv
		import ForwardDiff, Optim
		import ..wrc
		export ∂Se∂Ψ, ∂θ∂Ψ, ∂Ψ∂Se, ∂Ψ∂θ, Se_2_ΨDual, θ_2_ΨDual, Ψ_2_SeDual, Ψ_2_θDual

		θsθsMacMat =  0.001
	
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_θDual
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θDual(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0)) + θr

				if θs - θsMacMat > θsθsMacMat
					θ_Mac = 0.5 * (θs - θsMacMat) * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0))
				else
					θ_Mac = 0.0
				end

			return θ_Mac + θ_Mat
			end # function Ψ_2_θDual
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_SeDual
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				θ₂ = Ψ_2_θDual(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

			return wrc.θ_2_Se(θ₂, iZ::Int64, hydroParam)
			end # function Ψ_2_SeDual
		#-----------------------------------------------------------------
	

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θ_2_ΨDual
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_ΨDual(optionₘ, θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				function Of(Ψ₁, iZ, hydroParam)
					θmod = Ψ_2_θDual(optionₘ,Ψ₁, iZ, hydroParam)
				return OF = (θ₂ - θmod) ^ 4.0
				end # Of

				if θs - θsMacMat > θsθsMacMat
					Se = max( min((θ₂ - θr) / (θs - θr), 1.0) , 0.0)
					Ψ₀ = Ψm * exp(erfcinv(2.0 * Se) * σ * √2.0)
				else 
					Optimization = Optim.optimize(Ψ₁ -> Of(10.0 ^ Ψ₁, iZ, hydroParam), log10(0.001), log10(1.0E10), Optim.GoldenSection())
					Ψ₀ = 10.0 ^ Optim.minimizer(Optimization)[1]
				end
			return max(Ψ₀, 0.0)
			# return min( max(Ψ₀, 0.0) , 1.0E8)
			end # θ_2_ΨDual
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_ΨDual
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_ΨDual(optionₘ, Se, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				θ₂ = wrc.Se_2_θ(Se, iZ, hydroParam)

			return  θ_2_ΨDual(optionₘ, θ₂, iZ, hydroParam)
			end # Se_2_ΨDual
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				if Ψ₁ < eps(100.0)
					Ψ₁ += eps(10.0)
				end

				# If Ψ₁ is positive than ∂θ∂Ψ_Mat should be positive
				∂θ∂Ψ_Mat = - (θsMacMat - θr) * exp( -((log(Ψ₁ / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψ₁ * σ * √(π * 2.0))

				if θs - θsMacMat > θsθsMacMat
					∂θ∂Ψ_Mac = - (θs - θsMacMat) * exp( -((log(Ψ₁ / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψ₁ * σMac * √(π * 2.0))
				else
					∂θ∂Ψ_Mac = 0.0
				end

				# if isnan(∂θ∂Ψ₀)
				# 	println(" Ψ= $Ψ₁")
				# 	error("wrc.kg.∂θ∂Ψ = NaN")
				# else
				#	return ∂θ∂Ψ₀
				# end
			return ∂θ∂Ψ_Mat + ∂θ∂Ψ_Mac
			end # function ∂θ∂Ψ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ Mode
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ_Mode(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])
				
				Ψm_Mode = exp(log(Ψm)-σ^2)

				∂θ∂Ψ_Mode = ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=Ψm_Mode, σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ]) 

			return Ψm_Mode, ∂θ∂Ψ_Mode
			end # function ∂θ∂Ψ_Mode
		#-------------------------------------------------------------------
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Ψ∂θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Ψ∂θ(optionₘ, θ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				θ₁ = max(min(θs - eps(1000.0), θ₁), θr + eps(1000.0))
					
				θ₂ = fill(0.0::Float64, 1)
				
				∂Ψ∂θ_Numerical(θ₂::Vector) = θ_2_ΨDual(optionₘ, abs(θ₂[1]), iZ, hydroParam)
				θ₂[1] = θ₁
				Func_∂Ψ∂θ_Numerical = θ₂ -> ForwardDiff.gradient(∂Ψ∂θ_Numerical, θ₂)

				∂Ψ∂θ₀ = Func_∂Ψ∂θ_Numerical(θ₂)[1]

				if isnan(∂Ψ∂θ₀)
					println( "======== wrc.kg.∂Ψ∂θ = NaN ======== ")
					return 1.0
				else
					return ∂Ψ∂θ₀
				end
			end # function ∂Ψ∂θ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Ψ∂Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Ψ∂Se(optionₘ, Se₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

				Se₁ = max(min(Se₁, 1.0 - eps(1000.0)), eps(1000.0))
				
				Se₂ = fill(0.0::Float64, 1)
				∂Ψ∂Se_Numerical(Se₂::Vector) = Se_2_ΨDual(optionₘ, abs(Se₂[1]), iZ, hydroParam)
				Se₂[1] = Se₁
				Func_∂Ψ∂Se_Numerical = Se₂ -> ForwardDiff.gradient(∂Ψ∂Se_Numerical, Se₂)

				return Func_∂Ψ∂Se_Numerical(Se₂)[1]

				# if isnan(∂Ψ∂Se₀)
				# 	println(" ======== wrc.kg.∂Ψ∂Se = NaN ====")
				# 	return 1.0
				# else
					# return ∂Ψ∂Se₀
				# end
			end # function ∂Se∂Ψ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Se∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Se∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])
					
				Ψ₁ = max(Ψ₁, eps(1000.0))
				ψ = fill(0.0::Float64, 1)

				∂Ψ∂Se_Numerical(ψ::Vector) = Ψ_2_SeDual(optionₘ, abs(ψ[1]), iZ, hydroParam)
				ψ[1] = Ψ₁
				Func_∂Ψ∂Se_Numerical = ψ -> ForwardDiff.gradient(∂Ψ∂Se_Numerical, ψ)

				∂Se∂Ψ₀ = Func_∂Ψ∂Se_Numerical(ψ)[1]

				if isnan(∂Se∂Ψ₀)
					println("======== wrc.kg.∂Se∂Ψ₀ = NaN ======== ")
					return 1.0
				else
					return ∂Se∂Ψ₀
				end
			end # function ∂Ψ∂Se

	end # module kg # ...............................................




	# ===============================================================================================
	#		MODULE VAN GENUCHTEN
	# ===============================================================================================
	module vg
		import ..wrc
		export Ψ_2_θ, Ψ_2_Se, θ_2_Ψ, ∂θ∂Ψ
	
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
				Se = (1.0 + (Ψ₁ / Ψvg) ^ N ) ^ (-M)
			return θ₂ = wrc.Se_2_θ(Se, iZ, hydroParam)
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
			return Se = (1.0 + (Ψ₁ / Ψvg) ^ N ) ^ (-M)
			end # function Ψ_2_Se


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				θ₂ = max(min(θ₂, θs), θr)

				Se = wrc.θ_2_Se(θ₂, iZ, hydroParam) 

				M = 1. - Km / N
			return Ψ₁ = Ψvg * exp(log(exp(log(Se) / -M) - 1.) / N)
			end # Se_2_Ψ

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) 
				M = 1.0 - Km/N # van Genuchten
			return M * (θs-θr) / ( Ψvg*(1-M)) * ((Ψ₁/Ψvg) .^ (N * M)) .* (1.0 + (Ψ₁/Ψvg) .^ N) .^ (-M-1.)
			end # function ∂θ∂Ψ

	end # module vg # ..................................................................................................................



	# ======================================================================================
	#		MODULE BROOKS AND COREY
	# ======================================================================================
	module bc
		import ..wrc
		export Ψ_2_θ, ∂θ∂Ψ, Ψ_2_Se, θ_2_Ψ, GREEN_AMPT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ]) # Brooks & Corey WRC
				if Ψ₁ > Ψbc
					return θ₂ = θr + (θs - θr) * (Ψ₁ / Ψbc) ^ - λbc
				else
					return θ₂ = θs
				end
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				θ₂ = Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				return Se = wrc.θ_2_Se(θ₂, iZ::Int64, hydroParam)
			end # function Ψ_2_Se


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				Se = wrc.θ_2_Se(θ₂, iZ, hydroParam)
				return Ψ₁ = Ψbc * (Se ^ -λbc)
			end # θ_2_Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ]) 
				
				return (-λbc / Ψbc) * (θs-θr) * (Ψ₁/Ψbc) .^ (-λbc - 1.)
			end # function ∂θ∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : GREEN_AMPT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			""" Parameters derived by Brooks and Corey 1964
			Brooks RH, Corey AT. 1964. Hydraulic properties of porous media. Hydrology Papers 3, Colorado State University, Fort Collins, 27 p. 3"""
			function GREEN_AMPT(optionₘ, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], A=2.0, B=3.0, C=1.0)
				return Ψga = Ψbc * (A + B *λbc) / (C + B * λbc)
			end # GREEN_AMPT

		end # module brooksCorey # ...............................................


	# ======================================================================================
	#		MODULE CLAPP AND HORNBERGER
	# ======================================================================================
	module ch
		import ..wrc
		export Ψ_2_θ, ∂θ∂Ψ, Ψ_2_Se, θ_2_Ψ, GREEN_AMPT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch: Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ]) # Clapp & Hornberger WRC
				if Ψ₁ > Ψch
					return θ₂ = θs * (Ψ₁ / Ψch) ^ - λch
				else
					return θ₂ = θs
				end
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

				θ₂ = Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

				return Se = wrc.θ_2_Se(θ₂, iZ::Int64, hydroParam)
			end # function Ψ_2_Se

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

				Se = wrc.θ_2_Se(θ₂, iZ, hydroParam)
				return Ψ₁ = Ψch * (Se ^ -λch)
			end # θ_2_Ψ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ]) 
					
				return (-λch / Ψch) * θs * (Ψ₁/Ψch) .^ (-λch - 1.)
			end # function ∂θ∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : GREEN_AMPT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			""" Parameters derived by Brooks and Corey 1964
			Brooks RH, Corey AT. 1964. Hydraulic properties of porous media. Hydrology Papers 3, Colorado State University, Fort Collins, 27 p. 3"""
			function GREEN_AMPT(optionₘ, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], A=2.0, B=3.0, C=1.0)
				return Ψga = Ψch * (A + B *λch) / (C + B * λch)
			end # GREEN_AMPT

	end # module CLAPP AND HORNBERGER # ...............................................


	
	# ===============================================================================================
	#		MODULE VAN GENUCHTEN JULES
	# ===============================================================================================
	module vgJules
		import ..wrc
		export Ψ_2_θ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		vg FUNCTION : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
				if Ψ₁ > Ψvg
					Se = ( (Ψ₁ / Ψvg) ^ N ) ^ (-M)
				else
					Se =1.0
				end
			return θ₂ = wrc.Se_2_θ(Se, iZ, hydroParam)
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
				if Ψ₁ > Ψvg
					Se = (1+ (Ψ₁ / Ψvg) ^ N ) ^ (-M)
				else
					Se = 1.0
				end
			end # function Ψ_2_See

	end # module vgJules ...............................................

end # module wrc # ...............................................
