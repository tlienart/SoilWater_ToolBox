# =============================================================
#		module: distribution
# =============================================================
module distribution
	# using CairoMakie
	export PLOTINGS, DISTRIBUTION
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DISTRIBUTION
	#     Normal or LogNormal
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# """ DISTRIBUTION Normal or LogNormal """
		function DISTRIBUTION(X, μ, σ; Distribution⍰="Normal", Normalise=false, Invert=false)

			N =length(X)
			Y = fill(0.0::Float64, N)
			for i=1:N
				if Distribution⍰ == "LogNormal"
					Y[i] = LOGNORMAL_DISTRIBUTION(X[i], μ, σ, Normalise, Invert)

				elseif  Distribution⍰ == "Normal"
					Y[i] = NORMAL_DISTRIBUTION(X[i], μ, σ, Normalise, Invert)
				
				else
					error("DISTRIBUTION: only option LogNormal or Normal supported")
				end
			end
		return Y
		end  # function: DISTRIBUTION
		# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LOGNORMAL DISTRIBUTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function LOGNORMAL_DISTRIBUTION(X, μ, σ, Normalise, Invert)
			Y = exp( -((log(X / μ)) ^ 2.0) / (2.0 * σ ^ 2.0)) / (X * σ * √(π * 2.0))

			if Normalise || Invert
				Xmode = exp(log(μ) - σ^2)
				Y = Y / ( exp( -((log(Xmode / μ)) ^ 2.0) / (2.0 * σ ^ 2.0)) / (Xmode * σ * √(π * 2.0)))
			end

			if Invert
				Y = 1 - Y
			end
		return Y
		end  # function: LOGNORMAL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LOGNORMAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NORMAL_DISTRIBUTION(X, μ, σ, Normalise, Invert)
			Y = exp(( -(X - μ) ^ 2) / (2.0 * σ ^ 2.0)) / (σ * √(π * 2.0))
			if Normalise || Invert
				Y = Y * σ * √(π * 2.0)
			end

			if Invert
				Y = 1 - Y
			end
		return Y
		end  # function: LOGNORMAL
		# ------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function PLOTINGS()
		# 	X = range(0 , stop=4.0, length=100) 

		# 	Fig = Figure(backgroundcolor=RGBf0(0.98, 0.98, 0.98), resolution = (2500, 1000),  font="Sans", fontsize=16)

		# 	Axis1 = Axis(Fig[1,1], title="Distributions", titlesize=24, xlabel="X", ylabel="Y", xlabelsize=35,  ylabelsize=35, backgroundcolor=:white)

		# 	# Y = DISTRIBUTION(X, 2, 4.0,  Distribution⍰="Normal", Normalise =true, Invert =true)
		# 	# Fig1 = CairoMakie.lines!(Fig[1,1], X, Y )

		# 	Y = DISTRIBUTION(X, 2, log(2)/3.0,  Distribution⍰="LogNormal", Normalise =true, Invert=false)
		# 	Fig1 = CairoMakie.lines!(Fig[1,1], X, Y )

		# 	Fig[1, 1] = Axis1

		# 	Path ="D:\\Temp\\NormalDist\\Dist_1.svg"
		# 	CairoMakie.save(Path, Fig)
			
		# return nothing
		# end  # function: PLOT
	
end  # module: distribution
# ............................................................