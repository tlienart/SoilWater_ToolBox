module plot
   include("Path.jl")
   include("Wrc.jl")
   include("Kunsat.jl")
   include("Param.jl")
   include("Option.jl")
   include("Cst.jl")

   using PGFPlots

   export  PLOTTING_θr, PLOTTING_σ_Ψkg, PLOTTING, PLOTTING_OF

   function PLOTTING_θr(θr, θr_Psd, ∑Psd)
		Path = path.Plots_Psd * "ALL\\" * "ThetaR.svg"

		# # Sorting out with∑Psd
		Array = zeros(Float64, 3, length(∑Psd))
		Array[1,:] =∑Psd
		Array[2,:] = θr_Psd
		Array[3,:] = θr
		Array = sortslices(Array, dims=2)
		Psd = Array[1,:]
		θr_Psd = Array[2,:]
		θr = Array[3,:]

		# Minimum and maximum value
		θr_Min = 0.001
		θr_Max = param.θr_Max + 0.05
		Psd_Min = minimum(Psd)
		Psd_Max = maximum(Psd)

		Plot_θr = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_θr, Axis([
			PGFPlots.Plots.Scatter(Psd, θr, style="violet, very thick", onlyMarks=true, mark="diamond", markSize=4, legendentry=L"$\theta_{r} \ _{kg}$" ),
			PGFPlots.Plots.Linear(Psd, θr_Psd, style= "smooth, cyan, very thick", mark="none", legendentry= L"$\theta_{r} \ _{psd}$"),
			],
			style="width=12cm, height=8cm", title =" ", xlabel=L"$Clay \ [g \ g^{-1}]$", ylabel=L"$\theta_{r} \ [cm^3 \ cm^{-3}]$", xmin=Psd_Min, xmax=Psd_Max, ymin=θr_Min, ymax=θr_Max, legendPos="south east"))
		
		push!(Plot_θr, Axis([
			PGFPlots.Plots.Scatter(θr, θr_Psd, style="teal, very thick", onlyMarks=true, mark="o", markSize=4),
			PGFPlots.Plots.Linear(x-> x, (0., θr_Max), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=8cm, height=8cm", title =" ", xlabel=L"$\theta_{r} \ _{kg} \ [cm^3 \ cm^{-3}]$", ylabel=L"$\theta_{r} \ _{psd} \ [cm^3 \ cm^{-3}]$", xmin=θr_Min, xmax=θr_Max, ymin=θr_Min, ymax=θr_Max))
			

			PGFPlots.save(Path, Plot_θr)
	end # POTTING_θr =========================================



	function PLOTTING_∑Psd_2_ξ2(ξ2_Obs, ξ2, ∑Psd)
		Path = path.Plots_Psd * "ALL\\" * "T2_Psd.svg"

		# Sorting out with Psd
		Array = zeros(Float64, 3, length(∑Psd))
		Array[1,:] = ∑Psd
		Array[2,:] = ξ2
		Array[3,:] = ξ2_Obs
		Array = sortslices(Array, dims=2)
		∑Psd = Array[1,:]
		ξ2 = Array[2,:]
		ξ2_Obs = Array[3,:]

		ξ2_Min = 0.0 #min(minimum(ξ2_Obs), minimum(ξ2)) 
		ξ2_Max = max(maximum(ξ2_Obs), maximum(ξ2))
		∑Psd_Min = 0.0 # minimum(∑Psd)
		∑Psd_Max = maximum(∑Psd)

		Plot_ξ = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_ξ, Axis([
			PGFPlots.Plots.Scatter(ξ2_Obs, ξ2, style="teal, very thick", onlyMarks=true, mark="o", markSize=4),
			PGFPlots.Plots.Linear(x-> x, (ξ2_Min, ξ2_Max), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=8cm, height=8cm", title =" ", xlabel=L"$\tau _2 \ Obs \ [-]$", ylabel=L"$\tau _2 \ Sim \ [-]$", xmin=ξ2_Min , xmax=ξ2_Max, ymin=ξ2_Min, ymax=ξ2_Max))

		push!(Plot_ξ, Axis([
			PGFPlots.Plots.Scatter(∑Psd,ξ2_Obs, style="violet, very thick", onlyMarks=true, mark="diamond", markSize=4, legendentry=L"$\tau _2 \ obs$"),
			PGFPlots.Plots.Linear(∑Psd,ξ2, style="smooth, cyan, very thick", mark="none", legendentry=L"$\tau _2 \ sim$"),
			],
			style="width=12cm, height=8cm", title =" ", xlabel=L"$Psd \ 0.002mm \ [-]$", ylabel=L"$\tau _2 \ Obs \ [-]$", xmin=∑Psd_Min, xmax=∑Psd_Max, ymin=ξ2_Min, ymax=ξ2_Max, legendStyle ="{at={(-0.,-0.4)}, anchor=south west, legend columns=-1}"))
			PGFPlots.save(Path, Plot_ξ)
		return
	end



	function PLOTTING_σ_Ψkg(σMat, σMat_Psd, ΨkgMat, ΨkgMat_Psd, ΔθsMacMat, ΔθsMacMat_Psd)
		Path = path.Plots_Psd * "ALL\\" * "Sigma_Hkg_" * ".svg"

		# Limits
		σ_Min = min(minimum(σMat), minimum(σMat_Psd))
		σ_Max = max(maximum(σMat_Psd), maximum(σMat))

		Ψkg_Min = min(minimum(ΨkgMat), minimum(ΨkgMat_Psd))
		Ψkg_Max = max(maximum(ΨkgMat_Psd), maximum(ΨkgMat))

		# ΔθsMacMat_Min = min(minimum(ΔθsMacMat), minimum(ΔθsMacMat_Psd))
		ΔθsMacMat_Max = max(maximum(ΔθsMacMat), maximum(ΔθsMacMat_Psd))

		Plot_σ_Ψkg = GroupPlot(5, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_σ_Ψkg, Axis([
			PGFPlots.Plots.Scatter(σMat, σMat_Psd, style="violet, very thick", onlyMarks=true, mark="o", markSize=4),
			PGFPlots.Plots.Linear(x-> x, (0., σ_Max), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=8cm, height=8cm", title =" ", xlabel=L"$\sigma \_ Mat [-]$", ylabel=L"$\sigma \_ Mat \ Psd \ [-]$", xmin=σ_Min, xmax=σ_Max, ymin=σ_Min, ymax=σ_Max))

		push!(Plot_σ_Ψkg, Axis([
			PGFPlots.Plots.Scatter(ΨkgMat/10., ΨkgMat_Psd/10., style="blue, very thick", onlyMarks=true, mark="pentagon", markSize=4),
			PGFPlots.Plots.Linear(x-> x, (Ψkg_Min/10, Ψkg_Max/10.), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=12cm, height=8cm", title =" ", xlabel=L"$Hkg \ \_ Mat \ Obs\ [cm]$", ylabel=L"$Hkg \ \_ Mat \ Psd \ [cm]$", xmin=Ψkg_Min/10., xmax=Ψkg_Max/10., xmode="log", ymin=Ψkg_Min/10., ymax=Ψkg_Max/10., ymode="log", axisEqualImage=true, legendStyle ="{at={(-0.2,-0.4)}, anchor=south west, legend columns=-1}"))

		push!(Plot_σ_Ψkg, Axis([
			PGFPlots.Plots.Scatter(σMat, ΨkgMat/10., style="blue, very thick", onlyMarks=true, mark="pentagon", markSize=4),
			# PGFPlots.Plots.Linear(x-> x, (Ψkg_Min/10, Ψkg_Max/10.), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=12cm, height=8cm", title =" ", xlabel=L"$\sigma$", ylabel=L"$\psi [cm]$", xmin=0.0, xmax=σ_Max, ymin=Ψkg_Min/10., ymax=Ψkg_Max/10., ymode="log", axisEqualImage=true, legendStyle ="{at={(-0.2,-0.4)}, anchor=south west, legend columns=-1}"))

		push!(Plot_σ_Ψkg, Axis([
			PGFPlots.Plots.Scatter(σMat_Psd, ΨkgMat_Psd/10., style="blue, very thick", onlyMarks=true, mark="pentagon", markSize=4),
			# PGFPlots.Plots.Linear(x-> x, (Ψkg_Min/10, Ψkg_Max/10.), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=12cm, height=8cm", title =" ", xlabel=L"$\sigma \ Psd$", ylabel=L"$\psi \ Psd [cm]$", xmin=0.0, xmax=σ_Max, ymin=Ψkg_Min/10., ymax=Ψkg_Max/10., ymode="log", axisEqualImage=true, legendStyle ="{at={(-0.2,-0.4)}, anchor=south west, legend columns=-1}"))

		push!(Plot_σ_Ψkg, Axis([
			PGFPlots.Plots.Scatter(ΔθsMacMat, ΔθsMacMat_Psd, style="blue, very thick", onlyMarks=true, mark="pentagon", markSize=4),
			PGFPlots.Plots.Linear(x-> x, (0, ΔθsMacMat_Max), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			style="width=12cm, height=8cm", title =" ", xlabel=L"$\Delta \theta \ s $", ylabel=L"$\Delta \theta \ s \ Psd$", xmin=0.0, xmax=ΔθsMacMat_Max, ymin=0.0, ymax=ΔθsMacMat_Max, axisEqualImage=true, legendStyle ="{at={(-0.2,-0.4)}, anchor=south west, legend columns=-1}"))

		
		save(Path, Plot_σ_Ψkg)
		return
	end # PLOTTING_σ_Ψkg .................................


	function PLOTTING_OF(Nse_T1T2, Nse_Obs, Nse_Psd)
		Path = path.Plots_Psd * "ALL\\" * "Of_Nse.svg"

		Plot_Of = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_Of, Axis([
			# PGFPlots.Plots.Scatter(Nse_Obs, Nse_T1T2, style="red, very thick", onlyMarks=true, mark="Mercedes star", markSize=5, legendentry=L"$Single$"),
			# PGFPlots.Plots.Scatter(Nse_Obs, Nse_Psd, style="blue, very thick", onlyMarks=true, mark="Mercedes star flipped", markSize=5, legendentry=L"$All$"),
			PGFPlots.Plots.Scatter(Nse_Obs, Nse_Psd, style="blue, very thick", onlyMarks=true, mark="Mercedes star flipped", markSize=5),
			PGFPlots.Plots.Linear(x-> x, (0, 1), xbins=50, style="dashdotdotted, very thick, magenta"),
			],
			
			# style="width=8cm, height=8cm", title =" ", xlabel=L"$ NSE  _{\theta (\psi)} \ [-]$", xmin=0, xmax=1, ymin=0, ymax=1, ylabel=L"$NSE \ _{\theta (\psi)_{psd}} \ [-]$", legendPos="south west"))
			style="width=8cm, height=8cm", title =" ", xlabel=L"$ NSE_{\theta (\psi)} \ [-]$", xmin=0, xmax=1, ymin=0, ymax=1, ylabel=L"$NSE_{\theta (\psi)_{psd}} \ [-]$"))
		
		save(Path, Plot_Of)
	end


	# =========================================
	#       Plotting
	# ========================================
	function PLOTTING(iSoil, Nrpart, Rpart, Psd, ∑Psd, θsMac, θr, σMat, ΨkgMat, KsMat, KsMac, θsMat, σMac, ΨkgMac, Ψ_θΨ, θ_θΨ, N_θΨ, Ψ_Rpart, θ_Rpart, ξ; N_Kθ=1, K_Kθ=zeros(Float64, 1), Ψ_Kθ=zeros(Float64,1), θsMat_Psd=0., θr_Psd=0., σMat_Psd=0., ΨkgMat_Psd=0., θsMac_Psd=0., σMac_Psd=0., ΨkgMac_Psd=0.)

		Path = path.Plots_Psd * "PSD_Charac_" * param.Name * "_" *string(iSoil) * ".svg"
		println(Path)

		N_Se = 100
		Ψ = 10.0 .^ range(log(10.0 ^ -2.0), stop=log(10.0*param.Ψ_Max), length=N_Se)

		H_Uni = zeros(Float64, N_Se)
		H_Bim = zeros(Float64, N_Se)
		H_Uni_Psd = zeros(Float64, N_Se)
		H_Bim_Psd = zeros(Float64, N_Se)
		Kunsat_Uni = zeros(Float64, N_Se)
		Kunsat_Bim = zeros(Float64, N_Se)
		θ_Uni = zeros(Float64, N_Se)
		θ_Bim = zeros(Float64, N_Se)
		θ_Uni_Psd = zeros(Float64, N_Se)
		θ_Bim_Psd = zeros(Float64, N_Se)
		Pdf =zeros(Float64, N_Se)
		R_Sim = zeros(Float64, N_Se)
		H_Kunsat_Obs = zeros(Float64, N_Kθ)
		K_Kunsat_Obs = zeros(Float64, N_Kθ)

		# Observed θΨ
		H_θh_Obs = Ψ_θΨ
		θ_θh_Obs = θ_θΨ		

		# Observed k(h)
		for iH in 1:N_Kθ
			H_Kunsat_Obs[iH] = Ψ_Kθ[iH]
			K_Kunsat_Obs[iH] = K_Kθ[iH]
		end

		# Simulated
		@fastmath for i in 1:N_Se
			# Unimodal
			Se = wrc.kg.Ψ_2_Se(Ψ[i], ΨkgMat, σMat)
			θ_Uni[i] = wrc.se.Se_2_θ(Se, θsMat, θr)
			H_Uni[i] = Ψ[i]

			Kunsat_Uni[i] = kunsat.kg.Se_2_KUNSAT(Se, θsMat, θr, σMat, KsMac, θsMat, σMat)

			if option.Psd_2_HydrauParam
				Se_Psd = wrc.kg.Ψ_2_Se(Ψ[i], ΨkgMat_Psd, σMat_Psd)
				θ_Uni_Psd[i] = wrc.se.Se_2_θ(Se_Psd, θsMat_Psd, θr_Psd)
				H_Uni_Psd[i] = Ψ[i]
			end

			# Bimodal
			H_Bim[i] = Ψ[i]

			θ_Bim[i] = wrc.kg.Ψ_2_θdual(Ψ[i], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)

			Se_Bim =  wrc.se.θ_2_Se(θ_Bim[i], θsMac, θr)

			Kunsat_Bim[i] = kunsat.kg.Se_2_KUNSAT(Se_Bim, θsMac, θr, σMat, KsMac, θsMat, σMac)

			if option.Plot_Pdf
				if H_Bim[i] > 0.001
					R_Sim[i] = cst.Y / H_Bim[i]
					Pdf[i] = wrc.kg.DθDr_Dual(cst.Y / H_Bim[i], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)
				end
			end

			if option.Psd_2_HydrauParam
				θ_Bim_Psd[i] = wrc.kg.Ψ_2_θdual(Ψ[i], θsMac_Psd, θr_Psd, ΨkgMat_Psd, σMat_Psd, θsMat_Psd, ΨkgMac_Psd, σMac_Psd)
				H_Bim_Psd[i] = Ψ[i]
			end
		end

		if option.Plot_Pdf
			Pdf = Pdf ./ maximum(Pdf)
			Pdf = (Pdf .- minimum(Pdf)) ./ (1 - minimum(Pdf))
		end

		# INJECTING θsMac
		pushfirst!(θ_θh_Obs, θsMac)
		pushfirst!(θ_Rpart, θsMac)
		pushfirst!(H_θh_Obs, 1.0)
		pushfirst!(Ψ_Rpart, 1.0)
		N_θΨ += 1

		# LIMITS
		θ_Min = 0.
		θ_Max = 0.55 #θsMac

		if Ψ_Kθ[1] > 1.0 
			H_θh_Min = min(minimum(H_θh_Obs[1:N_θΨ]), minimum(Ψ_Rpart[1:Nrpart]), minimum(Ψ_Kθ))
			H_θh_Max = max(maximum(H_θh_Obs[1:N_θΨ]), maximum(Ψ_Rpart[1:Nrpart]), maximum(Ψ_Kθ))
		else
			H_θh_Min = min(minimum(H_θh_Obs[1:N_θΨ]), minimum(Ψ_Rpart[1:Nrpart]))
			H_θh_Max = max(maximum(H_θh_Obs[1:N_θΨ]), maximum(Ψ_Rpart[1:Nrpart]))
		end
				
		K_Min = 0.
		K_Max = maximum(K_Kθ)

		Plot_CharacUnsat = GroupPlot(4, 1, groupStyle = "horizontal sep = 2.cm, vertical sep = 2.5cm")

		if option.Plot_∑Pdf
			push!(Plot_CharacUnsat, Axis([
				# Plots.Linear(R_Sim, Pdf, style="smooth, blue, very thick", mark="none", legendentry=L"$Model$"),
				Plots.Linear(Rpart, ∑Psd, mark="square", markSize=4, onlyMarks=false, style="smooth, teal, very thick"),
				], 
				style="width=8cm, height=8cm", xlabel=L"$R_{part} \ [mm]$", xmode="log", ylabel=L"$\sum PSD $")
			)
		end

		if option.Plot_Pdf
			push!(Plot_CharacUnsat, Axis([
				# Plots.Linear(R_Sim, Pdf, style="smooth, blue, very thick", mark="none", legendentry=L"$Model$"),
				Plots.Linear(Rpart, Psd, mark="square", markSize=4, onlyMarks=false, style="smooth, cyan, very thick"),
				
				], 
				style="width=8cm, height=8cm", ymin=0.0, ymax=0.5, xlabel=L"$R_{part} \ [mm]$", xmode="log", ylabel=L"$PSD $")
			)
		end

		if option.Plot_IntergranularMixing 
			push!(Plot_CharacUnsat, Axis([
				Plots.Linear(Rpart, (Rpart .^ -ξ) ./ maximum(Rpart .^ -ξ), style="smooth, violet, very thick", mark="none"),
				], 
				style="width=8cm, height=8cm", xmin=0.0005, xmax=0.5, xlabel=L"$R_{part} \ [mm]$", xmode="log", ylabel=L"$Normalised \ R_{part}^{-\xi(R_{Part})} $")
			)
		end

		# Title = string("Se(Ψ) Sigma = ", σMat)
		if option.Psd_2_HydrauParam
			push!(Plot_CharacUnsat, Axis([
				Plots.Scatter(H_θh_Obs/10., θ_θh_Obs, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),

				# Plots.Linear(H_Uni/10., θ_Uni, mark="none", style="teal, very thick", legendentry=L"$Uni$"),
				Plots.Linear(H_Uni/10., θ_Uni, mark="none", style="smooth, red, very thick", legendentry=L"$KG$"),

				#Plots.Linear(H_Bim/10., θ_Bim, mark="none", style="smooth, red, very thick", legendentry=L"$Bim$"),
				
				Plots.Scatter(Ψ_Rpart/10., θ_Rpart, onlyMarks=true, style="smooth, blue, very thick", mark="o", markSize=4, legendentry=L"$PSD \ model$"),

				Plots.Linear(H_Uni_Psd/10., θ_Uni_Psd, mark="none", style="dashed, blue, very thick", legendentry=L"$PSD \ model \ KG$"),
				
				#Plots.Linear(H_Bim_Psd/10., θ_Bim_Psd, mark="none", style="dashed, blue, very thick", legendentry=L"$Bim \ Psd$"),
				], 
				style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=H_θh_Min/10., xmax=H_θh_Max/10., xmode="log", ymin=θ_Min, ymax=θ_Max, legendPos="south west" )
			)
		else
			push!(Plot_CharacUnsat, Axis([
				Plots.Scatter(H_θh_Obs/10., θ_θh_Obs,  mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),

				#Plots.Linear(H_Uni/10., θ_Uni, mark="none", style="teal, very thick", legendentry=L"$Uni$"),
				Plots.Linear(H_Uni/10., θ_Uni, mark="none", style="teal, very thick", legendentry=L"$KG$"),

				#Plots.Linear(H_Bim/10., θ_Bim, mark="none", style="smooth, red, very thick", legendentry=L"$Bim$"),
				Plots.Linear(H_Bim/10., θ_Bim, mark="none", style="smooth, red, very thick", legendentry=L"$KG \ Bim \ model $"),
				
				#Plots.Scatter(Ψ_Rpart/10., θ_Rpart, onlyMarks=true, style="smooth, orange, very thick", mark="o", markSize=4, legendentry=L"$Psd$"),
				Plots.Scatter(Ψ_Rpart/10., θ_Rpart, onlyMarks=true, style="smooth, orange, very thick", mark="o", markSize=4, legendentry=L"$PSD \ model$"),
				
				], 
				style="width=8cm, height=8cm", xlabel=L"$\Psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=H_θh_Min/10., xmax=H_θh_Max/10., xmode="log", ymin=θ_Min, ymax=θ_Max, legendStyle = "{at={(-0.1,-0.4)}, anchor=south west, legend columns=-1}")
			)

		end  
	
		if option.OptimizeKΨ && option.Plot_Kθ
			# push!(Plot_CharacUnsat, Axis([
			# 	Plots.Scatter(Ψ_Kθ/10., K_Kθ*360., mark="square", markSize=4, onlyMarks=true, style="red, very thick"),

			# 	# Plots.Linear(H_Uni/10., Kunsat_Uni*360., mark="none", style="smooth, teal, very thick"),

			# 	Plots.Linear(H_Bim/10., Kunsat_Bim*360., mark="none", style="smooth, red, very thick"),
			# 	], 
			# 	style="width=8cm, height=8cm", title=L"$K(\Psi)$", xlabel=L"$\Psi \ [cm]$", ylabel=L"$K(\Psi) \ [cm/h]$",  xmin=H_θh_Min/10., xmax=H_θh_Max/10., ymin=K_Min*360., ymax=K_Max*360., xmode="log")
			
			# )
		end
		
		
		try
			save(Path, Plot_CharacUnsat)
		catch

			Plot_CharacUnsat2 = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.cm, vertical sep = 2.5cm")


			if option.Psd_2_HydrauParam
				push!(Plot_CharacUnsat2, Axis([
					Plots.Scatter(H_θh_Obs/10., θ_θh_Obs, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),
	
					Plots.Linear(H_Uni/10., θ_Uni, mark="none", style="teal, very thick", legendentry=L"$Uni$"),
					
					Plots.Linear(H_Bim/10., θ_Bim, mark="none", style="smooth, red, very thick", legendentry=L"$Bim$"),
					
					Plots.Scatter(Ψ_Rpart/10., θ_Rpart, onlyMarks=true, style="smooth, blue, very thick", mark="o", markSize=4, legendentry=L"$Psd$"),
	
					Plots.Linear(H_Uni_Psd/10., θ_Uni_Psd, mark="none", style="dashed, orange, very thick", legendentry=L"$Uni \ Psd$"),
					
					Plots.Linear(H_Bim_Psd/10., θ_Bim_Psd, mark="none", style="dashed, blue, very thick", legendentry=L"$Bim \ Psd$"),
					], 
					style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=H_θh_Min/10., xmax=H_θh_Max/10., xmode="log", ymin=θ_Min, ymax=θ_Max, legendStyle = "{at={(-0.1,-0.4)}, anchor=south west, legend columns=-1}")
				)
			else
				push!(Plot_CharacUnsat2, Axis([
					Plots.Scatter(H_θh_Obs/10., θ_θh_Obs,  mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),
	
					Plots.Linear(H_Uni/10., θ_Uni, mark="none", style="teal, very thick", legendentry=L"$Uni$"),
					
					Plots.Linear(H_Bim/10., θ_Bim, mark="none", style="smooth, red, very thick", legendentry=L"$Bim$"),
					
					Plots.Scatter(Ψ_Rpart/10., θ_Rpart, onlyMarks=true, style="smooth, orange, very thick", mark="o", markSize=4, legendentry=L"$Psd$"),
					], 
					style="width=8cm, height=8cm", xlabel=L"$\Psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=H_θh_Min/10., xmax=H_θh_Max/10., xmode="log", ymin=θ_Min, ymax=θ_Max, legendStyle = "{at={(-0.1,-0.4)}, anchor=south west, legend columns=-1}")
				)
	
			end  
			println("Not able to plot Unsaturated hydraulic param")
			save(Path, Plot_CharacUnsat2)
		end
	end


	function PLOTTING_MIXING_PARAM(;ξ1=11.07, ξ2=0.12)
		Path = path.Plots_Psd * "ALL\\" * "Mixing_Param.svg"
		Rpart = range(0.001, stop=0.5/2.0, length=1000)
		
		ξ22 = ξ2-ξ2*0.1
		ξ13 = ξ1-ξ1*0.1
		ξ14 = ξ1+ξ1*0.5
		ξ24 = ξ2*2
		
		ξ_1 = ξ1 .* exp.(-(Rpart .^ .-ξ2))
		ξ_2 = ξ1 .* exp.(-(Rpart .^ .-0.13))# ξ_2 = ξ1 .* exp.(-(Rpart .^ .-ξ22))
		#ξ_3 = ξ1 .* exp.(-(Rpart .^ .-0.15))# ξ_3 = ξ13 .* exp.(-(Rpart .^ .-ξ2))
		# ξ_4 = ξ14 .* exp.(-(Rpart .^ .-ξ24))

		Plot_ξ = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_ξ, Axis([
			Plots.Linear(Rpart, (Rpart .^ -ξ_1) , style="smooth, blue, very thick",  mark="none", legendentry=L"$\xi_{1} = 11.07; \ \xi_{2} = 0.12 $" ),
			Plots.Linear(Rpart, (Rpart .^ -ξ_2) , style="smooth, orange, very thick", mark="none", legendentry=L"$\xi_{1} = 11.07; \ \xi_{2} = 0.13 $"),
			#Plots.Linear(Rpart, (Rpart .^ -ξ_3) , style="smooth, violet, very thick", mark="none", legendentry=L"$\xi_{1} = 11.07; \ \xi_{2} = 0.15 $"),
			#Plots.Linear(Rpart, Rpart .^ -ξ_4, style="smooth, green, very thick", mark="none", legendentry=L"$\xi_{1} = 9.84; \ \xi_{2} = 0.37 $"),	
			],
			
			style="width=12cm, height=8cm", xlabel=L"$R_{part} \ [mm]$", xmode="log", ylabel=L"$R_{part}^{-\xi(R_{part})} $", legendPos="north east")
		)
		
		save(Path, Plot_ξ)
	end # PLOTTING_MIXING_PARAM
	

	function PLOTTING_THETA_R_MODEL(;α1=15.95, α2=2.00) # se pueden introducir los parametros optimizados directamente
		Path = path.Plots_Psd * "ALL\\" * "Theta_r_Model.svg"
		Clay = range(0.0, stop=0.8, length=100)
		θr_max = 0.25

		α12 = α1+α1*0.5
		α23 = α2+α2*0.5
		α14 = α1+α1*0.5
		α24 = α2+α2*0.5
		
		θr1 = θr_max .* (1.0 .- exp.(-α1 .* Clay .^ α2)) 
		θr3 = θr_max .* (1.0 .- exp.(-α1 .* Clay .^ α23)) 
		θr2 = θr_max .* (1.0 .- exp.(-α12 .* Clay .^ α2)) 
		θr4 = θr_max .* (1.0 .- exp.(-α14 .* Clay .^ α24)) 

		Plot_θr = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_θr, Axis([
			Plots.Linear(Clay, θr1, style="smooth, blue, very thick",  mark="none", legendentry=L"$\alpha_{1} = 15.95; \ \alpha_{2} = 2.00 $" ),
			Plots.Linear(Clay, θr3, style="smooth, orange, very thick",  mark="none", legendentry=L"$\alpha_{1} = 15.95; \ \alpha_{2} = 3.00 $" ),
			Plots.Linear(Clay, θr2, style="smooth, cyan, very thick",  mark="none", legendentry=L"$\alpha_{1} = 23.92; \ \alpha_{2} = 2.00 $" ),
			Plots.Linear(Clay, θr4, style="smooth, violet, very thick",  mark="none", legendentry=L"$\alpha_{1} = 23.92; \ \alpha_{2} = 3.00 $" ),	
			],
			
			style="width=12cm, height=8cm", xlabel=L"$Clay \ [g \ g^{-1}]$", ylabel=L"$\theta _{r \ psd} \ [cm^3 \ cm^{-3}]$", legendPos="south east")   
		)
		
		save(Path, Plot_θr)
	end # PLOTTING_THETA_R_MODEL


	function PLOTTING_R_2_PSI_MODEL() # se pueden introducir los parametros directamente
		Path = path.Plots_Psd * "ALL\\" * "R_2_psi_Model.svg"
		Rmin = 0.001 # mm
		Rmax = 1.0 # mm
		Rpart = range(Rmin, stop=Rmax, length=10000)
		ψmax = 1600.0 * cst.kPa_2_mm #kPa to mm
		
		ψ = ψmax .* (((cst.Y./Rpart).-(cst.Y./Rmax))/((cst.Y./Rmin).-(cst.Y./Rmax))) .^ param.λ       # λ = 2
		ψ2 = ψmax .* (((cst.Y./Rpart).-(cst.Y./Rmax))/((cst.Y./Rmin).-(cst.Y./Rmax))) .^ (param.λ/2)  # λ = 1
		ψ3 = (cst.Y./(0.3.*Rpart)) 

		Plot_R_2_ψ = GroupPlot(1, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_R_2_ψ, Axis([
			Plots.Linear(Rpart, ψ ./cst.kPa_2_mm , style="smooth, blue, very thick",  mark="none", legendentry=L"$\lambda = 2 $" ),
			Plots.Linear(Rpart, ψ2./cst.kPa_2_mm, style="smooth, orange, very thick",  mark="none", legendentry=L"$\lambda = 1 $" ),
			Plots.Linear(Rpart, ψ3./cst.kPa_2_mm, style="smooth, cyan, very thick",  mark="none", legendentry=L"$Y \ / \ 0.3R $" ),	
			],
			
			style="width=12cm, height=8cm", xlabel=L"$R_{part} \ [mm]$", xmode="log", ylabel=L"$\psi \ [kPa]$", legendPos="north east")
		)
		
		save(Path, Plot_R_2_ψ)
	end # PLOTTING_R_2_PSI_MODEL


	 function PLOTTING_ξ2_MODEL(ξ2_Obs, ξ2, ∑Psd; β1=0.0858935363108931, β2=0.984704822032021) # se pueden introducir los parametros directamente
	#function PLOTTING_ξ2_MODEL(; β1=0.0858935363108931, β2=0.984704822032021) 
		Path = path.Plots_Psd * "ALL\\" * "Xi2_Model.svg"
		#Rpart = range(0.01, stop=1.0, length=1000)
		Rpart = param.Dpart ./ 2
		∑psd1 = [0.47, 0.68, 0.76, 0.89, 0.97, 0.99, 1.0,  1.0, 1.0, 1.0] # finer texture
		∑psd2 = [0.12, 0.19, 0.25, 0.49, 0.94, 1.0, 1.0, 1.0, 1.0, 1.0] # medium texture
		∑psd2 = [0.12, 0.19, 0.25, 0.49, 0.94, 0.99, 1.0, 1.0, 1.0, 1.0]
		∑psd3 = [0.10, 0.12, 0.13, 0.18, 0.4, 0.57, 0.8, 0.9, 0.95, 1.0] # coarser texture

		∑psd_upto_0003 = range(0.1, stop=0.7, length=1000)
		ξ2 = β1 .* exp.(β2 .* ∑psd_upto_0003 )

		Plot_ξ2_R = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

		push!(Plot_ξ2_R, Axis([
			Plots.Linear(Rpart, ∑psd1, style="smooth, blue, very thick",  mark="square", legendentry=L"$Fine \ texture$" ),
			# Plots.Linear(Rpart, ∑psd2, style="smooth, green, very thick",  mark="square", legendentry=L"$Medium \ texture$" ),	
			Plots.Scatter(Rpart[2], ∑psd1[2], style="blue",    mark="square*", markSize=4 ),
			Plots.Linear(Rpart, ∑psd3, style="smooth, orange, very thick",  mark="square", legendentry=L"$Coarse \ texture$" ),
			# Plots.Scatter(Rpart[2], ∑psd2[2], style="green",  mark="square*", markSize=4 ),
			Plots.Scatter(Rpart[2], ∑psd3[2], style="orange",  mark="square*", markSize=4 ),
			# PGFPlots.Plots.Linear(x-> x, (0., 1.), xbins=50, style="dashdotdotted, very thick, magenta")
			],
			
			style="width=12cm, height=8cm", xlabel=L"$R_{part} \ [mm]$", xmode="log", ylabel=L"$\sum \ PSD $", legendPos="south east")
		)

		push!(Plot_ξ2_R, Axis([
			#Plots.Linear(∑psd_upto_0003, ξ2, style="smooth, cyan, very thick", mark="none"  ),
			# PGFPlots.Plots.Linear(∑Psd, ξ2, style="smooth, cyan, very thick", mark="none" ),
			PGFPlots.Plots.Scatter(∑Psd,ξ2_Obs, style="violet, very thick", onlyMarks=true, mark="diamond", markSize=4 ),
			Plots.Linear(∑psd_upto_0003, ξ2, style="smooth, cyan, very thick", mark="none"  ),
			],

			style="width=12cm, height=8cm", title =" ", xlabel=L"$\sum \ PSD \ up \ to \ 0.003 \ mm \ R_{part} $", ylabel=L"$\xi_{2} \ [-]$", xmin=0.1, xmax=0.7, ymin=0.0, ymax=0.25, legendPos="south east")
		)
		
		save(Path, Plot_ξ2_R)
	end # PLOTTING_ξ2_MODEL


	


end # module plot