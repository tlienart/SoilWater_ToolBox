# =============================================================
#		MODULE: plot
# =============================================================
module plot

	using ...wrc, ...kunsat, ..path, ..cst
	using PGFPlots
	export HYDROPARAM

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROPARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROPARAM(Id_Select, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, N_SoilSelect, hydro)

		Ψ_θΨ_Min = 0.1
		Ψ_θΨ_Max = 3000.0 #kPa
		θ_θΨ_Min = 0.0
		θ_θΨ_Max = 0.58


		N_Se = 100
		Ψ = 10.0 .^ range(log(10.0 ^ -2.0), stop=log(10.0*Ψ_θΨ_Max), length=N_Se)

		Ψ_Uni = zeros(Float64, N_Se)
		Ψ_Bim = zeros(Float64, N_Se)
		θ_Uni = zeros(Float64, N_Se)
		θ_Bim = zeros(Float64, N_Se)
		Kunsat_Uni = zeros(Float64, N_Se)
		Kunsat_Bim = zeros(Float64, N_Se)



		for iSoil = 1:2 #N_SoilSelect

			Path = path.Plots_θΨK * "Theta_h_" *string(Id_Select[iSoil]) * ".svg"
			println(Path)
			println(Id_Select[iSoil], " ", θ_θΨ[iSoil,1:N_θΨ[iSoil]], " ", Ψ_θΨ[iSoil,1:N_θΨ[iSoil]], " ", N_θΨ[iSoil]) 
			println(Id_Select[iSoil], " ", K_KΨ[iSoil,1:N_KΨ[iSoil]], " ", Ψ_KΨ[iSoil,1:N_KΨ[iSoil]], " ", N_KΨ[iSoil]) 
			println(Id_Select[iSoil], " ", N_SoilSelect)
			println(hydro.θs[iSoil], " ", hydro.θr[iSoil], " ", hydro.Ks[iSoil], " ", hydro.σ[iSoil], " ", hydro.Ψm[iSoil], " ", hydro.θsMat[iSoil], " ", hydro.σMac[iSoil], " ", hydro.ΨmMac[iSoil])
			

			# Simulated
			@fastmath for i in 1:N_Se
				# Unimodal
				Se = wrc.kg.Ψ_2_SeDual(Ψ[i], iSoil::Int, hydro) 
				θ_Uni[i] = wrc.Se_2_θ(Se, iSoil::Int, hydro; θs=hydro.θs[iSoil], θr=hydro.θr[iSoil])   
				Ψ_Uni[i] = Ψ[i]

				Kunsat_Uni[i] = kunsat.kg.Se_2_KUNSAT(Se, iSoil, hydro; θs=hydro.θs[iSoil], θr=hydro.θr[iSoil], Ψm=hydro.Ψm[iSoil], σ=hydro.σ[iSoil], θsMat=hydro.θsMat[iSoil], ΨmMac=hydro.ΨmMac[iSoil], σMac=hydro.σMac[iSoil], Ks=hydro.Ks[iSoil])
							
				# # Bimodal
				# H_Bim[i] = Ψ[i]

				# θ_Bim[i] = wrc.kg.Ψ_2_θdual(Ψ[i], θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)

				# Se_Bim =  wrc.se.θ_2_Se(θ_Bim[i], θsMac, θr)

				# Kunsat_Bim[i] = kunsat.kg.Se_2_KUNSAT(Se_Bim, θsMac, θr, σMat, KsMac, θsMat, σMac)
				
			end





			println("Plotting soil $iSoil")
			Plot_CharacUnsat = GroupPlot(2, 1, groupStyle = "horizontal sep = 2.5cm, vertical sep = 1.5cm")

			push!(Plot_CharacUnsat, Axis([
				Plots.Scatter(Ψ_θΨ[iSoil,1:N_θΨ[iSoil]]*cst.kPa_2_mm/10.0, θ_θΨ[iSoil,1:N_θΨ[iSoil]], mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"),
				Plots.Linear(Ψ_Uni*cst.kPa_2_mm/10.0, θ_Uni, mark="none", style="smooth, red, very thick", legendentry=L"$KG$"),
				], 
				style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$\theta \ [cm^3 \ cm^{-3}]$", xmin=Ψ_θΨ_Min*cst.kPa_2_mm/10.0, xmax=Ψ_θΨ_Max*cst.kPa_2_mm/10.0, ymin=θ_θΨ_Min, ymax=θ_θΨ_Max, xmode="log", legendPos="south west")
			)

			push!(Plot_CharacUnsat, Axis([
				Plots.Scatter(Ψ_KΨ[iSoil,1:N_KΨ[iSoil]]*cst.kPa_2_mm/10.0, K_KΨ[iSoil,1:N_KΨ[iSoil]] * 360.0, mark="square", markSize=4, onlyMarks=true, style="red, very thick", legendentry=L"$Obs$"), 
				Plots.Linear(Ψ_Uni*cst.kPa_2_mm/10., Kunsat_Uni * 360.0, mark="none", style="smooth, red, very thick", legendentry=L"$KG$"),
				], 
				style="width=8cm, height=8cm", xlabel=L"$\psi \ [cm]$", ylabel=L"$K(\psi) \ [cm \ h^{-1}]$", xmin=Ψ_θΨ_Min*60.0*cst.kPa_2_mm/10.0, xmax=Ψ_θΨ_Max/15.0*cst.kPa_2_mm/10.0, xmode="log", legendPos="south west") # el 60.0 y el 15.0 son solo para limitar los ejes y se sacaran fuera
			
			)

			save(Path, Plot_CharacUnsat)
			
		end # for iSoil
		
		return
	end  # function: HYDROPARAM

	


end  # module plot
# ............................................................