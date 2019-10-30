
function HYDROPARAM()

	println("=== PSD MODEL RUNNING ===")
	println("... \n")

		if option.Psd == true
			σ_Psd = zeros(Float64, Nsample)
			Ψkg_Psd = zeros(Float64, Nsample)
			ξ1 = zeros(Float64, Nsample)
			ξ2 = zeros(Float64, Nsample)
			Of_Psd = zeros(Float64, Nsample)

			# Computing PSD data
				if option.OptimizeKΨ # If we have  OptimizeKΨ
					psd.PSD_START(Nsample, θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Ψ_θΨ, θ_θΨ, N_θΨ, Nse_θΨ_Bim_Mean, Flag_Good ; N_Kθ=N_Kθ, K_Kθ=K_Kθ, Ψ_Kθ=Ψ_Kθ)
				else
					psd.PSD_START(Nsample, θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Ψ_θΨ, θ_θΨ, N_θΨ, Nse_θΨ_Bim_Mean, Flag_Good)
				end
		end # option.Psd

	println("Particle Size Distribution Model Finished")

	return
   
 end #START_HYDROPARAM()


 HYDROPARAM()