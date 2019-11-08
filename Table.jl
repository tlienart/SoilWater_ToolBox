# =============================================================
#		MODULE: table
# =============================================================
module table
	using ..path
	import DelimitedFiles, CSV, DataFrames
	export θΨK, SINGLEOPT_ξ1_ξ2, HYDRAULICPARAM_Psd, INTERGRANULARMIXING

	# ===================================================
    #          Discretization
    # ===================================================
	function θΨK(Id_Select, Of, Of_θΨ, Of_Kunsat, N_SoilSelect, hydro)
		Header =  ["Id" "Thetas" "Thetar" "Ks" "Sigma" "Hm" "ThetasMat" "SigmaMac" "HmMac" "Nse" "Nse_ThetaH" "Nse_Kunsat"]

		Nse = 1 .- Of
		Nse_θΨ = 1 .- Of_θΨ
		Nse_Kunsat = 1 .- Of_Kunsat

		DelimitedFiles.writedlm(path.Table_θΨK, [Header; Id_Select hydro.θs hydro.θr hydro.Ks hydro.σ hydro.Ψm hydro.θsMat hydro.σMac hydro.ΨmMac Nse Nse_θΨ Nse_Kunsat] , ",")
	 end # Table DISCRETIZATION

	

		# =================================================================================================================
		#			Psd
		# =================================================================================================================
		

			# ===================================================
			#          Single
			# ===================================================
			function SINGLEOPT_ξ1_ξ2(θr_Psd, ξ1, ξ2, Nse_Psd, Subclay, hydro)  
				Header = ["THETAs_Obs", "THETAr_Obs", "THETAr_Psd", "SIGMA_Obs", "Hkg_Obs", "THETAs_Mac", "SIGMA_Mac", "HkgMac", "xi1", "xi2", "Subclay", "NSE_Psd"]
		
				Path_Table = path.Table * "IndividualSamples_T1_T2.csv"
				println("Path_Table $Path_Table \n")
				if isfile(Path_Table) # If the file is there than delete it before saving figure
					rm(Path_Table)
				end
				
				Data = [hydro.θsMac, hydro.θr, θr_Psd, hydro.σMat, hydro.ΨkgMat, hydro.θsMat, hydro.σMac, hydro.ΨkgMac, ξ1, ξ2, Subclay, Nse_Psd]
		
				CSV.write(Path_Table, DataFrames.DataFrame(Data); header=Header)
			end # Table write single optimization


			# ===================================================
			#          Psd Hydraulic Parameters
			# ===================================================
			function HYDRAULICPARAM_Psd(θsMat_Psd, θr_Psd_Kg, σMat_Psd, ΨkgMat_Psd, θsMac, σMac_Psd, ΨkgMac_Psd, Nse_θh_Uni_Psd, Nse_θh_Bim_Psd, ∑Psd)
			
				Header =  ["THETAs_Mat_Psd" "THETAr_Mat_Psd_Kg" "SIGMA_Mat_Psd" "HkgMat_Psd" "THETAs_Mac" "SIGMA_Mac_Psd" "Hkg_Mac_Psd" "THETAmac_Psd" "Nse_THETAh_Uni_Psd"  "Nse_THETAh_Bim_Psd"  "Nse_THETAh_Improv"  "PSD_0002_mm"  "PSD_0006_mm"  "PSD_001_mm"  "PSD_002_mm" "PSD_0063_mm"  "PSD_0125_mm"  "PSD_025_mm"  "PSD_05_mm"  "PSD_1_mm"  "PSD_2_mm"]
		
				Path_Table = path.Table * "PsdModelHydraulicParam.csv"
				println("Path_Table $Path_Table \n")
				if isfile(path.Table) # If the file is there than delete it before saving figure
					rm(path.Table)
				end
				
				DelimitedFiles.writedlm(Path_Table, [Header; θsMat_Psd θr_Psd_Kg σMat_Psd ΨkgMat_Psd θsMac σMac_Psd ΨkgMac_Psd θsMac-θsMat_Psd Nse_θh_Uni_Psd Nse_θh_Bim_Psd Nse_θh_Bim_Psd-Nse_θh_Uni_Psd ∑Psd], ",")
				
			end # Table hydraulic Param

			# ===================================================
			#          Single
			# ===================================================

			function INTERGRANULARMIXING(θr_Psd, ξ1, ξ2, Of_Psd, Subclay, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, P_∑Psd_2_ξ2_C, P_ξ2_2_ξ1_A, P_ξ2_2_ξ1_B, hydro)
				Header = ["THETAs_Obs", "THETAr_Obs", "THETAr_Psd", "SIGMA_Obs", "Hkg_Obs", "THETAs_Mac", "SIGMA_Mac", "HkgMac", "xi1", "xi2", "Subclay", "Of_Psd", "Psd_2_xi2_A", "Psd_2_xi2_B", "Psd_2_xi2_C", "xi2_2_xi1_A", "xi2_2_xi1_B"]
		
				Path_Table = path.Table * "IndividualSamples_PsdParam.csv"
				println("Path_Table $Path_Table \n")
				if isfile(Path_Table) # If the file is there than delete it before saving figure
					rm(Path_Table)
				end
				
				Nsample = length(θsMac)
				∑Psd_2_ξ2_β1 = fill(∑Psd_2_ξ2_β1, Nsample)
				∑Psd_2_ξ2_β2 = fill(∑Psd_2_ξ2_β2, Nsample)
				P_∑Psd_2_ξ2_C = fill(P_∑Psd_2_ξ2_C, Nsample)
				P_ξ2_2_ξ1_A = fill(P_ξ2_2_ξ1_A , Nsample)
				P_ξ2_2_ξ1_B = fill(P_ξ2_2_ξ1_B , Nsample)
		
				Data = [hydro.θsMac, hydro.θr, θr_Psd, hydro.σMat, hydro.ΨkgMat, hydro.θsMat, hydro.σMac, hydro.ΨkgMac, ξ1, ξ2, Subclay, Of_Psd, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, P_∑Psd_2_ξ2_C, P_ξ2_2_ξ1_A, P_ξ2_2_ξ1_B]
		
				CSV.write(Path_Table, DataFrames.DataFrame(Data); header=Header)
			end # Table write single optimization





	
end  # module table
# ............................................................