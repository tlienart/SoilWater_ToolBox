module table
   import DataFrames, CSV, DelimitedFiles

   export HYDRAULICPARAM, INTERGRANULARMIXING
   include("Path.jl")

    # ===================================================
    #          Hydraulic Parameters
    # ===================================================
   function HYDRAULICPARAM(Path_Table, θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Nse_θΨ_Uni, Nse_θΨ_Bim, Nse_Kθ_Uni, Nse_Kθ_Bim)
     
      Header =  ["ThetasMat", "ThetarMat", "SigmaMat", "HkgMat", "KsMat", "ThetasMac", "SigmaMac", "HkgMac", "KsMac", "Nse_THETAh_Uni", "Nse_THETAh_Bim", "Nse_Kh_Uni", "Nse_Kh_Bim", "Nse_THETAh_Improv", "Nse_Kh_Improv", "THETAmac"]

      Path_Table = Path_Table * "LabHydraulicParam.csv"
      println("Path_Table $Path_Table")
      if isfile(Path_Table) # If the file is there than delete it before saving figure
			rm(Path_Table)
		end
      
      Data = [θsMat, θr, σMat, ΨkgMat, KsMat, θsMac, σMac, ΨkgMac, KsMac, Nse_θΨ_Uni, Nse_θΨ_Bim, Nse_Kθ_Uni, Nse_Kθ_Bim, Nse_θΨ_Bim-Nse_θΨ_Uni, Nse_Kθ_Bim-Nse_Kθ_Uni,  θsMac-θsMat]

      CSV.write(Path_Table, DataFrames.DataFrame(Data); header=Header)
   end # Table hydraulic Param



   	# ===================================================
    #          Single
	# ===================================================
	function SINGLEOPT_T1_T2(θsMac, θr, θr_Psd, σMat, ΨkgMat, θsMat, σMac, ΨkgMac, ξ1, ξ2, Nse_Psd, Subclay)
		Header = ["THETAs_Obs", "THETAr_Obs", "THETAr_Psd", "SIGMA_Obs", "Hkg_Obs", "THETAs_Mac", "SIGMA_Mac", "HkgMac", "Tau1", "Tau2", "Subclay", "NSE_Psd"]
  
		Path_Table = path.Table * "IndividualSamples_T1_T2.csv"
		println("Path_Table $Path_Table \n")
		if isfile(Path_Table) # If the file is there than delete it before saving figure
			  rm(Path_Table)
		end
		
		Data = [θsMac, θr, θr_Psd, σMat, ΨkgMat, θsMat, σMac, ΨkgMac, ξ1, ξ2, Subclay, Nse_Psd]
  
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

	function INTERGRANULARMIXING(θsMac, θr, θr_Psd, σMat, ΨkgMat, θsMat, σMac, ΨkgMac, ξ1, ξ2, Of_Psd, Subclay, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, P_∑Psd_2_ξ2_C, P_ξ2_2_ξ1_A, P_ξ2_2_ξ1_B)
		Header = ["THETAs_Obs", "THETAr_Obs", "THETAr_Psd", "SIGMA_Obs", "Hkg_Obs", "THETAs_Mac", "SIGMA_Mac", "HkgMac", "Tau1", "Tau2", "Subclay", "Of_Psd", "Psd_2_T2_A", "Psd_2_T2_B", "Psd_2_T2_C", "T2_2_T1_A", "T2_2_T1_B"]
  
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
  
		Data = [θsMac, θr, θr_Psd, σMat, ΨkgMat, θsMat, σMac, ΨkgMac, ξ1, ξ2, Subclay, Of_Psd, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2, P_∑Psd_2_ξ2_C, P_ξ2_2_ξ1_A, P_ξ2_2_ξ1_B]
  
		CSV.write(Path_Table, DataFrames.DataFrame(Data); header=Header)
	 end # Table write single optimization



	# ===================================================
    #          SIMULATION
	# ===================================================
	function SIMULATION(Nse_Psd_σ, Nse_Psd_Ψkg, Nse_Psd_θs, Nse_Psd_Mean,Nse_θΨ_Bim_Mean)
		Header =  ["Nse_Psd_Sigma" ,"Nse_Psd_Hkg", "Nse_Psd_THETAs" ,"Nse_Psd","Nse_ThetaH_Bim_Mean"]

		Nse_Psd_σ_Arr = zeros(Float64,1)
		Nse_Psd_Ψkg_Arr = zeros(Float64,1)
		Nse_Psd_θs_Arr = zeros(Float64,1)
		Nse_Psd_Mean_Arr = zeros(Float64,1)
		Nse_θΨ_Bim_Mean_Arr = zeros(Float64,1)

		Nse_Psd_σ_Arr[1] = Nse_Psd_σ
		Nse_Psd_Ψkg_Arr[1] = Nse_Psd_Ψkg
		Nse_Psd_θs_Arr[1] = Nse_Psd_θs
		Nse_Psd_Mean_Arr[1] = Nse_Psd_Mean
		Nse_θΨ_Bim_Mean_Arr[1] = Nse_θΨ_Bim_Mean

		Path_Table = path.TableSimulation
		println("Path_Table $Path_Table \n")
		if isfile(Path_Table) # If the file is there than delete it before saving figure
			  rm(Path_Table)
		end

		Data = [Nse_Psd_σ_Arr, Nse_Psd_Ψkg_Arr, Nse_Psd_θs_Arr, Nse_Psd_Mean_Arr, Nse_θΨ_Bim_Mean_Arr]
  
		CSV.write(Path_Table, DataFrames.DataFrame(Data); header=Header)
	end

end # module table