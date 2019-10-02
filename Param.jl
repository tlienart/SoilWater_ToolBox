

module param
	DIR_Working = pwd() # Saving the current working directory
	push!(LOAD_PATH, DIR_Working) # All the subroutines are in the current working directory
	include("Cst.jl")

    # Kosugi params minimum and maximum

    i_Sample_Start = 1
    i_Sample_End = 15 #15

    Header = ["Id", "Se_Ini",	"THETAs",	"THETAr",	"Ks", "N",	"Hvg",	"Km",	"Sigma",	"Hkg",	"Err_Ks_BestG_Kg", "Err_Ks_BestG_Vg", "Err_Ks_BestGi_Kg","Err_Ks_BestGi_Vg", "Err_Ks_Qe_Kg", "Err_Ks_Qe_Vg", "Err_Qe", "Err_Sorpt_BestG_Kg", "Err_Sorpt_BestG_Vg", "Err_Sorpt_BestGi_Kg", "Err_Sorpt_BestGi_Vg", "Err_Sorpt_Qe_Kg", "Err_Sorpt_Qe_Vg", "Hkg_Qe",	"Hvg_BestG",	"Hvg_BestGi",	"Hvg_Qe",	"Hkg_BestG",	"Kr_THETAini_Vg",	"Ks_BestG_Kg",	"Ks_BestG_Vg",	"Ks_BestGi_Kg",	"Ks_BestGi_Vg",	"Ks_Qe_Kg",	"Ks_Qe_Vg", "N_BestG",	"N_BestGi",	"N_Qe",	"NSE_Hydro_Inf_BestG_Kg",	"NSE_Hydro_Inf_BestG_Vg",	"NSE_Hydro_Inf_BestG_Vg",	"NSE_Hydro_Inf_BestGi_Kg",	"NSE_Hydro_Inf_BestGi_Vg",	"NSE_Hydro_Inf_Qe_Kg",	"NSE_Hydro_Inf_Qe_Vg",	"NSE_Inf_BestG_Vg",	"NSE_Inf_BestGi_Kg",	"NSE_Inf_BestGi_Kg",	"NSE_Inf_BestGi_Vg",	"NSE_Inf_Qe_Kg",	"NSE_Inf_Qe_Vg",	"NSE_BestGi_Kg","NSE_BestGi_Vg", "NSE_BestG_Kg","NSE_BestG_Vg", "NSE_Qe_Kg",	"NSE_Qe_Vg",	"Sorpt_BestG_Kg",	"Sorpt_BestGi_Kg",	"Sorpt_BestGi_Vg",	"Sorpt_Hydro_Kg",	"Sorpt_Hydro_Vg",	"Sorpt_Qe_Kg",	"Time_TransStead_BestG_Vg",	"Time_TransStead_Hydro_Kg", "Time_TransStead_Hydro_Vg", "Hkg_BestGi", "NSE_Hydro_Inf_BestGi_Vg",	"NSE_Inf_BestG_Kg",	"Sorpt_BestG_Kg",	"Sigma_BestG",	"Sigma_BestGi", "Sigma_Qe"]

    ΔInf_SteadyTransit = 0.07 #0.05 #[mm / T] Maximum error of not meeting the slope	

    # 45
   TransStead_Multiply = 10000

    σ_Min = 0.7
    σ_Max = 4.0

    B_Min = (2.0- cst.β) / 3.
    B_Max = B_Min + (1.0+ cst.β) /3.

    Sorptivity_Min = 0.0001 # [mm s-0.5]
    Sorptivity_Max = 5.0# [mm s-0.5]

    Ks_Min = 0.00002 #[mm s-1]
    Ks_Max = 0.5 #[mm.s-1]

    Ks_Log_Min = log10(Ks_Min) #[mm s-1]
    Ks_Log_Max = log10(Ks_Max) #[mm.s-1]

    Kr_θini_Min = 0.0001
    Kr_θini_Max = 0.9999

    Hkg_Min = 30.0#[mm]
    Hkg_Max =  300000.0#10.0 ^7 #[mm]

    Hvg_Min= 20.0#[mm]
    Hvg_Max =10^7 #[mm]
    
    N_Km1_Min = 1.01
    N_Km1_Max = 2.8

    N_Km2_Min = 2.001
    N_Km2_Max = 3.0
    
    P_σ1 = 0.38188577734776497 #0.5715138086872856 # 0.5920 0.38188577734776497
    P_σ2 = 1.0804608847580066 #0.84660977230930 # 0.7679 1.0804608847580066

    H_Min_Kg = 0.001
    H_Min_Vg = 6.7

    ϵ = 0.000001

    ZinfMax_Min = 1.0#mm
    ZinfMax_Max = 50000.0#mm

    # Output of time series of hydraulic params
    ΔSe_Inf_Max = 0.15 # The largest timestep 
    Se_Inf_Max = 0.95 # The largest θ_inf which is worth computing the hydraulic params
end
