# =============================================================
#		MODULE: param
# =============================================================
module params

	using Configurations

	@option mutable struct GLOBALPARAM
		N_iZ_Plot_Start
		N_iZ_Plot_End
	end
	@option mutable struct KG
		Ψσ_Min
		Ψσ_Max
		Ψσ
	end
	@option mutable struct SMAPS
		Ψ_Table
	end
	@option mutable struct HYDROS
		Coeff_Φ_2_θs
		θs_θsMacMat
		ΨmacMat
		Ψ_Max
		TableComplete_θΨ
		K_Table
		kg::KG
	end

	@option mutable struct KSMODEL
		σₛₚₗᵢₜ::Float64
		WeightKsSlow::Float64
	end
	
	@option mutable struct IMP
		Ψ_Max
		λ
		ξ_Max
		ξ1
		ξ1_Min
		ξ1_Max
		ξ2_Max
		∑Psd_2_ξ2_β1
		∑Psd_2_ξ2_β1_Min
		∑Psd_2_ξ2_β1_Max
		∑Psd_2_ξ2_β2
		∑Psd_2_ξ2_β2_Min
		∑Psd_2_ξ2_β2_Max
		∑Psd_2_ξ2_Size
		Subclay
		Subclay_Min
		Subclay_Max
	end
	@option mutable struct CHANG
		ξ1
		ξ1_Min
		ξ1_Max
	end
	@option mutable struct PSDS
		Psd_2_θr_α1
		Psd_2_θr_α1_Min
		Psd_2_θr_α1_Max
		Psd_2_θr_α2
		Psd_2_θr_α2_Min
		Psd_2_θr_α2_Max
		Psd_2_θr_Size
		Ψ_Table
		imp::IMP
		chang::CHANG
	end
	@option mutable struct INFILTS
		SeIni_Output
		Npoint_Infilt
		ΔSlope_Err_SteadyState
	end

	@option mutable struct OBSΘ
		NmaxFuncEvals::Int64
		Year_Start::Int64
		Month_Start::Int64
		Day_Start::Int64
		Hour_Start::Int64
		Minute_Start::Int64
		Second_Start::Int64
		Year_End::Int64
		Month_End::Int64
		Day_End::Int64
		Hour_End::Int64
		Minute_End::Int64
		Second_End::Int64
		θobs_Uncert
	end
	@option mutable struct PLOT
		Year_Start
		Month_Start
		Day_Start
		Hour_Start
		Minute_Start
		Second_Start
		Year_End
		Month_End
		Day_End
		Hour_End
		Minute_End
		Second_End
		θΨ_Table
	end
	@option mutable struct HYPIXS
		iOpt_Start::Int64
		iOpt_End::Int64
		Year_Start::Int64
		Month_Start::Int64
		Day_Start::Int64
		Hour_Start::Int64
		Minute_Start::Int64
		Second_Start::Int64
		Year_End::Int64
		Month_End::Int64
		Day_End::Int64
		Hour_End::Int64
		Minute_End::Int64
		Second_End::Int64
		ZfineCoarse::Float64
		ΔZfine::Float64
		ΔZcoarse::Float64
		Cosα::Float64
		ΔHpondMax::Float64
		Ψ_MinMin::Float64
		Ψ_MaxMax::Float64
		Ψ_Top::Float64
		Ψ_Botom::Float64
		ΔT_Min::Float64
		ΔT_Max::Float64
		N_Iter::Int64
		ΔT_Rerun::Float64
		Δθ_Max::Float64
		NewtonStep_Min::Float64
		NewtonStep_Mean::Float64
		NewtonStep_Max::Float64
		NewtonStep_Power::Float64
		WaterBalanceResidual_Max::Float64
		ΔT_Output::Float64
		obsTheta::OBSΘ
		ploting::PLOT
	end

	@option mutable struct PARAM
		globalparam::GLOBALPARAM
		hydro::HYDROS
		ksModel::KSMODEL
		psd::PSDS
		infilt::INFILTS
		hyPix::HYPIXS
		smap::SMAPS 
	end

	function PARAM(Path_Data, SiteName)
		 # PARSING TOML FILE
		 Path = Path_Data * "/ParamOptionPath/" * SiteName * "_Param.toml"
		 return param = Configurations.from_toml(PARAM, Path)
	end # param
end # module param