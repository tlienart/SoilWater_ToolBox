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

	@option mutable struct PARAM
		globalparam::GLOBALPARAM
		hydro::HYDROS
		ksModel::KSMODEL
		psd::PSDS
		infilt::INFILTS
		smap::SMAPS 
	end

	function PARAM(Path_Data, SiteName)
		 # PARSING TOML FILE
		 Path = Path_Data * "/ParamOptionPath/" * SiteName * "_Param.toml"
		 return param = Configurations.from_toml(PARAM, Path)
	end # param
end # module param