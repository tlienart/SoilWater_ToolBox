# =============================================================
#		MODULE: paramHypix
# =============================================================
module paramsHypix

	using Configurations

	@option mutable struct OBSΘ
		NmaxFuncEvals::Int64
		ΔHpondMax::Float64
		θobs_Uncert::Float64
	end
	@option mutable struct PLOT
		θprofile_Time::Vector{Float64}
		θΨ_Table::Vector{Float64}
	end
	@option mutable struct PARAMHYPIX
		iOpt_Start::Int64
		iOpt_End::Int64
		ZfineCoarse::Float64
		ΔZfine::Float64
		ΔZcoarse::Float64
		Cosα::Float64
		Ψ_MinMin::Float64
		Ψ_MaxMax::Float64
		Ψ_Top::Float64
		Ψ_Botom::Float64
		Q_Botom::Float64
		ΔT_Min::Float64
		ΔT_Max::Float64
		N_Iter::Int64
		ΔT_Rerun::Float64
		Δθ_Max::Float64
		NewtonStep_Min::Float64
		NewtonStep_Max::Float64
		NewtonStep_Power::Float64
		WaterBalanceResidual_Max::Float64
		ΔT_Output::Float64
		obsTheta::OBSΘ
		ploting::PLOT
	end

	function PARAM_HYPIX(PathParamHypix::String)
		return Configurations.from_toml(PARAMHYPIX, PathParamHypix)
	end # paramHypix
end # module paramHypix