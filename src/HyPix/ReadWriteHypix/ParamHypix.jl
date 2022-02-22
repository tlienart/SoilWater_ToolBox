# =============================================================
#		MODULE: param
# =============================================================
module paramHypix

	using Configurations

	@option mutable struct OBSΘ
		NmaxFuncEvals::Int64
		ΔHpondMax::Float64
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
		θobs_Uncert::Float64
	end
	@option mutable struct PLOT
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
		θprofile_Time::Vector{Float64}
		θΨ_Table::Vector{Float64}
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

	@option mutable struct PARAMHYPIX
		hyPix::HYPIXS
	end

	function PARAM_HYPIX(Path)
		return Configurations.from_toml(PARAMHYPIX, Path)
	end # param
end # module param