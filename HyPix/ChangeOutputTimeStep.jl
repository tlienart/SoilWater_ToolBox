# =============================================================
#		MODULE: tincrease
# =============================================================
module Δtchange
	import ..cst, ..tool, ..param, ..interpolate
	import Dates: value, DateTime, Day, Second, Hour

	function CHANGE_OUTPUT_ΔT(∑Pet, ∑Pr, ∑T, ∑WaterBalance_η, ∑ΔSink, calibr, clim, N_iT::Int64, N_iZ::Int64, Q, veg, ΔEvaporation, ΔHpond, ΔT, θ, Ψ, ∑T_Climate)

		# PREPROCESSING ∑Evaporation, ∑ΔQ
         ∑Evaporation = Vector{Float64}(undef, N_iT)
         ∑ΔQ          = Array{Float64}(undef, N_iT, N_iZ+1)

			∑Evaporation[1]    = 0.0
         ∑ΔQ[1, 1:N_iZ+1]  .= 0.0
			for iT=2:N_iT
				∑Evaporation[iT] = ∑Evaporation[iT-1] + ΔEvaporation[iT]
				for iZ=1:N_iZ+1
					∑ΔQ[iT,iZ] = ∑ΔQ[iT-1, iZ] + ΔT[iT] * Q[iT,iZ]
				end
			end

		# ∑Pr_Gross
			∑Pr_Gross = Vector{Float64}(undef, clim.N_Climate)
			∑Pr_Gross[1] = 0.0
			for iT= 2:clim.N_Climate
				∑Pr_Gross[iT] = ∑Pr_Gross[iT-1] + clim.Pr[iT]
			end

		# PREPARING DATA FOR PLOTS
			Date_Start = clim.Date[2]
					
			Date_Start_Calibr = DateTime(param.hypix.calibr.Year_Start, param.hypix.calibr.Month_Start, param.hypix.calibr.Day_Start, param.hypix.calibr.Hour_Start, param.hypix.calibr.Minute_Start, param.hypix.calibr.Second_Start)
				
			Date_End_Calibr = DateTime(param.hypix.calibr.Year_End, param.hypix.calibr.Month_End, param.hypix.calibr.Day_End, param.hypix.calibr.Hour_End, param.hypix.calibr.Minute_End, param.hypix.calibr.Second_End)

			ΔT_Sim = value(Date_End_Calibr - Date_Start_Calibr) / 1000

			∑T_Plot = collect(range(0.0, step=cst.Day_2_Second, stop=ΔT_Sim)) 
			
			# Take account that we are starting at Date_Start_Calibr
			ΔT_Start_Calibr = value(Date_Start_Calibr - Date_Start) / 1000

			@. ∑T_Plot = ∑T_Plot + ΔT_Start_Calibr

			N_∑T_Plot = length(∑T_Plot)

		# PREPARING DATES WITH INTERVAL:
			∑T_Date_Plot = range(Date_Start_Calibr, step=Second(cst.Day_2_Second), Date_End_Calibr)	

		# INTERPOLATING DATA
			θ_Plot = fill(0.0::Float64, N_∑T_Plot, N_iZ)	
			θ_Plot = interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, N_iZ, θ_Plot, θ)

			θobs_Plot = fill(0.0::Float64, N_∑T_Plot, N_iZ)	
			θobs_Plot = interpolate.INTERPOLATE_2D_LOOP(calibr.∑T[1:calibr.N_iT], ∑T_Plot, N_∑T_Plot, calibr.N_iT, calibr.Ndepth, θobs_Plot, calibr.θobs[1:calibr.N_iT,1:calibr.Ndepth])

			Ψ_Plot = fill(0.0::Float64, N_∑T_Plot, N_iZ)
			Ψ_Plot =  interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, N_iZ, Ψ_Plot, Ψ) .* cst.Mm_2_Cm

			∑Flux_Plot = fill(0.0::Float64, N_∑T_Plot, N_iZ+1)
			∑Flux_Plot = interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, N_iZ+1, ∑Flux_Plot, ∑ΔQ)

			ΔT_Plot = fill(0.0::Float64, N_∑T_Plot)
			ΔT_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ΔT_Plot, ΔT)

			∑Evaporation_Plot = fill(0.0::Float64, N_∑T_Plot)
			∑Evaporation_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ∑Evaporation_Plot, ∑Evaporation)

			∑Pet_Plot = fill(0.0::Float64, N_∑T_Plot)
			∑Pet_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ∑Pet_Plot, ∑Pet)

			∑Pr_Plot = fill(0.0::Float64, N_∑T_Plot)
			∑Pr_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ∑Pr_Plot, ∑Pr)

			∑WaterBalance_η_Plot = fill(0.0::Float64, N_∑T_Plot)
			∑WaterBalance_η_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ∑WaterBalance_η_Plot, ∑WaterBalance_η)

			∑∑Sink = fill(0.0::Float64, N_∑T_Plot)
			∑∑Sink = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ∑∑Sink, ∑ΔSink)

			∑PrGross_Plot = fill(0.0::Float64, N_∑T_Plot)	
			∑PrGross_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T_Climate[1:clim.N_Climate], ∑T_Plot, N_∑T_Plot, clim.N_Climate, ∑PrGross_Plot, ∑Pr_Gross)

			# From ∑ to Δ
            ΔEvaporation_Plot = fill(0.0::Float64, N_∑T_Plot)
            ΔFlux_Plot        = fill(0.0::Float64, N_∑T_Plot, N_iZ+1)
            ΔPet_Plot         = fill(0.0::Float64, N_∑T_Plot)
				ΔPr_Plot          = fill(0.0::Float64, N_∑T_Plot)
				ΔPrGross_Plot     = fill(0.0::Float64, N_∑T_Plot)
            ΔSink_Plot        = fill(0.0::Float64, N_∑T_Plot)
            Date_Plot         = Vector{DateTime}(undef, N_∑T_Plot)

			# Root Water Uptake daily 
				# Initial condition
            Date_Plot[1]                = Date_Start_Calibr
            ΔEvaporation_Plot[1]        = 0.0
            ΔFlux_Plot[1,1:N_iZ+1]     .= 0.0
            ΔPet_Plot[1]                = 0.0
				ΔPr_Plot[1]                 = 0.0
				ΔPrGross_Plot[1]=0.0
            ΔSink_Plot[1]               = 0.0

				∑T_Plot[1] = ∑T_Plot[1] * cst.Second_2_Hour
				for iT=2:N_∑T_Plot
					∑T_Plot[iT] = ∑T_Plot[iT] .* cst.Second_2_Hour

					Date_Plot[iT] = Date_Plot[iT-1] + Second(Int(ceil((∑T_Plot[iT]-∑T_Plot[iT-1]) * cst.Hour_2_Second)))

					for iZ=1:N_iZ+1
						ΔFlux_Plot[iT,iZ] = ∑Flux_Plot[iT,iZ] - ∑Flux_Plot[iT-1,iZ]
					end

               ΔSink_Plot[iT]        = ∑∑Sink[iT] - ∑∑Sink[iT-1]

               ΔPet_Plot[iT]         = ∑Pet_Plot[iT] - ∑Pet_Plot[iT-1]

					ΔPr_Plot[iT]          = max(∑Pr_Plot[iT] - ∑Pr_Plot[iT-1], 0.0)

					ΔPrGross_Plot[iT]      = max(∑PrGross_Plot[iT] - ∑PrGross_Plot[iT-1], 0.0)
					
					ΔEvaporation_Plot[iT] = max(∑Evaporation_Plot[iT] -  ∑Evaporation_Plot[iT-1], 0.0) 
				end  # for iT=1:N_iT
		

			# PONDING: Finding when Ponding occures (non zero). Note that not corrected for Δt
				iΔPond_NonZero = findall(!iszero, ΔHpond[1:N_iT])

			ΔPond_Plot = fill(0.0, N_∑T_Plot)
			ΔPond_Plot = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Plot, N_∑T_Plot, N_iT, ΔPond_Plot, ΔHpond)
			if !(isempty(iΔPond_NonZero))
				Flag_Plot_Pond = true
			else
				Flag_Plot_Pond = false
			end
		
		return ∑T_Plot, ∑T_Date_Plot, ∑WaterBalance_η_Plot, Date_Plot, Flag_Plot_Pond, N_∑T_Plot, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot
	end # function: CHANGE_OUTPUT_ΔT
	
end  # module: tincrease
# ............................................................