# =============================================================
#		MODULE: tincrease
# =============================================================
module Δtchange
	import ..interpolate, ..reading, ..tool, ..cst
	import Dates: value, DateTime, Day, Second, Hour, now

	function CHANGE_OUTPUT_ΔT(∑Pet, ∑Pr, ∑T, ∑WaterBalance_η, ∑ΔSink, obsTheta, clim, Nit::Int64, NiZ::Int64, param, Q, ΔEvaporation, ΔHpond, ΔT, θ, Ψ, ∑T_Climate)

		# PREPROCESSING ∑Evaporation, ∑ΔQ
         ∑Evaporation = fill(0.0::Float64, Nit)
         ∑ΔQ          = fill(0.0::Float64, Nit, NiZ+1)

			∑Evaporation[1]    = 0.0
         ∑ΔQ[1, 1:NiZ+1]  .= 0.0
			for iT=2:Nit
				∑Evaporation[iT] = ∑Evaporation[iT-1] + ΔEvaporation[iT]
				for iZ=1:NiZ+1
					∑ΔQ[iT,iZ] = ∑ΔQ[iT-1, iZ] + ΔT[iT] * Q[iT,iZ]
				end
			end

		# ∑Pr_Gross
			∑Pr_Gross = fill(0.0::Float64, clim.N_Climate)
			∑Pr_Gross[1] = 0.0
			for iT= 2:clim.N_Climate
				∑Pr_Gross[iT] = ∑Pr_Gross[iT-1] + clim.Pr[iT]
			end

		# PREPARING DATA FOR PLOTS
			ΔT_Sim = value(obsTheta.Date[end] - obsTheta.Date[1]) / 1000

			∑T_Reduced = collect(range(0.0, step=param.hyPix.ΔT_Output, stop=ΔT_Sim)) 
			
			# Take account that we are starting at Date_Start_Calibr
			ΔT_Start_Calibr = value(obsTheta.Date[1] - clim.Date[1]) / 1000

			@. ∑T_Reduced = ∑T_Reduced + ΔT_Start_Calibr

			Nit_Reduced = length(∑T_Reduced)

		# PREPARING DATES WITH INTERVAL:
			∑T_Date_Reduced = range(obsTheta.Date[1], step=Second(param.hyPix.ΔT_Output), obsTheta.Date[end])

		# INTERPOLATING DATA
			θ_Reduced = fill(0.0::Float64, Nit_Reduced, NiZ)	
				θ_Reduced = interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, NiZ, θ_Reduced, θ)

			θobs_Reduced = fill(0.0::Float64, Nit_Reduced, obsTheta.Ndepth)
				θobs_Reduced = interpolate.INTERPOLATE_2D_LOOP(obsTheta.∑T[1:obsTheta.Nit], ∑T_Reduced, obsTheta.Nit, obsTheta.Ndepth, θobs_Reduced, obsTheta.θobs[1:obsTheta.Nit,1:obsTheta.Ndepth])

			Ψ_Reduced = fill(0.0::Float64, Nit_Reduced, NiZ)
				Ψ_Reduced =  interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, NiZ, Ψ_Reduced, Ψ) 

			∑ΔQ_Reduced = fill(0.0::Float64, Nit_Reduced, NiZ+1)
				∑ΔQ_Reduced = interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, NiZ+1, ∑ΔQ_Reduced, ∑ΔQ)

		# .<>.<>.<>
			ΔT_Reduced = fill(0.0::Float64, Nit_Reduced)
				ΔT_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ΔT_Reduced, ΔT)

			∑Evaporation_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑Evaporation_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑Evaporation_Reduced, ∑Evaporation)

			∑Pet_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑Pet_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑Pet_Reduced, ∑Pet)

			∑Pr_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑Pr_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑Pr_Reduced, ∑Pr)

			∑WaterBalanceη_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑WaterBalanceη_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑WaterBalanceη_Reduced, ∑WaterBalance_η)

			∑∑Sink = fill(0.0::Float64, Nit_Reduced)
				∑∑Sink = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑∑Sink, ∑ΔSink)

			∑PrGross_Reduced = fill(0.0::Float64, Nit_Reduced)	
				∑PrGross_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T_Climate[1:clim.N_Climate], ∑T_Reduced, Nit_Reduced, clim.N_Climate, ∑PrGross_Reduced, ∑Pr_Gross)

			ΔPond_Reduced = fill(0.0, Nit_Reduced)
				ΔPond_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ΔPond_Reduced, ΔHpond)

			# From ∑ to Δ
            ΔEvaporation_Reduced = fill(0.0::Float64, Nit_Reduced)
            ΔQ_Reduced        = fill(0.0::Float64, Nit_Reduced, NiZ+1)
            ΔPet_Reduced         = fill(0.0::Float64, Nit_Reduced)
				ΔPr_Reduced          = fill(0.0::Float64, Nit_Reduced)
				ΔPrGross_Reduced     = fill(0.0::Float64, Nit_Reduced)
            ΔSink_Reduced        = fill(0.0::Float64, Nit_Reduced)
				Date_Reduced         = fill(now()::DateTime, Nit_Reduced)

			# Root Water Uptake daily 
				# Initial condition
            Date_Reduced[1]                 = obsTheta.Date[1]
            ΔEvaporation_Reduced[1]         = 0.0
            ΔQ_Reduced[1,1:NiZ+1]      .= 0.0
            ΔPet_Reduced[1]                 = 0.0
            ΔPr_Reduced[1]                  = 0.0
            ΔPrGross_Reduced[1]             = 0.0
            ΔSink_Reduced[1]                = 0.0

			# ∑T_Reduced[1] = ∑T_Reduced[1]
			for iT=2:Nit_Reduced
				# ∑T_Reduced[iT] = ∑T_Reduced[iT]

				# Date_Reduced[iT] = Date_Reduced[iT-1] + Second(Int(ceil((∑T_Reduced[iT]-∑T_Reduced[iT-1]) * cst.Hour_2_Second)))

				for iZ=1:NiZ+1
					ΔQ_Reduced[iT,iZ] = ∑ΔQ_Reduced[iT,iZ] - ∑ΔQ_Reduced[iT-1,iZ]
				end

            ΔSink_Reduced[iT]        = ∑∑Sink[iT] - ∑∑Sink[iT-1]

            ΔPet_Reduced[iT]         = ∑Pet_Reduced[iT] - ∑Pet_Reduced[iT-1]

            ΔPr_Reduced[iT]          = max(∑Pr_Reduced[iT] - ∑Pr_Reduced[iT-1], 0.0) # Numerical stability

            ΔPrGross_Reduced[iT]     = max(∑PrGross_Reduced[iT] - ∑PrGross_Reduced[iT-1], 0.0)
				
            ΔEvaporation_Reduced[iT] = max(∑Evaporation_Reduced[iT] -  ∑Evaporation_Reduced[iT-1], 0.0)
			end  # for iT=1:Nit
		
	return ∑T_Date_Reduced, ∑T_Reduced, ∑WaterBalanceη_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔQ_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced
	end # function: CHANGE_OUTPUT_ΔT
	
end  # module: tincrease
# ............................................................