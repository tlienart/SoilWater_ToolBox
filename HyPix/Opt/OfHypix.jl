# =============================================================
#		MODULE: ofHypix
# =============================================================
module ofHypix

	# =============================================================
	#		MODULE: θ
	# =============================================================
	module θof
		import ...tool, ...param, ...interpolate
		import Statistics
		export WOF_θ, RMSE_θ

		const PowerGood = 2

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : WOF_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function WOF_θ(∑T, calibr, N_iT::Int, N_iZ::Int, ΔHpond, θ, θSim; θobs_Uncert=param.hypix.calibr.θobs_Uncert)

			θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, calibr.∑T[1:calibr.N_iT], calibr.N_iT, N_iT, N_iZ, θSim, θ)
			
			Wof = 0.0
			iCount = 0

			for iZ=1:calibr.Ndepth

				Wdepth = 2.0 * (Float64(calibr.Ndepth) + 1.0 - Float64(iZ) ) / (Float64(calibr.Ndepth) * (Float64(calibr.Ndepth) + 1.0))

				for iT=1:calibr.N_iT
				
					if θSim[iT,iZ] > 0.0 # avoiding no data

					# # 	# Error
							θerror = abs(calibr.θobs[iT,iZ] - θSim[iT, calibr.iZobs[iZ]])

					# 	# Taking into account uncertainty bar
					# 		θobs_Min = calibr.θobs[iT,iZ] - θobs_Uncert

					# 		θobs_Max = calibr.θobs[iT,iZ] + θobs_Uncert

					# 	if θobs_Min ≤ θSim[iT, calibr.iZobs[iZ]] ≤ θobs_Max

					# 		# Wof += Wdepth * θobs_Uncert  / 10.0

					# 		Wof += Wdepth * θobs_Uncert * (θerror / θobs_Uncert) ^ PowerGood
							
					# 	else
					# 		Wof += Wdepth * θerror
	
							Wof += Wdepth * θerror ^ PowerGood
						# end # if θobs_Min ≤ θSim[iT, calibr.iZobs[iZ]] ≤ θobs_Max

						iCount += 1

					end # if: calibr.θobs[iT,iZ] > 0.0
				end # for iT
			end # for iZ

			 # Penalty if we have too much ponding
			 Wof_Pond = max(ΔHpond[N_iZ] - param.hypix.ΔHpondMax, 0.0) / 1000.0

		return Wof = √(Wof / Float64(iCount)) + Wof_Pond
		end # function WOF_θ



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : RMSE
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function RMSE_θ(∑T, calibr, N_iT::Int, N_iZ::Int, θ, θSim)

				θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, calibr.∑T[1:calibr.N_iT], calibr.N_iT, N_iT, N_iZ, θSim, θ)
				
				Rmse = 0.0
				iCount = 0
				for iZ=1:calibr.Ndepth

					Wdepth = 2.0 * (Float64(calibr.Ndepth) + 1.0 - Float64(iZ) ) / (Float64(calibr.Ndepth) * (Float64(calibr.Ndepth) + 1.0))

					for iT=1:calibr.N_iT 	
						if θSim[iT,iZ] > 0.0 # avoiding no data
							Rmse +=  Wdepth * (calibr.θobs[iT,iZ] - θSim[iT, calibr.iZobs[iZ]]) ^ 2
							iCount += 1
						end # if: calibr.θobs[iT,iZ] > 0.0
					end # for it

				end # for iZ

			return Rmse = √(Rmse / (Float64(iCount)))		
			end # function WOF_θ


	end  # module θof

	# ............................................................

end  # module ofHypix
# ............................................................