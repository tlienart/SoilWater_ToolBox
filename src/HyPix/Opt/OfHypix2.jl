# =============================================================
#		MODULE: ofHypix
# =============================================================
module ofHypix

	# =============================================================
	#		MODULE: θ
	# =============================================================
	module θof
		import ...tool, ...interpolate
		import Statistics
		export WOF_θ, RMSE_θ

		const OFmodel ="CCC"

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : WOF_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function WOF_θ(∑T, Nit::Int, NiZ::Int, obsTheta, param, ΔHpond, θ, θSim; θobs_Uncert=param.hyPix.obsTheta.θobs_Uncert)

			θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, obsTheta.∑T[1:obsTheta.Nit], Nit, NiZ, θSim, θ)
			
			Wof = 0.0
			iCount = 0

			for iZ=1:obsTheta.Ndepth

				Wdepth = 2.0 * (Float64(obsTheta.Ndepth) + 1.0 - Float64(iZ) ) / (Float64(obsTheta.Ndepth) * (Float64(obsTheta.Ndepth) + 1.0))

				for iT=1:obsTheta.Nit
				
					if θSim[iT,iZ] > 0.0  && obsTheta.θobs[iT,iZ] > 0.0# avoiding no data

					# # 	# Error
							θerror = abs(obsTheta.θobs[iT,iZ] - θSim[iT, obsTheta.ithetaObs[iZ]])

					# 	# Taking into account uncertainty bar
					# 		θobs_Min = obsTheta.θobs[iT,iZ] - θobs_Uncert

					# 		θobs_Max = obsTheta.θobs[iT,iZ] + θobs_Uncert

					# 	if θobs_Min ≤ θSim[iT, obsTheta.ithetaObs[iZ]] ≤ θobs_Max

					# 		# Wof += Wdepth * θobs_Uncert  / 10.0

					# 		Wof += Wdepth * θobs_Uncert * (θerror / θobs_Uncert) ^ 2.0
							
					# 	else
					# 		Wof += Wdepth * θerror
	
							Wof += Wdepth * θerror ^ 2.0
						# end # if θobs_Min ≤ θSim[iT, obsTheta.ithetaObs[iZ]] ≤ θobs_Max

						iCount += 1

					end # if: obsTheta.θobs[iT,iZ] > 0.0
				end # for iT
			end # for iZ

			 # Penalty if we have too much ponding
			 Wof_Pond = max(ΔHpond[NiZ] - param.hyPix.obsTheta.ΔHpondMax, 0.0) / 1000.0

		return Wof = √(Wof / Float64(iCount)) + Wof_Pond
		end # function WOF_θ



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : RMSE_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function RMSE_θ(∑T, obsTheta, Nit::Int, NiZ::Int, θ, θSim)


				θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, obsTheta.∑T[1:obsTheta.Nit], Nit, NiZ, θSim, θ)
				
				Rmse = 0.0
				iCount = 0
				for iZ=1:obsTheta.Ndepth

					Wdepth = 2.0 * (Float64(obsTheta.Ndepth) + 1.0 - Float64(iZ) ) / (Float64(obsTheta.Ndepth) * (Float64(obsTheta.Ndepth) + 1.0))

					for iT=1:obsTheta.Nit 	
						if θSim[iT,iZ] > 0.0 && obsTheta.θobs[iT,iZ] > 0.0 # avoiding no data
							Rmse +=  Wdepth * (obsTheta.θobs[iT,iZ] - θSim[iT, obsTheta.ithetaObs[iZ]]) ^ 2
							iCount += 1
						end # if: obsTheta.θobs[iT,iZ] > 0.0
					end # for it

				end # for iZ

			return Rmse = √(Rmse / (Float64(iCount)))		
			end # function WOF_θ

	end  # module θof

	# ............................................................

end  # module ofHypix
# ............................................................