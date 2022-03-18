# =============================================================
#		MODULE: waterbalance
# =============================================================
module waterBalance

	export WATERBALANCE

	function WATERBALANCE(∑T, obsTheta, discret, hydro, N_iRoot::Int64, Nit::Int64, NiZ::Int64, Q, ΔSink, ΔT, θ, Ψ)
      ∑WaterBalance_η = fill(0.0::Float64, Nit)
      ∑ΔSink          = fill(0.0::Float64, Nit)
      ∑ΔQtop          = 0.0 ::Float64
      ∑ΔQbottom       = 0.0 ::Float64
      ∑∑WaterBalance  = 0.0 ::Float64
		ΔStorage        = 0.0 ::Float64
		ΔStorageSo      = 0.0 ::Float64
      i∑T_CalibrStart = 1::Int64

		# Starting to compute the waterbalance after the warmup period
			while ∑T[i∑T_CalibrStart] < obsTheta.∑T[1]
				i∑T_CalibrStart += 1
			end
		
		for iT=i∑T_CalibrStart:Nit
			# Computing ΔStorage
				ΔStorage = 0.0::Float64

				@fastmath @inbounds for iZ = 1:NiZ
					ΔStorage += discret.ΔZ[iZ] * ((θ[iT,iZ] - θ[i∑T_CalibrStart-1,iZ]))
	
					ΔStorageSo += discret.ΔZ[iZ] * hydro.So[iZ] * (Ψ[iT,iZ] - Ψ[iT-1,iZ]) * (θ[iT,iZ] / hydro.θs[iZ])
				end # for iT=1:NiZ

			# Sink term
				∑ΔSink[iT] = ∑ΔSink[iT-1]
				@fastmath @inbounds @simd for iZ = 1:N_iRoot
					∑ΔSink[iT] += ΔSink[iT,iZ]
				end # for: iZ = 1:NiZ

			# Cumulative water entering top cell
				∑ΔQtop += ΔT[iT] * Q[iT,1]

			# Cumulative water leaving bottom cell
				∑ΔQbottom += ΔT[iT] * Q[iT,NiZ+1]

			∑∑WaterBalance = ΔStorage - (∑ΔQtop - ∑ΔQbottom) + ∑ΔSink[iT] - ΔStorageSo

			# Normalised water balance, performed min algorithme to avoid issues when ∑ΔQtop is small at the beginning of time
				∑WaterBalance_η[iT] = min(abs(∑∑WaterBalance / ∑ΔQtop), abs(∑∑WaterBalance))
		end  # for iT=1:Nit

	return ∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage
	end  # function: WATERBALANCE
	#-------------------------------------------------------------------------------

end  # module waterbalance
# ............................................................