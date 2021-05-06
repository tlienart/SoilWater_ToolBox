# =============================================================
#		module: discret
# =============================================================
module discretization
	export DISCRETIZATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  DISCRETIZATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	struct DISCRET
		N_iZ
		Z_CellUp
		Znode
		ΔZ
		ΔZ_⬓
		ΔZ_Aver
		ΔZ_W
	end
	
	function DISCRETIZATION(N_iZ::Int64, Z)
      Z_CellUp = Vector{Float64}(undef, N_iZ)
      Znode    = Vector{Float64}(undef, N_iZ)
      ΔZ       = Vector{Float64}(undef, N_iZ)
      ΔZ_⬓     = Vector{Float64}(undef, N_iZ)
      ΔZ_Aver  = Vector{Float64}(undef, N_iZ+1)
      ΔZ_W     = Vector{Float64}(undef, N_iZ+1)
		
		# Layer 1
         ΔZ[1]       = Z[1]
         ΔZ_⬓[1]     = ΔZ[1] * 0.5
         Z_CellUp[1] = 0.0
         Znode[1]    = ΔZ_⬓[1]
         ΔZ_Aver[1]  = ΔZ_⬓[1]
         ΔZ_W[1]     = 1.0

		# All soils
		for iZ = 2:N_iZ
			ΔZ[iZ]       = Z[iZ] - Z[iZ-1]

			ΔZ_⬓[iZ]    = ΔZ[iZ] * 0.5

			Znode[iZ]    = Z[iZ] - ΔZ_⬓[iZ]

			ΔZ_Aver[iZ]  = (ΔZ[iZ] + ΔZ[iZ-1]) * 0.5

			ΔZ_W[iZ]     = ΔZ[iZ] /  (ΔZ[iZ] + ΔZ[iZ-1])

			Z_CellUp[iZ] = Z[iZ] - ΔZ[iZ]
		end # for

		ΔZ_Aver[N_iZ+1] = ΔZ[N_iZ] * 0.5

		ΔZ_W[N_iZ+1] = 0.0

		return discret = DISCRET(N_iZ, Z_CellUp, Znode, ΔZ, ΔZ_⬓, ΔZ_Aver, ΔZ_W)

	end # function DISCRETIZATION
	
end  # module: discret
# ............................................................