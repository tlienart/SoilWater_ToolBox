# =============================================================
#		module: discret
# =============================================================
module discretization
	export DISCRETIZATION, DISCRETIZATION_AUTO

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  DISCRETIZATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		struct DISCRET
			NiZ
			Z_CellUp
			Znode
			ΔZ
			ΔZ_⬓
			ΔZ_Aver
			ΔZ_W
		end
		
		function DISCRETIZATION(NiZ::Int64, Z)
			Z_CellUp = fill(0.0::Float64, NiZ)
			Znode    = fill(0.0::Float64, NiZ)
			ΔZ       = fill(0.0::Float64, NiZ+1)
			ΔZ_⬓     = fill(0.0::Float64, NiZ)
			ΔZ_Aver  = fill(0.0::Float64, NiZ+1)
			ΔZ_W     = fill(0.0::Float64, NiZ+1)
			
			# Cell 1
				ΔZ[1]       = Z[1]
				ΔZ_⬓[1]     = ΔZ[1] * 0.5
				Z_CellUp[1] = 0.0
				Znode[1]    = ΔZ_⬓[1]
				ΔZ_Aver[1]  = ΔZ_⬓[1]
				ΔZ_W[1]     = 1.0

			# All Cells
			for iZ = 2:NiZ
				ΔZ[iZ]       = Z[iZ] - Z[iZ-1]

				ΔZ_⬓[iZ]    = ΔZ[iZ] * 0.5

				Znode[iZ]    = Z[iZ] - ΔZ_⬓[iZ]

				ΔZ_Aver[iZ]  = (ΔZ[iZ] + ΔZ[iZ-1]) * 0.5

				ΔZ_W[iZ]     = ΔZ[iZ] /  (ΔZ[iZ] + ΔZ[iZ-1])

				Z_CellUp[iZ] = Z[iZ] - ΔZ[iZ]
			end # for

			ΔZ_Aver[NiZ+1] = ΔZ[NiZ]

			ΔZ_W[NiZ+1] = 0.5

			ΔZ[NiZ+1]  = ΔZ[NiZ]

		return discret = DISCRET(NiZ, Z_CellUp, Znode, ΔZ, ΔZ_⬓, ΔZ_Aver, ΔZ_W)
		end # function DISCRETIZATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DISCRETIZATION_AUTO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	""" 
					DISCRETIZATION_AUTO(N_Layer, Zlayer Zroot)

	Automatically performs the discretisatio of the HyPix model wheh you enter the depth of the layers
	"""
		function DISCRETIZATION_AUTO(param; Flag_θΨini, N_Layer, Zlayer, θini, Ψini)

			# Determine if we selected to input θini or Ψini
				if Flag_θΨini ==:Ψini
					θΨini = Ψini
				elseif Flag_θΨini ==:θini
					θΨini = θini
				end

			ΔZlayer = fill(0.0::Float64, N_Layer)

			# Computing ΔZlayer
				ΔZlayer[1]= Zlayer[1]
				for iZ = 2:N_Layer
					ΔZlayer[iZ] = Zlayer[iZ] - Zlayer[iZ-1]
				end # for

			# Computing the number of discretization
            ΔZcell     = []
            Layer      = []
            θΨini_Cell = []

				for iLayer = 1:N_Layer
					if  Zlayer[iLayer] < param.hyPix.ZfineCoarse
						ΔZ_Max = param.hyPix.ΔZfine
					else
						ΔZ_Max = param.hyPix.ΔZcoarse
					end

					Nsplit = ceil(ΔZlayer[iLayer] / ΔZ_Max) # Number of splitting from Layer->Cell
					ΔZcell₀ = ΔZlayer[iLayer] / Float64(Nsplit)

					for iDiscret=1:Nsplit
						append!(ΔZcell, ΔZcell₀)
                  append!(Layer, iLayer)
						append!(θΨini_Cell, θΨini[iLayer])
					end
				end # ilayer
				N = length(ΔZcell)

			# Computing the ∑ΔZcell
				Z = fill(0.0::Float64, N)

				Z[1] = ΔZcell[1]
				for iZ=2:N
					Z[iZ] = Z[iZ-1] + ΔZcell[iZ]
				end
			
		return Layer, Z, θΨini_Cell
		end  # function: DISCRETIZATION_AUTO


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DISCRETISATION_AUTO_θini
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	""" 
					DISCRETIZATION_AUTO(N_Layer, Zlayer)

	Automatically performs the discretisatio of the HyPix model wheh you enter the depth of the layers
	"""
		function DISCRETISATION_AUTO_θini(; N_Layer, Zlayer, Zroot, θᵢₙᵢ)
			HydrostaticEquilibrium=false

			ΔZlayer = fill(0.0::Float64, N_Layer)

			# Computing ΔZlayer
				ΔZlayer[1]= Zlayer[1]
				for iZ = 2:N_Layer
					ΔZlayer[iZ] = Zlayer[iZ] - Zlayer[iZ-1]
				end # for

			# Computing the number of discretization
            ΔZcell    = []
            Layer     = []
            θᵢₙᵢ_Cell = []
				for iLayer =1:N_Layer

					if  Zlayer[iLayer] < Zroot
						ΔZ_Max = param.hyPix.ΔZfine
					else
						ΔZ_Max = param.hyPix.ΔZcoarse
					end

					Nsplit = ceil(ΔZlayer[iLayer] / ΔZ_Max) # Number of splitting from Layer->Cell
					ΔZcell₀ = ΔZlayer[iLayer] / Float64(Nsplit)

					for iDiscret=1:Nsplit
						append!(ΔZcell, ΔZcell₀)
                  append!(Layer, iLayer)
						append!(θᵢₙᵢ_Cell, θᵢₙᵢ[iLayer])
					end
				end # ilayer
				N = length(ΔZcell)

			# Computing the ∑ΔZcell
				Z = fill(0.0::Float64, N)

				Z[1] = ΔZcell[1]
				for iZ=2:N
					Z[iZ] = Z[iZ-1] + ΔZcell[iZ]
				end

				if HydrostaticEquilibrium
					for iZ = 1:N
						θᵢₙᵢ_Cell = Z[N] - Z[iZ]
					end
				end
	
		return Layer, Z, θᵢₙᵢ_Cell
		end  # function: DISCRETIZATION_AUTO
	
end  # module: discret
# ............................................................