		# =============================================================
		#		module: plotKsModel
		# =============================================================
		module plotKsModel
		import ..cst
			using GLMakie

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : KsModel_3D
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KsModel_3D(hydro, NiZ, path)

				# GLMakie.activate!()
				# ; shading=false, colormap = :deep, axis = (show_axis = false,)

				
				Fig = Figure()
				Axis1= Axis(Fig[1,1])

				Z1 = Float64.(hydro.θs[1:NiZ] .- hydro.θr[1:NiZ])
				Z2 = Float64.(hydro.σ[1:NiZ])
				Z3 = Float64.(cst.MmS_2_MmH.*hydro.Ks[1:NiZ])

				surface!(Axis1, Z1, Z2, Z3; shading=true, colormap = :deep, axis = (show_axis = true,))

				Axis2 = Axis(Fig[2,1])
				Z1b = Float64.(hydro.θsMacMat[1:NiZ] .- hydro.θr[1:NiZ])
				surface!(Axis2, Z1b, Z2, Z3; shading=true, colormap = :deep, axis = (show_axis = true,))

				Axis3 = Axis(Fig[3,1])
				Z1c = Float64.(hydro.θsMacMat[1:NiZ] .- hydro.θr[1:NiZ])
				surface!(Axis3, Z1c, Z2, Z3; shading=true, colormap = :deep, axis = (show_axis = true,))

				Path = "D:/Main/MODELS/SoilWater_ToolBox/data/OUTPUT/SoilWater/Int/Plots/KsModel/KsDModel.png"
				save(Path, Fig)

			return nothing
			end  # function: KsModel_3D
			# ------------------------------------------------------------------
			
		end  # module: plotKsModel
		# ............................................................
		
