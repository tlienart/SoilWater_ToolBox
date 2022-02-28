# =============================================================
#		MODULE: checkerror

# =============================================================
module checkError
	import Dates: value, DateTime
	export CHECK_IFOPEN

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CHECK_ERROR
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function CHECK_ERROR(clim, hydroHorizon, N_Layer, N_iRoot, NiZ, paramHypix, veg, Z)

		# DETERMENING IF PATH IS OPEN
			CHECK_IFOPEN(pathHyPix.Table_Discretisation)
			CHECK_IFOPEN(pathHyPix.Table_TimeSerie)
			CHECK_IFOPEN(pathHyPix.Table_θ)
			CHECK_IFOPEN(pathHyPix.Table_Q)
			CHECK_IFOPEN(pathHyPix.Table_Ψ)
			CHECK_IFOPEN(pathHyPix.Table_WaterBalance)
			CHECK_IFOPEN(pathHyPix.Table_Se)

		# CHECKING IF THE OPTIONS ARE VALID	
			if optionHypix.Model⍰ ≠ "Kosugi" && optionoptionHypix ≠ "Vangenuchten"
				error("\n Hypix error: HydroModel⍰ optionHypix = $HydroModel⍰ not yet supported. HydroModel⍰ must = either [Vangenuchten] or [Kosugi]")
			end

			if optionHypix.BottomBoundary⍰ ≠ "Free" && optionHypix.BottomBoundary⍰ ≠ "Ψ"
				error("\n Hypix error: BottomBoundary⍰ optionHypix = $BottomBoundary⍰ not yet supported. BottomBoundary⍰ must = either [Free] or [Pressure]")
			end

			Date_Start = DateTime(paramHypix.Year_Start, paramHypix.Month_Start, paramHypix.Day_Start, paramHypix.Hour_Start, paramHypix.Minute_Start, paramHypix.Second_Start)

			Date_End = DateTime(paramHypix.Year_End, paramHypix.Month_End, paramHypix.Day_End, paramHypix.Hour_End, paramHypix.Minute_End, paramHypix.Second_End)

			if Date_End < Date_Start
				error("\n Hypix error: End Run Data = $(Date_End) before Start Run Data = $(Date_Start) !!!")
			end

		# CHECKING HYDRO PARAMETERS
			if optionHypix.Model⍰ == "Kosugi"
				for iHorizon in 1:N_Layer
					if hydroHorizon.θs[iHorizon] <  hydroHorizon.θsMacMat[iHorizon]
						error("\n Hypix error: at iHorizon = $iHorizon θs must be ≥ θsMacMat : $(pathHyPix.Hydraulic)")
					end
				end # for iHorizon in 1:N_Layer
			end # optionHypix.Model⍰
		
		# CHECKING THE ROOT DENSITY PARAMETERS
			if optionHypix.RootWaterUptake
				CHECK_ROOTDISTRIBUTION(veg, N_iRoot, Z)
			end



		return
	end  # function CHECK_ERROR
	

	
	# =====================================
	# 		CHECK IF FILE IS OPEN
	# =====================================
	function CHECK_IFOPEN(Path)
		try
			if isfile(Path) # If the file is there than delete it before saving figure
				rm(Path, force=true, recursive=true)
			end
			return
		catch
			error("\n Hypix ERROR: File open please close file before running = : ", Path, "\n \n")
		end
	end # function CHECK_IFOPEN



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CHECK_ROOTDISTRIBUTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function CHECK_ROOTDISTRIBUTION(veg, N_iRoot, Z)
		RootTop_Perc_Min = veg.Zroot_Top / Z[N_iRoot]
		if veg.ΔRdf_Top < RootTop_Perc_Min
			error("\n \n Hypix ERROR: ΔRdf_Top must be >   $RootTop_Perc_Min   for Zroot = $(Z[N_iRoot]) and Zroot_Top = $(veg.Zroot_Top)
			\n \n \n")
		end
		return
	end  # function CHECK_ROOTDISTRIBUTION



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CHECK_DATES_PLOTS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function CHECK_DATES_PLOTS(clim, paramHypix)

	# 	if !(DateTime(clim.Year[1], clim.Month[1], clim.Day[1], clim.Hour[1], clim.Minute[1], clim.Second[1]) ≤ DateTime(paramHypix.plot.Year_Start, paramHypix.plot.Month_Start, paramHypix.plot.Day_Start, paramHypix.plot.Hour_Start, paramHypix.plot.Minute_Start, paramHypix.plot.Second_Start) < DateTime(paramHypix.plot.Year_End, paramHypix.plot.Month_End, paramHypix.plot.Day_End, paramHypix.plot.Hour_End, paramHypix.plot.Minute_End, paramHypix.plot.Second_End) ≤ DateTime(clim.Year[clim.N_Climate], clim.Month[clim.N_Climate], clim.Day[clim.N_Climate], clim.Hour[clim.N_Climate], clim.Minute[clim.N_Climate], clim.Second[clim.N_Climate]))
				
	# 		return error("\n \n Hypix ERROR: Dates of plot Year_Start=$(paramHypix.plot.Year_Start), Month_Start=$(paramHypix.plot.Month_Start), Day_Start=$(paramHypix.plot.Day_Start), Hour_Start$(paramHypix.plot.Hour_Start), Minute_Start= $(paramHypix.plot.Minute_Start), Second= $(paramHypix.plot.Second_Start) .OR. Year_End=$(paramHypix.plot.Year_End), Month_End=$(paramHypix.plot.Month_End), Day_End=$(paramHypix.plot.Day_End), Hour_End$(paramHypix.plot.Hour_End), Minute_End= $(paramHypix.plot.Minute_End), Second= $(paramHypix.plot.Second_End)")
	# 	end # if 
	return
	end  # function CHECK_DATES_PLOTS

end  # module checkerror