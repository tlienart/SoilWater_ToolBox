# =============================================================
#		MODULE: checkerror

# =============================================================
module checkerror
	import ..option, ..param, ..path
	import Dates: value, DateTime
	export CHECK_IFOPEN

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CHECK_ERROR
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function CHECK_ERROR(N_iRoot, N_iZ, N_iHorizon, Z, clim, veg, hydroHorizon)

		# DETERMENING IF PATH IS OPEN
			CHECK_IFOPEN(path.Table_Discretisation)
			CHECK_IFOPEN(path.Table_TimeSerie)
			CHECK_IFOPEN(path.Table_θ)
			CHECK_IFOPEN(path.Table_Q)
			CHECK_IFOPEN(path.Table_Ψ)
			CHECK_IFOPEN(path.Table_WaterBalance)
			CHECK_IFOPEN(path.Table_RootWaterUptake)
			CHECK_IFOPEN(path.Table_Se)

		# CHECKING IF THE OPTIONS ARE VALID	
			if option.hypix.θΨKmodel ≠ :Kosugi && option.hypix.θΨKmodel ≠ :vanGenuchten
				error("\n Hypix error: θΨKmodel option = $θΨKmodel not yet supported. θΨKmodel must = either [vanGenuchten] or [Kosugi]")
			end

			if option.hypix.BottomBoundary ≠ :Free && option.hypix.BottomBoundary ≠ :Pressure
				error("\n Hypix error: BottomBoundary option = $BottomBoundary not yet supported. BottomBoundary must = either [Free] or [Pressure]")
			end

			Date_Start = DateTime(param.hypix.Year_Start, param.hypix.Month_Start, param.hypix.Day_Start, param.hypix.Hour_Start, param.hypix.Minute_Start, param.hypix.Second_Start)

			Date_End = DateTime(param.hypix.Year_End, param.hypix.Month_End, param.hypix.Day_End, param.hypix.Hour_End, param.hypix.Minute_End, param.hypix.Second_End)

			if Date_End < Date_Start
				error("\n Hypix error: End Run Data = $(Date_End) before Start Run Data = $(Date_Start) !!!")
			end

		# CHECKING HYDRO PARAMETERS
			if option.hypix.θΨKmodel == :Kosugi
				for iHorizon in 1:N_iHorizon
					if hydroHorizon.θs[iHorizon] <  hydroHorizon.θsMacMat[iHorizon]
						error("\n Hypix error: at iHorizon = $iHorizon θs must be ≥ θsMacMat : $(path.Hydraulic)")
					end
				end # for iHorizon in 1:N_iHorizon
			end # option.hypix.θΨKmodel
		
		# CHECKING THE ROOT DENSITY PARAMETERS
			if option.hypix.RootWaterUptake
				CHECK_ROOTDISTRIBUTION(veg, N_iRoot, Z)
			end

		# CHECKING STARTING & ENDING DATES OF PLOTS 
			if option.Plot
				CHECK_DATES_PLOTS(clim, param)
			end

		# CHECKING 
			if maximum(param.hypix.plot.Cells_Plot) ≥ N_iZ 
				error("\n Hypix error:  param.hypix.plot.Cells_Plot = $(param.hypix.plot.Cells_Plot) must be ≤  N_iZ =  $N_iZ")
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
	function CHECK_DATES_PLOTS(clim, param)

	# 	if !(DateTime(clim.Year[1], clim.Month[1], clim.Day[1], clim.Hour[1], clim.Minute[1], clim.Second[1]) ≤ DateTime(param.hypix.plot.Year_Start, param.hypix.plot.Month_Start, param.hypix.plot.Day_Start, param.hypix.plot.Hour_Start, param.hypix.plot.Minute_Start, param.hypix.plot.Second_Start) < DateTime(param.hypix.plot.Year_End, param.hypix.plot.Month_End, param.hypix.plot.Day_End, param.hypix.plot.Hour_End, param.hypix.plot.Minute_End, param.hypix.plot.Second_End) ≤ DateTime(clim.Year[clim.N_Climate], clim.Month[clim.N_Climate], clim.Day[clim.N_Climate], clim.Hour[clim.N_Climate], clim.Minute[clim.N_Climate], clim.Second[clim.N_Climate]))
				
	# 		return error("\n \n Hypix ERROR: Dates of plot Year_Start=$(param.hypix.plot.Year_Start), Month_Start=$(param.hypix.plot.Month_Start), Day_Start=$(param.hypix.plot.Day_Start), Hour_Start$(param.hypix.plot.Hour_Start), Minute_Start= $(param.hypix.plot.Minute_Start), Second= $(param.hypix.plot.Second_Start) .OR. Year_End=$(param.hypix.plot.Year_End), Month_End=$(param.hypix.plot.Month_End), Day_End=$(param.hypix.plot.Day_End), Hour_End$(param.hypix.plot.Hour_End), Minute_End= $(param.hypix.plot.Minute_End), Second= $(param.hypix.plot.Second_End)")
	# 	end # if 
	return
	end  # function CHECK_DATES_PLOTS

end  # module checkerror