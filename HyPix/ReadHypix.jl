# =============================================================
#		module: readingHypix
# =============================================================
module readHypix
	import ..path, ..option, ..tool, ..param, ..horizonLayer
	import Dates: value, DateTime, hour, minute, month
	import DelimitedFiles
	export CLIMATE, DISCRETIZATION, HYPIX_PARAM, LOOKUPTABLE_LAI, LOOKUPTABLE_CROPCOEFICIENT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DISCRETIZATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DISCRETIZATION()
			# Read data
				Data = DelimitedFiles.readdlm(path.Discretization, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

			Z, N_iZ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Z")
			θ_Ini, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "θini")
			Horizon, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Horizon")

			N_iHorizon = maximum(Horizon)
			return Horizon, N_iHorizon, N_iZ, Z, θ_Ini
		end # function DISCRETIZATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYPIX_PARAM(Horizon, hydro, hydroHorizon, iSim::Int64, N_iZ::Int64, veg)

			# Read data
				Data = DelimitedFiles.readdlm(path.Hypix_Param, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

			# Readingt the type of data
				Type, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "TYPE")

			# Reading the names of the parameters
				Name, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "NAME")

				Name_Unique = unique(Name)

				N_NameUnique = length(Name_Unique)
		
			# Reading the values of the parameters for the simulation of interest
				Param, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SIM_$(iSim)")
			
			# Minimum value of the param
				Param_Min, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MIN")

			# Maximum value of the param
				Param_Max, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MAX")

			# Determening which param to optimize
				Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT_$(iSim)")

			# Maximum value of the param
				Opt_LogTransform, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "LogTransform")
				
			# Determine if we need to optimize
				if sum(Opt) ≥ 1
					Flag_Opt = true
				else
					Flag_Opt = false
				end

			"""Determening if multistep optimisation is performed (not the first step)
			This is such that the optimal values of the previous optimisation step is kept in memory
			We need to determine what next param to optimize"""
				if Flag_Opt && (iSim ≥ param.hypix.iSim_Start + 1)
					Flag_MultiStepOpt = true
				else
					Flag_MultiStepOpt = false 
				end
				
			# ====================================================

			# Does not matter if repeated in multistep optimisation
         ParamOpt              = []
         ParamOpt_HorizonEq    = []
         ParamOpt_Max          = []
         ParamOpt_Min          = []
         ParamOpt_Type         = []
         ParamOpt_LogTransform = []
         # θs_Min                = []
         # θs_Max                = []

			for i in eachindex(Name_Unique)
				# Finding the position of each Param name in .csv
					indexName = findall(isequal(Name_Unique[i]), Name)

				# Values of param for every Name to put in hydroHorizon
					Param_Vect = Float64.(Param[indexName])

				if Type[i] == "hydro" && !(Flag_MultiStepOpt)
					# Putting soil param in hydroHorizon

					# θsMacMat value depends on θs
					if  Symbol(Name_Unique[i]) == :θsMacMat_ƞ 
						for iZ =1:length(Param_Vect)
							hydroHorizon.θsMacMat[iZ] =  hydroHorizon.θs[iZ] * Param_Vect[iZ]
						end
					end
					
					setfield!(hydroHorizon, Symbol(Name_Unique[i]), Param_Vect)

					# Minimum and maximum value of the hydraulic parameters such as θs_Min and θs_Max
						setfield!(hydroHorizon, Symbol(Name_Unique[i] * "_Min"), Float64.(Param_Min[indexName]))
						setfield!(hydroHorizon, Symbol(Name_Unique[i] * "_Max"), Float64.(Param_Max[indexName]))

				elseif Type[i] == "veg" && !(Flag_MultiStepOpt)
					# Putting veg param in veg
					setfield!(veg, Symbol(Name_Unique[i]), Float64(Param[indexName][1]))
				end
				
				# Param to optimize. The complication is that there are different layers of hydraulic parameters which can be optimized.  
				if sum(Opt[indexName]) > 0

					# Type of parameters
						append!(ParamOpt_Type, [Type[i]]) 

					# Appending name of param to optimize by removing dublicates
						append!(ParamOpt, [Name_Unique[i]])

					# Appending the iHorizon of the param which will have = values. The horizon to be optimized must follow 
						iHorizonOpt_Start = findfirst(x->x==1, Opt[indexName])
						iHorizonOpt_End = findlast(x->x==1, Opt[indexName])

						append!(ParamOpt_HorizonEq, [[iHorizonOpt_Start ; iHorizonOpt_End]])

					# Minimum and Maximum value of the parameter to be optimized. If we have layers than we use the value of the top layer
					iNameOpt = findfirst(x->x==Name_Unique[i], Name) + iHorizonOpt_Start - 1

					iStart = iNameOpt
					iEnd = iNameOpt + iHorizonOpt_End - iHorizonOpt_Start

					# We take the minimum to be the minimum of iHorizonOpt_Start and iHorizonOpt_End and the same for maximum
					append!(ParamOpt_Min, minimum(Param_Min[iStart:iEnd]))
					append!(ParamOpt_Max, maximum(Param_Max[iStart:iEnd]))

					# Appending name of param to perform logTransform if optimized by removing dublicates
					if sum(Opt_LogTransform[indexName]) > 0
						append!(ParamOpt_LogTransform, [true])

						ParamOpt_Min[end] = log1p(ParamOpt_Min[end])
						ParamOpt_Max[end] = log1p(ParamOpt_Max[end])
					else
						append!(ParamOpt_LogTransform, [false])
					end

					if Param_Min[iNameOpt] > Param_Max[iNameOpt]
						error("HYPIX ERROR: $(Param_Min[iNameOpt]) < $(Name_Unique[i]) < $(Param_Max[iNameOpt]) !")
					end

				end
			end # for loop

			# # Getting the minimum and maximum value of θs for every iHorizonOpt_Start::iHorizonOpt_End
			# i=1
			# for iName in Name
			# 	if  Symbol(iName) == :θs 
			# 		append!(θs_Min, Param_Min[i])
			# 		append!(θs_Max, Param_Max[i])
			# 	end
			# 	i+=1
			# end

			if !(Flag_MultiStepOpt)
				# Hydraulic parameters per horizon to layers
				hydro = horizonLayer.HYDROHORIZON_2_HYDRO(N_iZ, Horizon, hydroHorizon)
			end

			NparamOpt = length(ParamOpt)

			# CHECKING FOR UNCONSISTENCY WITH OPTIONS	
			if Flag_Opt && option.hypix.σ_2_Ψm ≠ :No && "Ψm" ∈ ParamOpt
				iψm = findfirst(isequal("Ψm"), ParamOpt)[1]

				if option.hypix.σ_2_Ψm==:UniqueRelationship && "Ψm" ∈ ParamOpt
					error( "**** HyPix Error: combination of options which are not possible (option.hypix.σ_2_Ψm==:UniqueRelationship) && (Optimise=Ψm)!")

				elseif option.hypix.σ_2_Ψm==:Constrained && !("Ψm" ∈ ParamOpt)
					error("*** HyPix Error: combination of options which are not possible (option.hypix.σ_2_Ψm==:Constrained) && (not Optimising=Ψm)!")

				elseif option.hypix.σ_2_Ψm==:Constrained && ParamOpt_LogTransform[iψm]==1
					error("*** option.hypix.σ_2_Ψm==:Constrained CANNOT log transforme Ψm") 
				end
			end # Flag_Opt

			# Putting all the parameters in  NamedTuple
			optim = (ParamOpt_Min=ParamOpt_Min, ParamOpt_Max=ParamOpt_Max, ParamOpt_HorizonEq=ParamOpt_HorizonEq, ParamOpt_Type=ParamOpt_Type, ParamOpt=ParamOpt, NparamOpt=NparamOpt, Flag_Opt=Flag_Opt, ParamOpt_LogTransform=ParamOpt_LogTransform)
			# θs_Min=θs_Min, θs_Max=θs_Max

			if Flag_Opt == true
				println("		=== === Optimizing the following parameters === ===")
				println("			NparamOpt=" , NparamOpt)
				println("			ParamOpt= " , optim.ParamOpt_Type .* optim.ParamOpt)
				println("			Min_Value= " , optim.ParamOpt_Min)
				println("			Max_Value= " , optim.ParamOpt_Max)
				println("			Hydro_HorizonEq= " , optim.ParamOpt_HorizonEq)
				println("			LogTransform = " , optim.ParamOpt_LogTransform)
				println("		=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n \n")
			end

			return hydro, hydroHorizon, optim, veg
	end  # function: HYPIX_PARAM



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CLIMATE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		struct CLIMATEDATA
			Date      :: Vector{DateTime}
			Pr        :: Vector{Float64}
			Pet       :: Vector{Float64}
			Temp      :: Vector{Float64}
			N_Climate :: Int64
			Pr_Through :: Vector{Float64}
		end

		function CLIMATE()
			if option.hypix.ClimateDataTimestep == "Daily"
				Pr_Name          = "Rain(mm)"
				Pet_Name         = "PET(mm)"
				Temperature_Name = "Tmax(C)" # Maximum temperature which is not correct

			elseif option.hypix.ClimateDataTimestep == "Hourly"
				Pr_Name          = "Pr_mm"
				Pet_Name         = "Pet_mm"
				Temperature_Name = "Temp_c"
			end #  option.hypix.ClimateDataTimestep

			# Read data
				Data = DelimitedFiles.readdlm(path.Climate, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

			Year, N_Climate = tool.readWrite.READ_HEADER_FAST(Data, Header,"Year")
			Month, ~        = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month")
			Day, ~          = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day")
			Hour, ~         = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour")
			Minute, ~       = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute")
			Second, ~       = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second")
			Pr, ~           = tool.readWrite.READ_HEADER_FAST(Data, Header, Pr_Name)
			Pet, ~          = tool.readWrite.READ_HEADER_FAST(Data, Header, Pet_Name)
			Temp, ~         = tool.readWrite.READ_HEADER_FAST(Data, Header, Temperature_Name)
			
			# REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
				Date_Start = DateTime(param.hypix.Year_Start, param.hypix.Month_Start, param.hypix.Day_Start, param.hypix.Hour_Start, param.hypix.Minute_Start, param.hypix.Second_Start)
				
				Date_End = DateTime(param.hypix.Year_End, param.hypix.Month_End, param.hypix.Day_End, param.hypix.Hour_End, param.hypix.Minute_End, param.hypix.Second_End)

			# CHECKING THAT
				# Date_End is feasible
					Date_End_Maximum = DateTime(Year[N_Climate], Month[N_Climate], Day[N_Climate], Hour[N_Climate], Minute[N_Climate], Second[N_Climate]) 


					if Date_End_Maximum > Date_End
						Date_End = min(Date_End_Maximum, Date_End)
						println("		~ HyPix WARNING: Date_End not feasible so modified to match the data ")
					end #warning Date_End

				# Date_Start is feasible
					Date_Start_Minimum = DateTime(Year[2], Month[2], Day[2], Hour[2], Minute[2], Second[2]) 

					if Date_Start_Minimum > Date_Start
						Date_Start = max(Date_Start_Minimum , Date_Start)
						println("		~HyPix WARNING: Date_Start = $Date_Start not feasible so modified to match the data ")
					end #warning Date_Start

				# Initializing the Array
				True = falses(N_Climate)
				Date = Array{DateTime}(undef, N_Climate) 
				for iT=1:N_Climate
					Date[iT] = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])

					if (Date_Start ≤ Date[iT] ≤ Date_End)
						True[iT] = true
					end  # if: 
				end # iT=1:N_Climate

				# Need to include one date iT-1 at the beginning to compute ΔT
				iTrue_First = findfirst(True[1:N_Climate])
				True[iTrue_First-1] = true

				# New reduced number of simulations
					Date = Date[True[1:N_Climate]]
					Pr   = Pr[True[1:N_Climate]]
					Pet  = Pet[True[1:N_Climate]]
					Temp = Temp[True[1:N_Climate]]
					
					N_Climate = count(True[1:N_Climate]) # New number of data
			
			# To be used after interception model
			Pr_Through = zeros(Float64, N_Climate)


			# PUTTING CLIMATE DATA INTO clim STRUCTURE FOR CLEANESS
			return clim = CLIMATEDATA(Date, Pr, Pet, Temp, N_Climate, Pr_Through)

		end # function: CLIMATE


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θOBSERVATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct θOBSERVATION
			Date    :: Vector{DateTime}
			Z  	  :: Vector{Float64}
			iZobs   :: Vector{Int64}
			N_iT    :: Int64 # Number of time steps
			Ndepth  :: Int64 # Numver of soil profile with observed θ
			θobs 	  :: Array{Float64,2}
			∑T  	  :: Vector{Float64}
		end # mutable struct

		function TIME_SERIES()
		# Read data
			Data = DelimitedFiles.readdlm(path.calibr, ',')
		# Read header
			Header = Data[1,1:end]
		# Remove first READ_ROW_SELECT
			Data = Data[2:end,begin:end]

			Year, N_iT   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year")
			Month, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month")
			Day, ~     = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day")
			Hour, ~    = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour")
			Minute, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute")
			Second, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second")

			# READING THE DEPTH OF Θ MEASUREMENTS FROM HEADER: data having Z=
				θobs, Header = DelimitedFiles.readdlm(path.calibr, ','; header=true)

				Array_iHeader = []
				Ndepth = 0
				iCount = 0
				for iHeader in Header
					iCount += 1
					if occursin("Z=", iHeader) # Searching for 'Z=' in the header
						Ndepth += 1
						append!(Array_iHeader, iCount) 
					end # occursin
				end # iHeader

				# Isolating data with Z= measurements
				N_iT,~ = size(θobs)
				θobs = θobs[1:N_iT, minimum(Array_iHeader): maximum(Array_iHeader)]

				# The depths were we have θ measurements
				Z =  Vector{Float64}(undef, Ndepth)

				i = 0
				for iHeader in Header
					if occursin("Z=", iHeader)
						i += 1
						# Cleaning the header to get the integer
						iHeader = replace(iHeader, "Z=" => "")
						iHeader = replace(iHeader, "mm" => "")
						iHeader = replace(iHeader, " " => "")
						iHeader=  parse(Float64, iHeader)
						Z[i] = iHeader
					end # occursin("Z=", iHeader)
				end #  iHeader

			# REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
				Date_Start_Calibr = DateTime(param.hypix.calibr.Year_Start, param.hypix.calibr.Month_Start, param.hypix.calibr.Day_Start, param.hypix.calibr.Hour_Start, param.hypix.calibr.Minute_Start, param.hypix.calibr.Second_Start)
				
				Date_End_Calibr = DateTime(param.hypix.calibr.Year_End, param.hypix.calibr.Month_End, param.hypix.calibr.Day_End, param.hypix.calibr.Hour_End, param.hypix.calibr.Minute_End, param.hypix.calibr.Second_End)

			# ERROR CHECKING Assuring that Date_End ≤ Date_Clim_End
				Date_Clim_End = DateTime(param.hypix.Year_End, param.hypix.Month_End, param.hypix.Day_End, param.hypix.Hour_End, param.hypix.Minute_End, param.hypix.Second_End)

				Date_End_Calibr = min(Date_Clim_End, Date_End_Calibr)

			# SELECTING THE DATA WITHING FEASIBLE RANGE
				True = falses(N_iT) # Initiating with false
				Date = Array{DateTime}(undef, N_iT)
				iCount = 0 
				for iT=1:N_iT
					Date[iT] = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])

					if (Date_Start_Calibr ≤ Date[iT] ≤ Date_End_Calibr)
						iCount += 1
						True[iT] = true
					end  # if: 
				end # iT=1:N_Climate

				# New reduced number of simulations selected with dates
				Date = Date[True[1:N_iT]]
				θobs = θobs[True[1:N_iT],1:Ndepth]

				N_iT = iCount # New number of data

			# REDUCING THE AMOUNT OF DATA TO HOURLY
				if option.hypix.θobs_Hourly
					True = falses(N_iT)
					iCount = 0 
					for iT=1:N_iT
						if hour(Date[iT]) == 0 && minute(Date[iT]) == 0
							True[iT] = true
							iCount += 1
						end # if
					end # for
				
					# New reduced number of simulations selected with dates
					Date = Date[True[1:N_iT]]
					θobs = θobs[True[1:N_iT],1:Ndepth]

					N_iT = iCount # New reduced amount of data
				end # θobs_Hourly)

			# This will be computed at PrioProcess
				∑T    = Array{Float64}(undef, N_iT)
				iZobs = Array{Int64}(undef, Ndepth)

			return calibr = θOBSERVATION(Date, Z, iZobs, N_iT, Ndepth, θobs, ∑T)
		end  # function: TIME_SERIES

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LOOKUPTABLE
	#		Parameters as a function of time
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function LOOKUPTABLE_LAI(clim, veg)	
			if option.hypix.LookupTable_Lai == true
				LookUpTable_Lai, ~   = tool.readWrite.READ_HEADER(path.LookUpTable_Lai, "Lai")
			end
			
			i = 1
			Laiᵀ_Norm = Array{Float64}(undef, clim.N_Climate) 
			for Date in clim.Date
				Month = month(Date)
				if option.hypix.LookupTable_Lai == true
					Laiᵀ_Norm[i] = LookUpTable_Lai[Month]
				else
					Laiᵀ_Norm[i] = veg.Lai
				end
				i+=1
			end
		
		return Laiᵀ_Norm
		end  # function: LOOKUPTABLE_LAI


		function LOOKUPTABLE_CROPCOEFICIENT(clim, veg)
			if option.hypix.LookUpTable_CropCoeficient == true
				LookUpTable_CropCoeficient, ~   = tool.readWrite.READ_HEADER(path.LookUpTable_CropCoeficient, "CropCoeficient")
			end
			
			i = 1
			CropCoeficientᵀ_Norm = Array{Float64}(undef, clim.N_Climate) 
			for Date in clim.Date
				Month = month(Date)
				if option.hypix.LookUpTable_CropCoeficient == true
					CropCoeficientᵀ_Norm[i] = LookUpTable_CropCoeficient[Month]
				else
					CropCoeficientᵀ_Norm[i] = veg.CropCoeficient
				end
				i+=1
			end

		return CropCoeficientᵀ_Norm
		end  # function: LOOKUPTABLE_LAI


end  # module: readingHypix
# ............................................................