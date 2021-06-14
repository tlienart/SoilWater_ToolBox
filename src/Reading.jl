# =============================================================
#		MODULE: reading
# =============================================================
module reading
	import ..tool
	import  DelimitedFiles

	export ID, θΨ, KUNSATΨ, INFILTRATION, PSD, READ_STRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ID
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ID(;PathIdSlect, PathOptionSelect)
			println("    ~  $(PathIdSlect) ~")

			# Read data
				Data = DelimitedFiles.readdlm(PathIdSlect, ',')

				Header = Data[1,begin:end]

				Data = Data[2:end,begin:end]

				Data = sortslices(Data, dims=1)

				# Data, Header = DelimitedFiles.readdlm(Path, ',', Any, header=true)

				Id, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")
			
				Id_True, N_iZ_All = tool.readWrite.READ_HEADER_FAST(Data, Header, PathOptionSelect)

				Id = Int64.(Id)
			
				Id_True = Int64.(Id_True)

				Id_Select_True = convert.(Bool, Id_True)

			# Checking for errors
				for iZ=2:N_iZ_All
					if (Id[iZ] - Id[iZ-1]) < 1
						error("Id does not increase monotically at iD $(Id[iZ]) ")
					end
				end # for iZ=2:N_iZ_All
		
			N_SoilSelect = sum(Id_True)

			Id_Select = fill(0::Int64, N_SoilSelect)
			# Id_Select_True =  Array{Bool}(undef, N_iZ_All)
			iTrue = 1

			for iZ = 1:N_iZ_All
				if Id_True[iZ] == 1
					Id_Select[iTrue] = Id[iZ]
					iTrue += 1
				end	# Id_Soil == 1	
			end  # for: Id_Soil = Id_True
			
		return Id_Select, Id_Select_True, N_SoilSelect
		end  # function: ID



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BULKDENSITY_FINEEARTH_ROCKFRAGMENT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function BULKDENSITY_FINEEARTH_ROCKFRAGMENT(Id_Select_True, N_SoilSelect, Path)
		# 	# Load data
		# 		Data = JuliaDB.loadtable(Path)

		# 		# Getting the data
		# 			ρb_Fe_Rf = JuliaDB.select(Data, :BulkDensity_FineEarth_RockFragment_g_cm3)
		# 			ρb_Fe_Rf = ρb_Fe_Rf[Id_Select_True]

		# 			ρb_Fe = JuliaDB.select(Data,:BulkDensity_FineEarth_g_cm3)
		# 			ρb_Fe_Rf = ρb_Fe_Rf[Id_Select_True]

		# 			ρp_Fe = JuliaDB.select(Data,:ParticleDensity_FineEarth_g_cm3)
		# 			ρb_Fe_Rf = ρb_Fe_Rf[Id_Select_True]

		# 			ρ_Rf = JuliaDB.select(Data,:Density_RockFragment_g_cm3)
		# 			Rf_01 = JuliaDB.select(Data,:RockFragment_0_1)
		# 			Rf_Rlarge = JuliaDB.select(Data,:RockFragment_Rlarge_mm)
		# 			Rf_Rsmall = JuliaDB.select(Data,:RockFragment_Rsmall_mm)
		# 			Rf_Wetability = JuliaDB.select(Data, :RockFragment_Wettability_mm_mm2)

		# 	return
		# end  # function: name



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ρ_Ψθ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BULKDENSITY(Id_Select, N_SoilSelect, Path)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			ρbSoil, ~  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "BulkDensitySoil[g_cm-3]",  N_SoilSelect)
			
			ρp_Fine, ~ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "ParticleDensity_Fine[g_cm-3]",  N_SoilSelect)

			ρ_Rock, ~  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Density_Rock[g_cm-3]", N_SoilSelect)
			
			RockW, ~   = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Rock%", N_SoilSelect)

		return RockW, ρ_Rock, ρbSoil, ρp_Fine
		end # function: BulkDensity

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θψ_FILE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θψ_FILE(N_SoilSelect, N_θΨ::Int64, param, Path::String, θ_θΨ, Ψ_θΨ)
			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

				N_Ψ = Int64(length(param.hydro.Ψ_Table))

			# Writting the Header
				FieldName_String = Array{String}(undef, (N_Ψ))
				for iΨ =1:N_Ψ
					FieldName_String[iΨ] = string(Int64(param.hydro.Ψ_Table[iΨ]) ) * "mm"

					θobs, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, FieldName_String[iΨ])

					θ_θΨ =  [θobs[1:N_SoilSelect] θ_θΨ[1:N_SoilSelect,:] ]

					Ψ_Table = fill(Float64(param.hydro.Ψ_Table[iΨ]), N_SoilSelect)
				
					Ψ_θΨ = [Ψ_Table[1:N_SoilSelect] Ψ_θΨ[1:N_SoilSelect,:] ]
				end #for iΨ =1:N_Ψ

				for iZ=1:N_SoilSelect
					N_θΨ[iZ] += 2
				end
					
		return N_θΨ, θ_θΨ, Ψ_θΨ
		end  # function: θψ_FILE
	


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨ(Id_Select, N_SoilSelect, Path)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

			# Get the data of interest
				Ψ_θΨ, N_θΨ  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "H[mm]", N_SoilSelect)
			
				θ_θΨ, ~ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Theta", N_SoilSelect)

		return θ_θΨ, Ψ_θΨ, N_θΨ
		end  # function: θΨ



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSATΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KUNSATΨ(Id_Select, N_SoilSelect, path)
			# Determeining where to read the data
			if isfile(path.inputSoilwater.Kunsat)
				Path = path.inputSoilwater.Kunsat
			elseif isfile(path.inputSoilwater.Kunsat_Model)
				Path = path.inputSoilwater.Kunsat_Model
			else
				error("\n SoilWater-ToolBox input error: No Kunsat data. You coud derive K(θ) from Kosugi model with option.UsePointKosugiBimodal \n")
			end
			
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')

			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)
			
			Ψ_KΨ, N_KΨ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "H[mm]", N_SoilSelect)
			
         K_KΨ, ~    = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Kunsat[mm_s]", N_SoilSelect)
			
		return K_KΨ, Ψ_KΨ, N_KΨ 
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : READ_STRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function READ_STRUCT(structures, Path; iStart=1, iEnd=2^63 - 1)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,1:end]
			# Select data of interest
				N_SoilSelect = size(Data)[1] # Initial
				iEnd= min(N_SoilSelect, iEnd)
				Data = Data[iStart:iEnd,1:end]
				N_SoilSelect = iEnd - iStart + 1 # Final

			# Reading the Model data
			for iFieldname in propertynames(structures)

				# Putting the values of Output_Vector into structures					
				try
					Output_Vector, Ndata = tool.readWrite.READ_HEADER_FAST(Data, Header, string(iFieldname))

					# Putting the values of OutPutData into the the structures hydroData
						if length(Output_Vector) == 1 # Not a vector
							Output_Vector = Output_Vector[1]
						end

						setfield!(structures, Symbol(iFieldname), Float64.(Output_Vector))
				catch
					# @warn "SoilWater-ToolBox: cannong find $iFieldname"
					Output_Vector = fill(0.0::Float64, N_SoilSelect)
					if length(Output_Vector) == 1 # Not a vector
                  Output_Vector = Output_Vector[1]
               end
	
					setfield!(structures, Symbol(iFieldname), Float64.(Output_Vector))
				end
			end

		return structures, N_SoilSelect
		end  # function: READ_STRUCT



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct INFILT
			RingRadius
			θ_Ini
			γ
			β
		end # struct INFILT

		function INFILTRATION(Id_Select, N_SoilSelect, PathInfilt, PathInfiltParam)
			println("    ~  $(PathInfilt) ~")

			# Read data
				Data = DelimitedFiles.readdlm(PathInfilt, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			# Reading select data
				Tinfilt, N_Infilt = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Tinfilt[s]", N_SoilSelect)
				
				∑Infilt_Obs , ~ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Cumul_Infiltration[mm]", N_SoilSelect)
				
			println("    ~  $(PathInfiltParam) ~")

			# Read data
				Data = DelimitedFiles.readdlm(PathInfiltParam, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

			RingRadius , ~  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "RingRadius[mm]", N_SoilSelect; N_Point_Max=1)

			θ_Ini , ~       = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Theta_Ini[-]", N_SoilSelect; N_Point_Max=1)

			γ , ~           = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Lambda[-]", N_SoilSelect; N_Point_Max=1)

			β , ~           = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Beta[-]", N_SoilSelect; N_Point_Max=1)

			infiltParam = INFILT(RingRadius, θ_Ini, γ, β)

		return Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam
		end  # function: INFILTRATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PSD(Id_Select, N_SoilSelect, Path) # TODO make sure that the particles are ordered from smalest to largest
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			Diameter_Psd, N_Psd 	= tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Diameter[mm]",  N_SoilSelect)

			∑Psd , ~ 				= tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Cumul_Psd", N_SoilSelect)

			Rpart = @. Diameter_Psd / 2.0

		return Rpart, ∑Psd, N_Psd
		end  # function: PSD



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct OPTIM
			Param_Name :: Vector{String}
			ParamOpt_Min :: Vector{Float64}
			ParamOpt_Max :: Vector{Float64}
			Param_Min :: Vector{Float64}
			Param_Max :: Vector{Float64}
			ParamOpt :: Vector{String}
			NparamOpt :: Int64
			Flag_Opt :: Bool
			ParamOpt_LogTransform :: Vector{Bool}
		end

		function HYDRO_PARAM(optionₘ, hydro, N_SoilSelect, Path)
		# Read data
			Data = DelimitedFiles.readdlm(Path, ',')
		# Read header
			Header = Data[1,1:end]
		# Remove first READ_ROW_SELECT
			Data = Data[2:end,begin:end]

		# Reading the Model data
			HydroModel, Ndata   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MODEL")

		# Determening which parameters correspond to the selected model
		iSelectModel = [] 
		for i=1:Ndata
			if HydroModel[i] == string(optionₘ.HydroModel)
				append!(iSelectModel, i)
			end
		end

		# Reading the names of the parameters
			Param_Name, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "ParamName")
				# Selecing data
				Param_Name = Param_Name[iSelectModel]

		# Reading minimum value of the parameters
			Param_Min, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MIN")
				# Selecing data
				Param_Min = Param_Min[iSelectModel]

		# Reading maximum value of the parameters
			Param_Max, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MAX")
				# Selecing data
				Param_Max= Param_Max[iSelectModel]

		# Reading parameters requires log transformation [1 or 0]
			Opt_LogTransform, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "LogTransform")
				# Selecing data
				Opt_LogTransform= Opt_LogTransform[iSelectModel]

		# Reading the values of the parameters if they are not optimized
			ParamValue, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "VALUE")
				# Selecing data
				ParamValue = ParamValue[iSelectModel]

		# Reading which parameters to be optimized [1 or 0]
			Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT")
			# Selecing data
			Opt = Opt[iSelectModel]
			
		# Determine if we need to optimize
			if sum(Opt) ≥ 1
				Flag_Opt = true
			else
				Flag_Opt = false
			end

		# ====================================================
		ParamOpt              = []
		ParamOpt_Max          = []
		ParamOpt_Min          = []
		ParamOpt_LogTransform = []

		i = 1
		# For every hydraulic parameter
		for inParamValue in Param_Name
			# Putting the value of the parameters in hydro. Repeating the value of the parameter for all soils data: N_SoilSelect
			ParamValue_Vector = fill(Float64(ParamValue[i]), N_SoilSelect)
			setfield!(hydro, Symbol(inParamValue), ParamValue_Vector)

			# θsMacMat value depends on θs
			if  Symbol(inParamValue) == :θsMacMat_ƞ
				for iZ = 1:N_SoilSelect 
					hydro.θsMacMat[iZ] =  hydro.θs[iZ] * hydro.θsMacMat_ƞ[iZ]
				end
			end # Symbol(inParamValue) == :θsMacMat_ƞ

			# Putting the minimum value in the parameter
				ParamValue_Vector = fill(Float64(Param_Min[i]), N_SoilSelect)
				setfield!(hydro, Symbol(inParamValue * "_Min"), ParamValue_Vector)

			# Putting the maximum value in the parameter
				ParamValue_Vector = fill(Float64(Param_Max[i]), N_SoilSelect)
				setfield!(hydro, Symbol(inParamValue * "_Max"), ParamValue_Vector)
	
			# ParamValue to optimize. The complication is that there are different layers of hydraulic parameters which can be optimized.  
			if Opt[i] == 1

				# appending the values of the parameters
				append!(ParamOpt, [Param_Name[i]])

				append!(ParamOpt_Min, Param_Min[i])
				
				append!(ParamOpt_Max, Param_Max[i])

				# Appending name of param to perform logTransform if optimized
				if Opt_LogTransform[i] == 1
					append!(ParamOpt_LogTransform, [true])
				else
					append!(ParamOpt_LogTransform, [false])
				end

				if Param_Min[i] > Param_Max[i]
					error("LabOpt ERROR: $(Param_Min[i]) < $(ParamValue[i]) < $(Param_Max[i]) !")
				end
			end # if Flag_Opt

			i += 1
		end # for loop

		# Number of parameters to be optimised
			NparamOpt = length(ParamOpt)
	
		# Putting all the in mutable structure
			optim = OPTIM(Param_Name,ParamOpt_Min,ParamOpt_Max,Param_Min,Param_Max,ParamOpt,NparamOpt,Flag_Opt,ParamOpt_LogTransform)

		if Flag_Opt == true
			println("		=== === Optimizing the following parameters === ===")
			println("			Model=" , optionₘ.HydroModel)
			println("			NparamOpt=" , NparamOpt)
			println("			ParamOpt= " ,  optim.ParamOpt)
			println("			Min_Value= " , optim.ParamOpt_Min)
			println("			Max_Value= " , optim.ParamOpt_Max)
			println("			LogTransform = " , optim.ParamOpt_LogTransform)
			println("		=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n \n")
		end

	return hydro, optim
	end  # function: HydroParam_ThetaH


	# =============================================================
	#		module: reading Hypix
	# =============================================================
	module hyPix
		import  ...tool, ..horizonLayer
		import Dates: value, DateTime, hour, minute, month, now
		import DelimitedFiles
		export CLIMATE, DISCRETIZATION, HYPIX_PARAM, LOOKUPTABLE_LAI, LOOKUPTABLE_CROPCOEFICIENT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : DATES
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function DATES(param, pathHyPix)
				# Read data
					Data = DelimitedFiles.readdlm(pathHyPix.Dates, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]

				# Dates of climate data
					Year_Start₀ , ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year_Sim_Start")
						param.hyPix.Year_Start = Year_Start₀[1]
					Month_Start₀, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month_Sim_Start")
							param.hyPix.Month_Start = Month_Start₀[1]
					Day_Start₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day_Sim_Start")
						param.hyPix.Day_Start = Day_Start₀[1]

					Year_End₀, ~    = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year_Sim_End")
						param.hyPix.Year_End = Year_End₀[1]
					Month_End₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month_Sim_End")
						param.hyPix.Month_End = Month_End₀[1]
					Day_End₀, ~     = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day_Sim_End")
						param.hyPix.Day_End = Day_End₀[1]
					
				# Dates of observed data
					param.hyPix.obsθ.Year_Start  = param.hyPix.Year_Start
					param.hyPix.obsθ.Month_Start = param.hyPix.Month_Start
					param.hyPix.obsθ.Day_Start   = param.hyPix.Day_Start

					param.hyPix.obsθ.Year_End    = param.hyPix.Year_End
					param.hyPix.obsθ.Month_End   = param.hyPix.Month_End
					param.hyPix.obsθ.Day_End     = param.hyPix.Day_End

				# Dates of plots
				param.hyPix.ploting.Year_Start  = param.hyPix.Year_Start
				param.hyPix.ploting.Month_Start = param.hyPix.Month_Start
				param.hyPix.ploting.Day_Start   = param.hyPix.Day_Start

				param.hyPix.ploting.Year_End    = param.hyPix.Year_End
				param.hyPix.ploting.Month_End   = param.hyPix.Month_End
				param.hyPix.ploting.Day_End     = param.hyPix.Day_End

			return param
			end  # function: DATES

			
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : DISCRETIZATION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function DISCRETIZATION(pathHyPix)
				# Read data
					Data = DelimitedFiles.readdlm(pathHyPix.Discretization, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]

				Z, N_iZ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Z")
				θ_Ini, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "θini")
				Layer, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Layer")

				N_iHorizon = maximum(Layer)
			return Layer, N_iHorizon, N_iZ, Z, θ_Ini
			end # function DISCRETIZATION


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYPIX_PARAM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYPIX_PARAM(Layer, hydro, hydroHorizon, iSim::Int64, N_iZ::Int64, veg)
				# Read data
					Data = DelimitedFiles.readdlm(pathHyPix.Hypix_Param, ',')
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
					if Flag_Opt && (iSim ≥ param.hyPix.iSim_Start + 1)
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

				if !(Flag_MultiStepOpt)
					# Hydraulic parameters per horizon to layers
					hydro = horizonLayer.HYDROHORIZON_2_HYDRO(N_iZ, Layer, hydroHorizon)
				end

				NparamOpt = length(ParamOpt)

				# CHECKING FOR UNCONSISTENCY WITH OPTIONS	
				if Flag_Opt && option.hyPix.σ_2_Ψm ≠ :No && "Ψm" ∈ ParamOpt
					iψm = findfirst(isequal("Ψm"), ParamOpt)[1]

					if option.hyPix.σ_2_Ψm==:UniqueRelationship && "Ψm" ∈ ParamOpt
						error( "**** HyPix Error: combination of options which are not possible (option.hyPix.σ_2_Ψm==:UniqueRelationship) && (Optimise=Ψm)!")

					elseif option.hyPix.σ_2_Ψm==:Constrained && !("Ψm" ∈ ParamOpt)
						error("*** HyPix Error: combination of options which are not possible (option.hyPix.σ_2_Ψm==:Constrained) && (not Optimising=Ψm)!")

					elseif option.hyPix.σ_2_Ψm==:Constrained && ParamOpt_LogTransform[iψm]==1
						error("*** option.hyPix.σ_2_Ψm==:Constrained CANNOT log transforme Ψm") 
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

			Option_ReadTemperature = false

			function CLIMATE(option, param, pathHyPix)
				if option.hyPix.ClimateDataTimestep == "Daily"
					Pr_Name          = "Rain(mm)"
					Pet_Name         = "PET(mm)"
					Temperature_Name = "Tmax(C)" # Maximum temperature which is not correct

				elseif option.hyPix.ClimateDataTimestep == "Hourly"
					Pr_Name          = "Pr_mm"
					Pet_Name         = "Pet_mm"
					Temperature_Name = "Temp_c"
				end #  option.hyPix.ClimateDataTimestep

				# READ DATA
					Data = DelimitedFiles.readdlm(pathHyPix.Climate, ',')
					Header = Data[1,1:end]
					Data = Data[2:end,begin:end]

					Year, N_Climate = tool.readWrite.READ_HEADER_FAST(Data, Header,"Year")
					Month, ~        = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month")
					Day, ~          = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day")
					Hour, ~         = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour")
					Minute, ~       = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute")
					Second, ~       = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second")
					Pr, ~           = tool.readWrite.READ_HEADER_FAST(Data, Header, Pr_Name)
					Pet, ~          = tool.readWrite.READ_HEADER_FAST(Data, Header, Pet_Name)
					if Option_ReadTemperature 
						Temp, ~         = tool.readWrite.READ_HEADER_FAST(Data, Header, Temperature_Name)
					else
						Temp = fill(24.0::Float64, N_Climate)
					end

				# READING DATES FROM FILE
					param = DATES(param, pathHyPix)

				# REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
					Date_Start = DateTime(param.hyPix.Year_Start, param.hyPix.Month_Start, param.hyPix.Day_Start, param.hyPix.Hour_Start, param.hyPix.Minute_Start, param.hyPix.Second_Start)
					
					Date_End = DateTime(param.hyPix.Year_End, param.hyPix.Month_End, param.hyPix.Day_End, param.hyPix.Hour_End, param.hyPix.Minute_End, param.hyPix.Second_End)

				# CHECKING
					# End Date feasible
						Date_End_Maximum = DateTime(Year[N_Climate], Month[N_Climate], Day[N_Climate], Hour[N_Climate], Minute[N_Climate], Second[N_Climate]) 

						if Date_End_Maximum < Date_End
							Date_End = min(Date_End_Maximum, Date_End)
							println("		~ HyPix WARNING: Date_End not feasible so modified to match the data ")
						end #warning Date_End

					# Start Date feasible
						Date_Start_Minimum = DateTime(Year[2], Month[2], Day[2], Hour[2], Minute[2], Second[2]) 

						if Date_Start_Minimum > Date_Start
							Date_Start = max(Date_Start_Minimum , Date_Start)
							println("		~HyPix WARNING: Date_Start = $Date_Start not feasible so modified to match the data ")
						end #warning DATE_START

				# SELECTING DATES OF INTEREST
					True = falses(N_Climate)
					Date = fill(now()::DateTime,  N_Climate) 
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

			# STRUCTURE
				clim = CLIMATEDATA(Date, Pr, Pet, Temp, N_Climate, Pr_Through)

			# SAVING SPACE 
				Data = nothing
				True = nothing

			return clim
			end # function: CLIMATE


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θOBSERVATION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			mutable struct θOBSERVATION
				Date    :: Vector{DateTime}
				Z  	  :: Vector{Float64}
				ithetaObs   :: Vector{Int64}
				N_iT    :: Int64 # Number of time steps
				Ndepth  :: Int64 # Numver of soil profile with observed θ
				θobs 	  :: Array{Float64,2}
				∑T  	  :: Vector{Float64}
			end # mutable struct

			function TIME_SERIES(option, param, pathHyPix)
			# Read data
				Data = DelimitedFiles.readdlm(pathHyPix.obsθ, ',')
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
					θobs, Header = DelimitedFiles.readdlm(pathHyPix.obsθ, ','; header=true)

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
					Z = fill(0.0::Float64, Ndepth)

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

				# READING DATES FROM FILE
					param = DATES(param, pathHyPix)

				# REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
					Date_Start_Calibr = DateTime(param.hyPix.obsθ.Year_Start, param.hyPix.obsθ.Month_Start, param.hyPix.obsθ.Day_Start, param.hyPix.obsθ.Hour_Start, param.hyPix.obsθ.Minute_Start, param.hyPix.obsθ.Second_Start)
					
					Date_End_Calibr = DateTime(param.hyPix.obsθ.Year_End, param.hyPix.obsθ.Month_End, param.hyPix.obsθ.Day_End, param.hyPix.obsθ.Hour_End, param.hyPix.obsθ.Minute_End, param.hyPix.obsθ.Second_End)

				# ERROR CHECKING Assuring that Date_End ≤ Date_Clim_End
					Date_Clim_End = DateTime(param.hyPix.Year_End, param.hyPix.Month_End, param.hyPix.Day_End, param.hyPix.Hour_End, param.hyPix.Minute_End, param.hyPix.Second_End)

					Date_End_Calibr = min(Date_Clim_End, Date_End_Calibr)

				# SELECTING THE DATA WITHING FEASIBLE RANGE
					True = falses(N_iT) # Initiating with false
					Date = fill(now()::DateTime,  N_iT)
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
					ΔTimeStep = value(Date[5]-Date[4])/ 1000
					if option.hyPix.θobs_Hourly && ΔTimeStep < 86400
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
					∑T        = fill(0.0::Float64, N_iT)
					ithetaObs = fill(0::Int64, Ndepth)

				# STRUCTURE
					obsθ = θOBSERVATION(Date, Z, ithetaObs, N_iT, Ndepth, θobs, ∑T)

				# SAVING SPACE 
					Data = nothing
					True = nothing
			return obsθ
			end  # function: TIME_SERIES


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : LOOKUPTABLE
		#		Parameters as a function of time
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function LOOKUPTABLE_LAI(clim, option, pathHyPix, veg)	
				if option.hyPix.LookupTable_Lai == true
					LookUpTable_Lai, ~   = tool.readWrite.READ_HEADER(pathHyPix.LookUpTable_Lai, "Lai")
				end
				
				i = 1
				Laiᵀ_Norm = fill(0.0::Float64, clim.N_Climate) 
				for Date in clim.Date
					Month = month(Date)
					if option.hyPix.LookupTable_Lai == true
						Laiᵀ_Norm[i] = LookUpTable_Lai[Month]
					else
						Laiᵀ_Norm[i] = veg.Lai
					end
					i+=1
				end
			
			return Laiᵀ_Norm
			end  # function: LOOKUPTABLE_LAI


			function LOOKUPTABLE_CROPCOEFICIENT(clim, option, pathHyPix, veg)
				if option.hyPix.LookUpTable_CropCoeficient == true
					LookUpTable_CropCoeficient, ~   = tool.readWrite.READ_HEADER(pathHyPix.LookUpTable_CropCoeficient, "CropCoeficient")
				end
				
				i = 1
				CropCoeficientᵀ_Norm = fill(0.0::Float64, clim.N_Climate) 
				for Date in clim.Date
					Month = month(Date)
					if option.hyPix.LookUpTable_CropCoeficient == true
						CropCoeficientᵀ_Norm[i] = LookUpTable_CropCoeficient[Month]
					else
						CropCoeficientᵀ_Norm[i] = veg.CropCoeficient
					end
					i+=1
				end

			return CropCoeficientᵀ_Norm
			end  # function: LOOKUPTABLE_LAI

	end  # module: readingHypix
	# .........................................................

end  # module: reading
# ............................................................		