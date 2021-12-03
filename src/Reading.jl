# =============================================================
#		MODULE: reading
# =============================================================
module reading
	import ..tool, ..table
	import  DelimitedFiles
	export ID, θΨ, KUNSATΨ, INFILTRATION, PSD, READ_STRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ID
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ID(;PathIdSelect, PathOptionSelect, PathModelName)
			println("    ~  $(PathIdSelect) ~")

			# Read data
				Data = DelimitedFiles.readdlm(PathIdSelect, ',')
				Header = Data[1,begin:end]
				Data = Data[2:end,begin:end]
				Data = sortslices(Data, dims=1)

				Id, N_iZ_All  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")
			
				Id_True, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, PathOptionSelect)

				Id = Int64.(Id)

			# Soilname is optional
				Soilname = []
				try
					Soilname, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")
				catch # If not available
					Soilname = fill("", N_iZ_All)
					for i=1:N_iZ_All
						Soilname[i] = PathModelName  * "_" * string(Id[i])
					end
				end
		
				Id_True = Int64.(Id_True)

				IdSelect_True = convert.(Bool, Id_True)

			# Checking for errors
				for iZ=2:N_iZ_All
					if (Id[iZ] - Id[iZ-1]) < 1
						error("Id does not increase monotically at iD $(Id[iZ]) ")
					end
				end # for iZ=2:N_iZ_All
		
			NiZ = sum(Id_True)

			IdSelect = Id[IdSelect_True]
			Soilname = Soilname[IdSelect_True]
	
		return IdSelect, IdSelect_True, Soilname, NiZ
		end  # function: ID


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : bulk density
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BULKDENSITY(IdSelect, NiZ, Path)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			ρᵦ_Soil, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "BulkDensitySoil[g_cm-3]",  NiZ, N_Point_Max=1)

			ρₚ_Fine, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "ParticleDensity_Fine[g_cm-3]",  NiZ, N_Point_Max=1)

			ρₚ_Rock, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Density_Rock[g_cm-3]", NiZ, N_Point_Max=1)
			
			RockFragment, ~   = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header,"RockFragment[0-1]", NiZ, N_Point_Max=1)
		return RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil
		end # function: BulkDensity


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Base.@kwdef mutable struct INFILT
			RingRadius
			θini
			γ
			β
		end # struct INFILT

		function INFILTRATION(IdSelect, NiZ, PathInfilt, PathInfiltParam)
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
				Tinfilt, N_Infilt = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Tinfilt[s]", NiZ)
				
				∑Infilt_Obs , ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Cumul_Infiltration[mm]", NiZ)
				
			#-----------------------------------------------------------------------
			println("    ~  $(PathInfiltParam) ~")

			# Read data
				Data = DelimitedFiles.readdlm(PathInfiltParam, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

			RingRadius , ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "RingRadius[mm]", NiZ; N_Point_Max=1)

			θini , ~       = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header,"Theta_Ini[-]", NiZ; N_Point_Max=1)

			γ , ~           = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Lambda[-]", NiZ; N_Point_Max=1)

			β , ~           = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Beta[-]", NiZ; N_Point_Max=1)

			infiltParam = INFILT(RingRadius, θini, γ, β)
		return Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam
		end  # function: INFILTRATION


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨ(IdSelect, NiZ, path)
				println("    ~  $(path.inputSoilwater.Ψθ) ~")

				# Read data
					Data = DelimitedFiles.readdlm(path.inputSoilwater.Ψθ, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first row
					Data = Data[2:end,begin:end]
				# Sort data
					Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

				# Determeining if data has only 3 columns: Id, H and Theta
				if length(Header) ≤ 6
					# Get the data of interest
						Ψ_θΨobs, N_θΨobs  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "H[mm]", NiZ)
				
						θ_θΨobs, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Theta[0-1]", NiZ)
				
				# Data is in square [X=iZ, Y =iΨ]
				else
					N_θΨobs, θ_θΨobs, Ψ_θΨobs = tool.readWrite.READ_θΨK_2D(Data, Header, IdSelect, NiZ)

					table.convert.CONVERT_θΨ_2D_2_1D(IdSelect, NiZ, N_θΨobs, path.convertSoilwater.Table_Convert_θΨ_2D_2_1D, θ_θΨobs, Ψ_θΨobs)
				end # length(Header) == 3

			return θ_θΨobs, Ψ_θΨobs, N_θΨobs
			end  # function: θΨ
		#----------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KUNSATΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KUNSATΨ(IdSelect, NiZ, path, Path)
				println("    ~  $(Path) ~")

				# Read data
					Data = DelimitedFiles.readdlm(Path, ',')

				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]
				# Sort data
					Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

				# Determeining if data has only 3 columns: Id, H and Theta
				if length(Header) == 3
					Ψ_KΨobs, N_KΨobs = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "H[mm]", NiZ)
						
					K_KΨobs, ~    = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header,"Kunsat[mm_s]", NiZ)
				# Data is in square [X=iZ, Y =iΨ]
				else
					N_KΨobs, K_KΨobs, Ψ_KΨobs = tool.readWrite.READ_θΨK_2D(Data, Header, IdSelect, NiZ)

					table.convert.CONVERT_KΨ_2D_2_1D(IdSelect, NiZ, N_KΨobs, path.convertSoilwater.Table_Convert_KΨ_2D_2_1D, K_KΨobs, Ψ_KΨobs)
				end
			return K_KΨobs, Ψ_KΨobs, N_KΨobs 
			end  # function: θΨ
		#----------------------------------------------------------------------

		
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PSD(IdSelect, NiZ, Path) # TODO make sure that the particles are ordered from smalest to largest
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first RockWetability
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			Diameter_Psd, N_Psd = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Diameter[mm]", NiZ)

			∑Psd , ~            = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "Cumul_Psd", NiZ)

			Rpart = @. Diameter_Psd / 2.0
		return Rpart, ∑Psd, N_Psd
		end  # function: PSD
	#----------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :SOIL_INOFRMATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PEDOLOGICAL(IdSelect, NiZ, Path)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)
			
			IsTopsoil, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "IsTopsoil", NiZ, N_Point_Max=1)
			
			RockClass, ~ = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "RockClass", NiZ, N_Point_Max=1)
		return IsTopsoil, RockClass
		end # function: SOIL_INOFRMATION
	#----------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : bulk density
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Φ(IdSelect, NiZ, Path)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			Φ, ~  = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header, "TotalPorosity[0-1]", NiZ, N_Point_Max=1)
			
			RockFragment, ~   = tool.readWrite.READ_ROW_SELECT(IdSelect, Data, Header,"RockFragment[0-1]", NiZ, N_Point_Max=1)
		return RockFragment, Φ
		end # function: BulkDensity
	#----------------------------------------------------------------------

		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θψ_ADDPOINTS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θψ_ADDPOINTS(NiZ, N_θΨobs::Int64, param, Path::String, θ_θΨobs, Ψ_θΨobs)
			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

				N_Ψ = Int64(length(param.hydro.Ψ_Table))

			# Writting the Header
				FieldName_String = fill(""::String, (N_Ψ))
				for iΨ =1:N_Ψ
					FieldName_String[iΨ] = string(Int64(param.hydro.Ψ_Table[iΨ]) ) * "mm"

					θobs, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, FieldName_String[iΨ])

					θ_θΨobs =  [θobs[1:NiZ] θ_θΨobs[1:NiZ,:] ]

					Ψ_Table = fill(Float64(param.hydro.Ψ_Table[iΨ]), NiZ)
				
					Ψ_θΨobs = [Ψ_Table[1:NiZ] Ψ_θΨobs[1:NiZ,:] ]
				end #for iΨ =1:N_Ψ

				for iZ=1:NiZ
					N_θΨobs[iZ] += 2
				end		
		return N_θΨobs, θ_θΨobs, Ψ_θΨobs
		end  # function: θψ_ADDPOINTS+
	#----------------------------------------------------------------------


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
				NiZ = size(Data)[1] # Initial
				iEnd= min(NiZ, iEnd)
				Data = Data[iStart:iEnd,1:end]
				NiZ = iEnd - iStart + 1 # Final

			# Reading the Model data
			for iFieldname in propertynames(structures)

				# Putting the values of Output_Vector into structures
				Output_Vector = fill(0.0::Float64, NiZ)					
				try
					Output_Vector, Ndata = tool.readWrite.READ_HEADER_FAST(Data, Header, string(iFieldname))
				catch
					# @warn "SoilWater-ToolBox: cannong find $iFieldname"
					Output_Vector = fill(0.0::Float64, NiZ)
				end

				try
					setfield!(structures, Symbol(iFieldname), Float64.(Output_Vector))
				catch
					setfield!(structures, Symbol(iFieldname), Float64(Output_Vector[1]))
				end
			end

		return structures, NiZ
		end  # function: READ_STRUCT
	#----------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Base.@kwdef mutable struct OPTIM
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

		function HYDRO_PARAM(optionₘ, hydro, NiZ, Path)
		# Read data
			Data = DelimitedFiles.readdlm(Path, ',')
		# Read header
			Header = Data[1,1:end]
		# Remove first READ_ROW_SELECT
			Data = Data[2:end,begin:end]

		# Reading the Model data
			HydroModel⍰, Ndata   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MODEL")

		# Determening which parameters correspond to the selected model
		iSelectModel = [] 
		for i=1:Ndata
			if HydroModel⍰[i] == string(optionₘ.HydroModel⍰)
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
			# Putting the value of the parameters in hydro. Repeating the value of the parameter for all soils data: NiZ
			ParamValue_Vector = fill(Float64(ParamValue[i]), NiZ)
			setfield!(hydro, Symbol(inParamValue), ParamValue_Vector)

			# θsMacMat value depends on θs
			if  Symbol(inParamValue) == :θsMacMat_ƞ
				for iZ = 1:NiZ 
					hydro.θsMacMat[iZ] =  hydro.θs[iZ] * hydro.θsMacMat_ƞ[iZ]
				end
			end # Symbol(inParamValue) == :θsMacMat_ƞ

			# Putting the minimum value in the parameter
				ParamValue_Vector = fill(Float64(Param_Min[i]), NiZ)
				setfield!(hydro, Symbol(inParamValue * "_Min"), ParamValue_Vector)

			# Putting the maximum value in the parameter
				ParamValue_Vector = fill(Float64(Param_Max[i]), NiZ)
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
			println("	=== === Optimizing the following parameters === ===")
			println("		Model=" , optionₘ.HydroModel⍰)
			println("		NparamOpt=" , NparamOpt)
			println("		ParamOpt= " ,  optim.ParamOpt)
			println("		Min_Value= " , optim.ParamOpt_Min)
			println("		Max_Value= " , optim.ParamOpt_Max)
			println("		LogTransform = " , optim.ParamOpt_LogTransform)
			println("	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n")
		end

	return hydro, optim
	end  # function: GUI_HydroParam


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  KSMODEL_PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Base.@kwdef mutable struct OPTIMKS
         Param_Name   :: Array{String}
         ParamOpt_Min :: Array{Float64}
         ParamOpt_Max :: Array{Float64}
         ParamOpt     :: Array{String}
         NparamOpt    :: Vector{Int64}
         Flag_Opt     :: Bool
		end

		function KSMODEL_PARAM(ksmodelτ, option, Path)
			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Reading Model data
				 KₛModel⍰, Ndata = tool.readWrite.READ_HEADER_FAST(Data, Header, "MODEL")

			# Determening which parameters correspond to the selected model
			iSelectModel = [] 
			for i=1:Ndata
				if  KₛModel⍰[i] == string(option.ksModel. KₛModel⍰)
					append!(iSelectModel, i)
				end
			end # for i=1:Ndata

			# Reading names of the parameters
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
					Param_Max = Param_Max[iSelectModel]

			# Reading values of the default values of the parameters
				ParamValue, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "VALUE")
					# Selecing data
					ParamValue = ParamValue[iSelectModel]

			# Reading which parameters to be optimized [1 or 0]
				Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT")
				# Selecing data
				Opt = Opt[iSelectModel]
				
			N_Opt = sum(Opt)
			# Determine if we need to optimize
				if N_Opt ≥ 1
					Flag_Opt = true
				else
					Flag_Opt = false
				end

			# ====================================================
            ParamOpt     = fill(""::String, (2, N_Opt))
            ParamOpt_Min = fill(0.0::Float64, (2, N_Opt))
            ParamOpt_Max = fill(0.0::Float64, (2, N_Opt))
            NparamOpt    = fill(0::Int64, 2)

			# For every parameter
			# The τ[1] = TopLayer and   τ[3] = SubLayer			
			i = 1
			for ipParamName in Param_Name
				if occursin("_1", ipParamName)
					ipParamName = replace(ipParamName, "_1" => "" )
					iGroup = 1
					Param_Name[i] = ipParamName
					Flag_Top = true
					
				elseif occursin("_2", ipParamName)
					ipParamName = replace(ipParamName, "_2" => "" )
					iGroup = 2
					Param_Name[i] = ipParamName
					Flag_Top = false
				end

				# Getting the Vector values of the τ parameters
					ParamValue_Vector = getfield(ksmodelτ, Symbol(ipParamName))
					ParamValue_Vector[iGroup] = Float64(ParamValue[i])
					# Storing the value
					setfield!(ksmodelτ, Symbol(ipParamName), ParamValue_Vector)

				# Putting the minimum value in the parameter
					ParamValue_Vector = getfield(ksmodelτ, Symbol(ipParamName * "_Min"))
					ParamValue_Vector[iGroup] = Float64(Param_Min[i])
					# Storing the value
					setfield!(ksmodelτ, Symbol(ipParamName * "_Min"), ParamValue_Vector)

				# Putting the maximum value in the parameter
					ParamValue_Vector = getfield(ksmodelτ, Symbol(ipParamName * "_Max"))
					ParamValue_Vector[iGroup] = Float64(Param_Max[i])
					# Storing the value
					setfield!(ksmodelτ, Symbol(ipParamName * "_Max"), ParamValue_Vector)

				# ParamValue to optimize.  
				if Opt[i] == 1
					NparamOpt[iGroup] += 1

					if Flag_Top
                  ParamOpt[iGroup, NparamOpt[iGroup]]     = Param_Name[i]
                  ParamOpt_Min[iGroup, NparamOpt[iGroup]] = Param_Min[i]
                  ParamOpt_Max[iGroup, NparamOpt[iGroup]] = Param_Max[i]
					else
                  ParamOpt[iGroup, NparamOpt[iGroup]]     = Param_Name[i]
                  ParamOpt_Min[iGroup, NparamOpt[iGroup]] = Param_Min[i]
                  ParamOpt_Max[iGroup, NparamOpt[iGroup]] = Param_Max[i]
					end

					# Checking error
						if ParamOpt_Min[iGroup, NparamOpt[iGroup]] > ParamOpt_Max[iGroup, NparamOpt[iGroup]]
							error("SoilWater LabOpt ERROR: $(ParamOpt[iGroup, NparamOpt[iGroup]]) $(ParamOpt_Min[iGroup, NparamOpt[iGroup]] ) < $(ParamValue[i]) < $( ParamOpt_Max[iGroup, NparamOpt[iGroup]]) !")
						end
				end # if Flag_Opt
			i += 1
			end # for loop

			# Putting all the in mutable structure
				optimKsmodel = OPTIMKS(Param_Name, ParamOpt_Min, ParamOpt_Max, ParamOpt, NparamOpt, Flag_Opt)

			if Flag_Opt == true
				println("	=== === Optimizing the following τ parameters === === \n")
				println("		KsModel=" , option.ksModel. KₛModel⍰)
				println("		ksmodelτ=", Param_Name)
				# println("		ksmodelτ=", ksmodelτ)
				println("		NparamOpt_τ=" , optimKsmodel.NparamOpt)
				println("		ParamOpt_τ= " ,  optimKsmodel.ParamOpt)
				println("		Min_Value_τ= " , optimKsmodel.ParamOpt_Min)
				println("		Max_Value_τ = " , optimKsmodel.ParamOpt_Max)
				println("	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n")
			end
	return ksmodelτ, optimKsmodel
	end  # function: KSMODEL_PARAM
	# ............................................................


	# =============================================================
	#		module: reading Hypix
	# =============================================================
	module hyPix
		import  ...tool, ...horizonLayer
		import Dates: value, DateTime, hour, minute, month, now, Hour
		import DelimitedFiles
		export CLIMATE, DATES, DISCRETIZATION, HYPIX_PARAM, LOOKUPTABLE_CROPCOEFICIENT, LOOKUPTABLE_LAI

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

				# Dates of observed climate data
					Year_StartObs₀ , ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year_Obs_Start")
						param.hyPix.Year_Start = Int64(Year_StartObs₀[1])

					Month_StartObs₀, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month_Obs_Start")
						param.hyPix.Month_Start = Int64(Month_StartObs₀[1])
					
					Day_StartObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day_Obs_Start")
						param.hyPix.Day_Start = Int64(Day_StartObs₀[1])
					
					Hour_StartObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour_Obs_Start")
						param.hyPix.Hour_Start = Int64(Hour_StartObs₀[1])
					
					Minute_StartObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute_Obs_Start")
						param.hyPix.Minute_Start = Int64(Minute_StartObs₀[1])
					
					Second_StartObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second_Obs_Start")
						param.hyPix.Second_Start = Int64(Second_StartObs₀[1])


					Year_EndObs₀ , ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year_Obs_End")
						param.hyPix.Year_End = Int64(Year_EndObs₀[1])
					
					Month_EndObs₀, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month_Obs_End")
							param.hyPix.Month_End = Int64(Month_EndObs₀[1])
					
					Day_EndObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day_Obs_End")
						param.hyPix.Day_End = Int64(Day_EndObs₀[1])

					Hour_EndObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour_Obs_End")
						param.hyPix.Hour_End = Int64(Hour_EndObs₀[1])
					
					Minute_EndObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute_Obs_End")
						param.hyPix.Minute_End = Int64(Minute_EndObs₀[1])
					
					Second_EndObs₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second_Obs_End")
						param.hyPix.Second_End = Int64(Second_EndObs₀[1])

				# Dates of sim data
					Year_StartSim₀ , ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year_Sim_Start")
						param.hyPix.obsTheta.Year_Start = Int64(Year_StartSim₀[1])

					Month_StartSim₀, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month_Sim_Start")
							param.hyPix.obsTheta.Month_Start = Int64(Month_StartSim₀[1])

					Day_StartSim₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day_Sim_Start")
						param.hyPix.obsTheta.Day_Start = Int64(Day_StartSim₀[1])

					Hour_StartSim₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour_Sim_Start")
						param.hyPix.obsTheta.Hour_Start = Int64(Hour_StartSim₀[1])

					Minute_StartSim₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute_Sim_Start")
						param.hyPix.obsTheta.Minute_Start = Int64(Minute_StartSim₀[1])

					Second_StartSim₀, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second_Sim_Start")
						param.hyPix.obsTheta.Second_Start = Int64(Second_StartSim₀[1])

				# Dates of observed data
               param.hyPix.obsTheta.Year_End   = param.hyPix.Year_End
               param.hyPix.obsTheta.Month_End  = param.hyPix.Month_End
               param.hyPix.obsTheta.Day_End    = param.hyPix.Day_End
               param.hyPix.obsTheta.Hour_End   = param.hyPix.Hour_End
               param.hyPix.obsTheta.Minute_End = param.hyPix.Minute_End
               param.hyPix.obsTheta.Second_End = param.hyPix.Second_End


				# Dates of plots
               param.hyPix.ploting.Year_Start   = param.hyPix.obsTheta.Year_Start
               param.hyPix.ploting.Month_Start  = param.hyPix.obsTheta.Month_Start
               param.hyPix.ploting.Day_Start    = param.hyPix.obsTheta.Day_Start
               param.hyPix.ploting.Month_Start  = param.hyPix.obsTheta.Month_Start
               param.hyPix.ploting.Second_Start = param.hyPix.obsTheta.Second_Start

               param.hyPix.ploting.Year_End   = param.hyPix.Year_End
               param.hyPix.ploting.Month_End  = param.hyPix.Month_End
               param.hyPix.ploting.Day_End    = param.hyPix.Day_End
               param.hyPix.ploting.Hour_End   = param.hyPix.Hour_End
               param.hyPix.ploting.Minute_End = param.hyPix.Minute_End
               param.hyPix.ploting.Second_End = param.hyPix.Second_End

				# Checking for error
					Date_Start = DateTime(param.hyPix.Year_Start, param.hyPix.Month_Start, param.hyPix.Day_Start, param.hyPix.Hour_Start, param.hyPix.Minute_Start, param.hyPix.Second_Start)
		
					Date_SimStart = DateTime(param.hyPix.obsTheta.Year_Start, param.hyPix.obsTheta.Month_Start, param.hyPix.obsTheta.Day_Start, param.hyPix.obsTheta.Hour_Start, param.hyPix.obsTheta.Minute_Start, param.hyPix.obsTheta.Second_Start)
						
					Date_End = DateTime(param.hyPix.Year_End, param.hyPix.Month_End, param.hyPix.Day_End, param.hyPix.Hour_End, param.hyPix.Minute_End, param.hyPix.Second_End)

					if !(Date_End ≥ Date_SimStart && Date_SimStart ≥ Date_Start)
						error("HyPix: Date_Start=$Date_Start >  Date_SimStart= $Date_SimStart > Date_End=$Date_End")
					end
			return param
			end  # function: DATES

			
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : DISCRETIZATION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function DISCRETIZATION(PathDiscretization)
				# Read data
					Data = DelimitedFiles.readdlm(PathDiscretization, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]

				Z, NiZ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Z")
				Layer, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Layer")
				N_Layer = maximum(Layer)

				# Depending on the initial boundary condition 
					if "θini" ∈ Header
						θini, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "θini")
						Ψini=[]
						Flag_θΨini = :θini 

					elseif "Ψini" ∈ Header
						Ψini, ~ =  tool.readWrite.READ_HEADER_FAST(Data, Header, "Ψini")
						θini =[]
						Flag_θΨini = :Ψini

					else
						error("In $PathDiscretization cannot find <θini> or <Ψini> in $Header")
					end

			return Flag_θΨini, Layer, N_Layer, NiZ, Z, θini, Ψini
			end # function DISCRETIZATION
		#-------------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYPIX_PARAM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYPIX_PARAM(Layer, hydro, hydroHorizon, iOpt::Int64, NiZ::Int64, option, param, Path::String, veg)
				# Read data
					Data = DelimitedFiles.readdlm(Path, ',')
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
					Param, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SIM_$(iOpt)")
				
				# Minimum value of the param
					Param_Min, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MIN")

				# Maximum value of the param
					Param_Max, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MAX")

				# Determening which param to optimize
					Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT_$(iOpt)")

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
					if Flag_Opt && (iOpt ≥ param.hyPix.iOpt_Start + 1)
						Flag_MultiStepOpt = true
						println(veg)
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
							if Symbol(Name_Unique[i]) == :θsMacMat_ƞ
								for iZ =1:length(Param_Vect)
									hydroHorizon.θsMacMat[iZ] = hydroHorizon.θs[iZ] * Param_Vect[iZ]
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
					hydro = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ, option)
				end

				NparamOpt = length(ParamOpt)

				# CHECKING FOR UNCONSISTENCY WITH OPTIONS	
				if Flag_Opt && option.hyPix.σ_2_Ψm⍰ ≠ "No" && "Ψm" ∈ ParamOpt
					iψm = findfirst(isequal("Ψm"), ParamOpt)[1]

					if option.hyPix.σ_2_Ψm⍰=="UniqueRelationship" && "Ψm" ∈ ParamOpt
						error( "**** HyPix Error: combination of options which are not possible (option.hyPix.σ_2_Ψm⍰==UniqueRelationship) && (Optimise=Ψm)!")

					elseif option.hyPix.σ_2_Ψm⍰=="Constrained" && !("Ψm" ∈ ParamOpt)
						error("*** HyPix Error: combination of options which are not possible (option.hyPix.σ_2_Ψm⍰==Constrained) && (not Optimising=Ψm)!")

					elseif option.hyPix.σ_2_Ψm⍰=="Constrained" && ParamOpt_LogTransform[iψm]==1
						error("*** option.hyPix.σ_2_Ψm⍰==Constrained CANNOT log transforme Ψm") 
					end
				end # Flag_Opt

				# Putting all the parameters in  NamedTuple
				optim = (ParamOpt_Min=ParamOpt_Min, ParamOpt_Max=ParamOpt_Max, ParamOpt_HorizonEq=ParamOpt_HorizonEq, ParamOpt_Type=ParamOpt_Type, ParamOpt=ParamOpt, NparamOpt=NparamOpt, Flag_Opt=Flag_Opt, ParamOpt_LogTransform=ParamOpt_LogTransform)
				# θs_Min=θs_Min, θs_Max=θs_Max

				if Flag_Opt == true
					println("	=== === Optimizing the following parameters === ===")
					println("		NparamOpt=" , NparamOpt)
					println("		ParamOpt= " , optim.ParamOpt_Type .* optim.ParamOpt)
					println("		Min_Value= " , optim.ParamOpt_Min)
					println("		Max_Value= " , optim.ParamOpt_Max)
					println("		Hydro_HorizonEq= " , optim.ParamOpt_HorizonEq)
					println("		LogTransform = " , optim.ParamOpt_LogTransform)
					println("	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n")
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
				# if option.hyPix.ClimateDataTimestep⍰ == "Daily"
					Pr_Name          = "Rain(mm)"
					Pet_Name         = "PET(mm)"
					Temperature_Name = "Tmax(C)" # Maximum temperature which is not correct

				# elseif option.hyPix.ClimateDataTimestep⍰ == "Hourly"
				# 	Pr_Name          = "Pr_mm"
				# 	Pet_Name         = "Pet_mm"
				# 	Temperature_Name = "Temp_c"
				# end #  option.hyPix.ClimateDataTimestep⍰

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

				# REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
					Date_Start = DateTime(param.hyPix.Year_Start, param.hyPix.Month_Start, param.hyPix.Day_Start, param.hyPix.Hour_Start, param.hyPix.Minute_Start, param.hyPix.Second_Start)
					
					Date_End = DateTime(param.hyPix.Year_End, param.hyPix.Month_End, param.hyPix.Day_End, param.hyPix.Hour_End, param.hyPix.Minute_End, param.hyPix.Second_End)

				# CHECKING & CORRECTING
					# End Date feasible
						Date_End_Maximum = DateTime(Year[N_Climate], Month[N_Climate], Day[N_Climate], Hour[N_Climate], Minute[N_Climate], Second[N_Climate]) 

						Date_End = min(Date_End_Maximum, Date_End)

					# Start Date feasible
						Date_Start_Minimum = DateTime(Year[2], Month[2], Day[2], Hour[2], Minute[2], Second[2]) 

						Date_Start = max(Date_Start_Minimum , Date_Start)

				# SELECTING DATES OF INTEREST
					True = falses(N_Climate)
					Date = fill(Date_Start_Minimum::DateTime, N_Climate) 
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

					# Update N_Climate	
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
		#---------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θOBSERVATION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Base.@kwdef mutable struct θOBSERVATION
				Date    :: Vector{DateTime}
				Z  	  :: Vector{Float64}
				ithetaObs :: Vector{Int64}
				Nit    :: Int64 # Number of time steps
				Ndepth  :: Int64 # Numver of soil profile with observed θ
				θobs 	  :: Array{Float64,2}
				∑T  	  :: Vector{Float64}
			end # mutable struct

			function TIME_SERIES(clim, option, param, pathHyPix)
			# Read data
				Data = DelimitedFiles.readdlm(pathHyPix.obsTheta, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

				Year, Nit   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Year")
				Month, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "Month")
				Day, ~     = tool.readWrite.READ_HEADER_FAST(Data, Header, "Day")
				Hour, ~    = tool.readWrite.READ_HEADER_FAST(Data, Header, "Hour")
				Minute, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Minute")
				Second, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Second")

				# READING THE DEPTH OF Θ MEASUREMENTS FROM HEADER: data having Z=
					# Data, Header = DelimitedFiles.readdlm(pathHyPix.obsTheta, ','; header=true)

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
						Nit,~ = size(Data)
						θobs = Data[1:Nit, minimum(Array_iHeader): maximum(Array_iHeader)]

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

				# REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
					Date_Start_Calibr = DateTime(param.hyPix.obsTheta.Year_Start, param.hyPix.obsTheta.Month_Start, param.hyPix.obsTheta.Day_Start, param.hyPix.obsTheta.Hour_Start, param.hyPix.obsTheta.Minute_Start, param.hyPix.obsTheta.Second_Start)
					
					Date_End_Calibr = DateTime(param.hyPix.obsTheta.Year_End, param.hyPix.obsTheta.Month_End, param.hyPix.obsTheta.Day_End, param.hyPix.obsTheta.Hour_End, param.hyPix.obsTheta.Minute_End, param.hyPix.obsTheta.Second_End)

				# Compared to observed climate data
               Date_Start_Calibr = max(Date_Start_Calibr, clim.Date[2])
               Date_End_Calibr   = min(Date_End_Calibr, clim.Date[end])

				# SELECTING THE DATA WITHING FEASIBLE RANGE
					True = falses(Nit) # Initiating with false
					Date = fill(Date_Start_Calibr::DateTime,  Nit)
					iCount = 0 
					for iT=1:Nit
						Date[iT] = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])

						if (Date_Start_Calibr ≤ Date[iT] ≤ Date_End_Calibr)
							iCount += 1
							True[iT] = true
						end  # if
					end # for iT=1:Nit

					# New reduced number of simulations selected with dates
						Date = Date[True[1:Nit]]
						θobs = θobs[True[1:Nit], 1:Ndepth]

						Nit = count(True[1:Nit]) # New number of data

				# REDUCING THE AMOUNT OF DATA TO HOURLY TODO

				# ∑T_Reduced = collect(range(Date[1], step=Hour[1], stop=Date[end])) 
					# ΔTimeStep = param.hyPix.ΔT_Output
					# if option.hyPix.θobs_Hourly && ΔTimeStep < 86400
					# 	True = falses(Nit)
					# 	iCount = 0 
					# 	for iT=1:Nit
					# 		if hour(Date[iT]) == 0 && minute(Date[iT]) == 0
					# 			True[iT] = true
					# 			iCount += 1
					# 		end # if
					# 	end # for
					
					# 	# New reduced number of simulations selected with dates
					# 	Date = Date[True[1:Nit]]
					# 	θobs = θobs[True[1:Nit],1:Ndepth]
					# 	Nit = count(True[1:Nit]) # New reduced amount of data
					# end # θobs_Hourly

				# This will be computed at PrioProcess
					∑T        = fill(0.0::Float64, Nit)
					ithetaObs = fill(0::Int64, Ndepth)

				# SAVING SPACE 
					Data = nothing
					True = nothing

				# STRUCTURE
					return obsTheta = θOBSERVATION(Date, Z, ithetaObs, Nit, Ndepth, θobs, ∑T)
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
	
	
	module nsdr
	   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #		FUNCTION : θψLAB_2D_2_1D
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θψLAB_2D_2_1D(Path)
				println("    ~  $(Path) ~")

				# Read data
					Data = DelimitedFiles.readdlm(Path, ',')
				# Read header
					Header = Data[1,1:end]
				# Remove first READ_ROW_SELECT
					Data = Data[2:end,begin:end]
				# Sort data
					Data = sortslices(Data, dims=1)

				# Read data of interest
					Id₂, NiZ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")

					Soilname₂, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")

					Ψdata = []
					θData = []
					for iHeader in Header
						if occursin("wrc", iHeader)
							θ₀, NiZ = tool.readWrite.READ_HEADER_FAST(Data, Header, iHeader)

							iHeader = replace(iHeader, "wrc" => "")
							iHeader = replace(iHeader, "kpa" => "")
							iHeader = replace(iHeader, " " => "")
							iHeader_Float=  parse(Float64, iHeader)

							iHeader_Float = iHeader_Float * cst.kPa_2_Mm

							append!(Ψdata, iHeader_Float)

							try
								θData = hcat(θData[1:NiZ, :], θ₀[1:NiZ])
							catch
								θData = θ₀[1:NiZ]
							end
						end # occursin("wrc", iHeader)
					end # for iHeader in Header

					θ_θΨobs₂ = zeros(Float64, NiZ, length(Ψdata))
					Ψ_θΨobs₂ = zeros(Float64, NiZ, length(Ψdata))
					N_θΨobs₂ = zeros(Int64, NiZ)
	
					for iZ=1:NiZ
						iΨ_Count = 1
						for iΨ=1:length(Ψdata)
							if !isnan(θData[iZ, iΨ])
								Ψ_θΨobs₂[iZ, iΨ_Count] = Ψdata[iΨ]
								θ_θΨobs₂[iZ, iΨ_Count] = θData[iZ, iΨ]
								N_θΨobs₂[iZ] += 1
								iΨ_Count += 1
							end #  !isnan(θData[iZ, iΨ])
						end # iΨ
					end # iZ

			return Id₂, N_θΨobs₂, Soilname₂, θ_θΨobs₂, Ψ_θΨobs₂
		end  # function: θψLAB_2D_2_1D
		
	end
end  # module: reading
# ............................................................		