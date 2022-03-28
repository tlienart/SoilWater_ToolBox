# =============================================================
#		MODULE: reading
# =============================================================
module reading
	import ..tool, ..table
	import  DelimitedFiles
	export ID, θΨ, KUNSATΨ, INFILTRATION, PSD

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

				Id, N_Scenario  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")
			
				Id_True, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, PathOptionSelect)

				Id = Int64.(Id)

			# Soilname is optional
				Soilname = []
				try
					Soilname, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")
				catch # If not available
					Soilname = fill("", N_Scenario)
					for i=1:N_Scenario
						Soilname[i] = PathModelName  * "_" * string(Id[i])
					end
				end
		
				Id_True = Int64.(Id_True)

				IdSelect_True = convert.(Bool, Id_True)

			# Checking for errors
				for iZ=2:N_Scenario
					if (Id[iZ] - Id[iZ-1]) < 1
						error("Id does not increase monotically at iD $(Id[iZ]) ")
					end
				end # for iZ=2:N_Scenario
		
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