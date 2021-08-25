Base.@kwdef mutable struct KSMODELτ
   τ₁        :: Vector{Float64}
   τ₂        :: Vector{Float64}
   τ₃        :: Vector{Float64}
   τ₁Mac     :: Vector{Float64}
   τ₂Mac     :: Vector{Float64}
   τ₃Mac     :: Vector{Float64}

   τ₁_Min    :: Vector{Float64}
   τ₂_Min    :: Vector{Float64}
   τ₃_Min    :: Vector{Float64}
   τ₁Mac_Min :: Vector{Float64}
   τ₂Mac_Min :: Vector{Float64}
   τ₃Mac_Min :: Vector{Float64}
   
	τ₁_Max    :: Vector{Float64}
   τ₂_Max    :: Vector{Float64}
   τ₃_Max    :: Vector{Float64}
   τ₁Mac_Max :: Vector{Float64}
   τ₂Mac_Max :: Vector{Float64}
   τ₃Mac_Max :: Vector{Float64}
 end # mutable struct KSMODEL


 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 #		FUNCTION : STRUCT_KSMODEL
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function STRUCT_KSMODEL(;Nτ_Layer = 2)
		τ₁       = fill(0.0::Float64, Nτ_Layer)
		τ₂       = fill(0.0::Float64, Nτ_Layer)
		τ₃       = fill(0.0::Float64, Nτ_Layer)
		τ₁Mac    = fill(0.0::Float64, Nτ_Layer)
		τ₂Mac    = fill(0.0::Float64, Nτ_Layer)
		τ₃Mac    = fill(0.0::Float64, Nτ_Layer)

		τ₁_Min   = fill(0.0::Float64, Nτ_Layer)
		τ₂_Min   = fill(0.0::Float64, Nτ_Layer)
		τ₃_Min   = fill(0.0::Float64, Nτ_Layer)
		τ₁Mac_Min= fill(0.0::Float64, Nτ_Layer)
		τ₂Mac_Min= fill(0.0::Float64, Nτ_Layer)
		τ₃Mac_Min= fill(0.0::Float64, Nτ_Layer)
		
		τ₁_Max   = fill(0.0::Float64, Nτ_Layer)
		τ₂_Max   = fill(0.0::Float64, Nτ_Layer)
		τ₃_Max   = fill(0.0::Float64, Nτ_Layer)
		τ₁Mac_Max= fill(0.0::Float64, Nτ_Layer)
		τ₂Mac_Max= fill(0.0::Float64, Nτ_Layer)
		τ₃Mac_Max= fill(0.0::Float64, Nτ_Layer)

		ksmodelτ = KSMODELτ(
			τ₁,
			τ₂,
			τ₃,
			τ₁Mac,
			τ₂Mac,
			τ₃Mac,
			τ₁_Min,
			τ₂_Min,
			τ₃_Min,
			τ₁Mac_Min,
			τ₂Mac_Min,
			τ₃Mac_Min,
			τ₁_Max,
			τ₂_Max,
			τ₃_Max,
			τ₁Mac_Max,
			τ₂Mac_Max,
			τ₃Mac_Max)

	return ksmodelτ 
	end  # function: STRUCT_KSMODEL


	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  KSMODEL_PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Base.@kwdef mutable struct OPTIMKS
		Param_Name :: Vector{String}
		ParamOpt_Min :: Vector{Float64}
		ParamOpt_Max :: Vector{Float64}
		Param_Min :: Vector{Float64}
		Param_Max :: Vector{Float64}
		ParamOpt :: Vector{String}
		NparamOpt :: Int64
		Flag_Opt :: Bool
	end

	function KSMODEL_PARAM(optionₘ, ksmodelτ, Path)
	# Read data
		Data = DelimitedFiles.readdlm(Path, ',')
	# Read header
		Header = Data[1,1:end]
	# Remove first READ_ROW_SELECT
		Data = Data[2:end,begin:end]

	# Reading the Model data
		KsModel⍰, Ndata   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MODEL")

	# Determening which parameters correspond to the selected model
	iSelectModel = [] 
	for i=1:Ndata
		if KsModel⍰[i] == string(option.ksModel.KsModel⍰)
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

	Nlayer = 2
	i = 1
	# For every hydraulic parameter
	for inParamValue in Param_Name
		# Putting the value of the parameters in ksmodelτ. Repeating the value of the parameter for all soils data: Nlayer

		if occursin("_Top", inParamValue)
			inParamValue = replace(inParamValue, "_Top" => "" )
			iLayer = 1

		elseif occursin("_Sub", inParamValue)
			inParamValue = replace(inParamValue, "_Sub" => "" )
			iLayer = 2
		end

		ParamValue_Vector = getfield(ksmodelτ, Symbol(inParamValue))
		ParamValue_Vector[iLayer] = inParamValue

		# ParamValue_Vector = fill(Float64(ParamValue[i]), Nlayer)
		setfield!(ksmodelτ, Symbol(inParamValue), ParamValue_Vector)

		# Putting the minimum value in the parameter
			ParamValue_Vector = fill(Float64(Param_Min[i]), Nlayer)
			setfield!(ksmodelτ, Symbol(inParamValue * "_Min"), ParamValue_Vector)

		# Putting the maximum value in the parameter
			ParamValue_Vector = fill(Float64(Param_Max[i]), Nlayer)
			setfield!(ksmodelτ, Symbol(inParamValue * "_Max"), ParamValue_Vector)

		# ParamValue to optimize. The complication is that there are different layers of hydraulic parameters which can be optimized.  
		if Opt[i] == 1
			# appending the values of the parameters
			append!(ParamOpt, [Param_Name[i]])

			append!(ParamOpt_Min, Param_Min[i])
			
			append!(ParamOpt_Max, Param_Max[i])

			if Param_Min[i] > Param_Max[i]
				error("LabOpt ERROR: $(Param_Min[i]) < $(ParamValue[i]) < $(Param_Max[i]) !")
			end
		end # if Flag_Opt

		i += 1
	end # for loop

	# Number of parameters to be optimised
		NparamOpt = length(ParamOpt)

	# Putting all the in mutable structure
		optim = 	OPTIMKS(Param_Name, ParamOpt_Min, ParamOpt_Max, Param_Min, Param_Max, ParamOpt, NparamOpt, Flag_Opt)

	if Flag_Opt == true
		println("	=== === Optimizing the following parameters === ===")
		println("		Model=" , option.ksModel.KsModel⍰)
		println("		NparamOpt=" , NparamOpt)
		println("		ParamOpt= " ,  optim.ParamOpt)
		println("		Min_Value= " , optim.ParamOpt_Min)
		println("		Max_Value= " , optim.ParamOpt_Max)
		println("	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n")
	end

return ksmodelτ, optim
end  # function: KSMODEL_PARAM
