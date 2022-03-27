# =============================================================
#		module: readHypix
# =============================================================
module readHypix

   import ..climate, ..discretisation, ..horizonLayer, ..hydroStruct, ..memory, ..optionsHypix, ..paramsHypix, ..pathsHypix, ..tableHypix, ..thetaObs, ..tool, ..vegStruct

   import Dates: value, DateTime, hour, minute, month, now, Hour
   import DelimitedFiles
   import CSV, Tables
   export READ_START, HYPIX_PARAM_OPT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : READ_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function READ_START(dateHypix, Id, iScenario, N_Scenario, Path_Hypix, pathInputHypix, ProjectHypix, SiteName)

         paramHypix = paramsHypix.PARAM_HYPIX(pathInputHypix.ParamHypix[iScenario])

         optionHypix = optionsHypix.OPTION_HYPIX(pathInputHypix.OptionHypix[iScenario])

         pathOutputHypix = pathsHypix.PATH_HYPIX(Path_Hypix, pathInputHypix.PathHypix[iScenario], ProjectHypix, SiteName[iScenario])
         
         # READING CLIMATE DATA ====
            clim = readHypix.CLIMATE(dateHypix, iScenario, pathInputHypix.Climate[iScenario])

            # Process climate
               ∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp = climate.CLIMATE(clim, optionHypix)

         # READING DISCRETISATION ===
            # Create Discretisation.csv from SoilLayer.csv
            if optionHypix.Discretisation_File_Auto⍰ == "Auto"
               # Read SoilLayer, could be either θini, Ψini
               Flag_θΨini, Layer, N_Layer, ~, Zlayer, θini_or_Ψini = readHypix.DISCRETISATION(pathInputHypix.SoilLayer[iScenario])

               # Performing auto discretisation			
               Layer, NiZ, Z, θini_or_Ψini = discretisation.DISCRETISATION_AUTO(optionHypix, paramHypix; N_Layer=N_Layer, Zlayer=Zlayer, θini_or_Ψini=θini_or_Ψini)
         
               tableHypix.DISCRETISATION_AUTO(Flag_θΨini, Layer, pathInputHypix.Discretisation[iScenario], Z, θini_or_Ψini)
            else
               # Read discretisation
               Flag_θΨini, Layer, N_Layer, NiZ, Z, θini_or_Ψini = readHypix.DISCRETISATION(pathInputHypix.Discretisation[iScenario])
            end # if optionHypix.Discretisation_File_Auto⍰ == "Auto" 

            # Process discretisation of the soil profile ~~~~~
               discret = discretisation.DISCRETISATION(NiZ, Z)

         # READING OBSERVED θ ===
            if optionHypix.θobs
               # Read observed θ
                  obsTheta = readHypix.θDATA(clim, dateHypix, iScenario, pathInputHypix.θdata[iScenario])
               # Process observed θ
                  obsTheta = thetaObs.ΘOBS(obsTheta, clim, discret, Z)
            end #  optionHypix.θobs
      
         # LOOKUP TABLE ===
            # Initialiozing vegetation parameters into veg structure
               veg = vegStruct.VEGSTRUCT() 
               Laiᵀ_η            = readHypix.LOOKUPTABLE_LAI(clim, optionHypix, pathInputHypix.LookUpTable_Lai[iScenario], veg)
               CropCoeficientᵀ_η = readHypix.LOOKUPTABLE_CROPCOEFICIENT(clim, optionHypix, pathInputHypix.LookUpTable_Crop[iScenario], veg)
               
         # INITIALIZING THE STRUCTURE ===
            # Initializing hydroHorizon structure
               hydroHorizon = hydroStruct.HYDROSTRUCT(optionHypix, N_Layer)

            # Initializing hydraulic paramHypix into structure 
               hydro = hydroStruct.HYDROSTRUCT(optionHypix, NiZ)

            # Optimisation
               hydroHorizon_best = hydroStruct.HYDROSTRUCT(optionHypix, N_Layer)
               hydro_best        = hydroStruct.HYDROSTRUCT(optionHypix, NiZ)
               veg_best          = vegStruct.VEGSTRUCT()

         # VEGETATION PARAMETERS
            if ! (optionHypix.opt.Optimisation)
               veg, ~ = tool.readWrite.READ_STRUCT_SIMPLE(veg, pathInputHypix.Vegetation[iScenario])

         # HYDRAULIC PARAMETERS
            hydroHorizon, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroHorizon, pathInputHypix.HydroInput[iScenario])

            @inbounds @simd for iZ=1:N_Layer
               hydroHorizon.So[iZ] = paramHypix.So # 1.0E-8
            end

            hydro = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ, optionHypix)

            end # optionHypix.Optimisation

         # MEMORY 
            # Multistep optimisation
               ∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, NseBest, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average = memory.MEMORY_MULTISTEPOPTIMISATION(paramHypix)

            # Memory	
            ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, Hpond, iNonConverge_iOpt, Laiᵀ, Q, Residual, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPr, ΔRunoff, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = memory.MEMORY(clim, N_∑T_Climate, NiZ, obsTheta, paramHypix)
         
   return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑∑ΔSink, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, ∑ΔQ_Bot, CccBest, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Efficiency, Flag_θΨini, Flag_θΨini, Global_WaterBalance, Global_WaterBalance_NormPr, Hpond, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iNonConverge_iOpt, Laiᵀ, Laiᵀ_η, Layer, N_∑T_Climate, N_Layer, NiZ, NseBest, obsTheta, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, Q, Residual, SwcRoots, Temp, veg, veg_best, WilmotBest, WofBest, Z, Zlayer, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPr, ΔRunoff, ΔRunTimeHypix, ΔSink, ΔT, ΔT_Average, θ, θini_or_Ψini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest
   end  # function: READ_START
   # ------------------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : DISCRETISATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function DISCRETISATION(Path::String)
         # Read data
            Data = CSV.File(Path, header=true)
            Header = string.(Tables.columnnames(Data))

            Z       = convert(Vector{Float64}, Tables.getcolumn(Data, :Z))
            Layer   = convert(Vector{Int64}, Tables.getcolumn(Data, :Layer))
            N_Layer = Int64(maximum(Layer))

            NiZ = length(Z)

         # Depending on the initial boundary condition 
            if "θini" ∈ Header
               θini_or_Ψini = convert(Vector{Float64}, Tables.getcolumn(Data, :θini))
               Flag_θΨini = :θini 

            elseif "Ψini" ∈ Header
               θini_or_Ψini = convert(Vector{Float64}, Tables.getcolumn(Data, :Ψini))
               Flag_θΨini = :Ψini

            else
               error("In $Path cannot find <θini> or <Ψini> in $Header")
            end

            θini_or_Ψini 

      return Flag_θΨini, Layer, N_Layer, NiZ, Z, θini_or_Ψini
      end # function DISCRETISATION
   #-------------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYPIX_PARAM_OPT
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function HYPIX_PARAM_OPT(Layer, hydro, hydroHorizon, iMultistep::Int64, NiZ::Int64, optionHypix, paramHypix, Path::String, veg)
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

            # N_NameUnique = length(Name_Unique)
      
         # Reading the values of the parameters for the simulation of interest
            Param, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SIM_$(iMultistep)")
         
         # Minimum value of the param
            Param_Min, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MIN")

         # Maximum value of the param
            Param_Max, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MAX")

         # Determening which param to optimize
            Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT_$(iMultistep)")

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
            if Flag_Opt && (iMultistep ≥ paramHypix.opt.iOptMultiStep_Start + 1)
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
            hydro = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, NiZ, optionHypix)
         end

         NparamOpt = length(ParamOpt)

         # CHECKING FOR UNCONSISTENCY WITH OPTIONS	
         if Flag_Opt && optionHypix.opt.σ_2_Ψm⍰ ≠ "No" && "Ψm" ∈ ParamOpt
            iψm = findfirst(isequal("Ψm"), ParamOpt)[1]

            if optionHypix.opt.σ_2_Ψm⍰=="UniqueRelationship" && "Ψm" ∈ ParamOpt
               error( "**** HyPix Error: combination of options which are not possible (optionHypix.opt.σ_2_Ψm⍰==UniqueRelationship) && (Optimise=Ψm)!")

            elseif optionHypix.opt.σ_2_Ψm⍰=="Constrained" && !("Ψm" ∈ ParamOpt)
               error("*** HyPix Error: combination of options which are not possible (optionHypix.opt.σ_2_Ψm⍰==Constrained) && (not Optimising=Ψm)!")

            elseif optionHypix.opt.σ_2_Ψm⍰=="Constrained" && ParamOpt_LogTransform[iψm]==1
               error("*** optionHypix.opt.σ_2_Ψm⍰==Constrained CANNOT log transforme Ψm") 
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
   end  # function: HYPIX_PARAM_OPT



   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : CLIMATE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      struct CLIMATEDATA
         Date      ::  Vector{DateTime}
         Pr        ::  Vector{Float64}
         Pet       ::  Vector{Float64}
         Temp      ::  Vector{Float64}
         N_Climate ::  Int64
         Pr_Through :: Vector{Float64}
      end

      function CLIMATE(dateHypix, iScenario::Int64, PathClimate::String; Option_ReadTemperature=false)

         # READ DATA
            Data = CSV.File(PathClimate, header=true)

            Year   = convert(Vector{Int64}, Tables.getcolumn(Data, :Year))
            Month  = convert(Vector{Int64}, Tables.getcolumn(Data, :Month))
            Day    = convert(Vector{Int64}, Tables.getcolumn(Data,  :Day))
            Hour   = convert(Vector{Int64}, Tables.getcolumn(Data,  :Hour))
            Minute = convert(Vector{Int64}, Tables.getcolumn(Data, :Minute))
            Second = convert(Vector{Int64}, Tables.getcolumn(Data,  :Second))
            Pr     = convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("Rain(mm)")))
            Pet    = convert(Vector{Float64}, Tables.getcolumn(Data, Symbol("PET(mm)")))
            
            N_Climate = length(Year)

            if Option_ReadTemperature 
               Temp=  convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("Tmax(C)")))
            else
               Temp = fill(24.0::Float64, N_Climate)
            end


         # REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
            Date_Start = DateTime(dateHypix.Year_Start[iScenario], dateHypix.Month_Start[iScenario], dateHypix.Day_Start[iScenario], dateHypix.Hour_Start[iScenario], dateHypix.Minute_Start[iScenario], dateHypix.Second_Start[iScenario])
            
            Date_End = DateTime(dateHypix.Year_End[iScenario], dateHypix.Month_End[iScenario], dateHypix.Day_End[iScenario], dateHypix.Hour_End[iScenario], dateHypix.Minute_End[iScenario], dateHypix.Second_End[iScenario])

         # CHECKING & CORRECTING
            # End Date feasible
               Date_End_Maximum = DateTime(Year[end], Month[end], Day[end], Hour[end], Minute[end], Second[end]) 

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
               iTrue_First = findfirst(True[:])
               True[iTrue_First-1] = true

            # New reduced number of simulations
               Date = Date[True[:]]
               Pr   = Pr[True[:]]
               Pet  = Pet[True[:]]
               Temp = Temp[True[:]]

            # Update N_Climate	
               N_Climate = count(True[:]) # New number of data
         
         # To be used after interception model
            Pr_Through = fill(0.0::Float64, N_Climate)

      # STRUCTURE
         return clim = CLIMATEDATA(Date, Pr, Pet, Temp, N_Climate, Pr_Through)

      end # function: CLIMATE
   #---------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θOBSERVATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Base.@kwdef mutable struct θOBSERVATION
         Date      :: Vector{DateTime}
         Z  	    :: Vector{Float64}
         ithetaObs :: Vector{Int64}
         Nit       :: Int64 # Number of time steps
         Ndepth    :: Int64 # Numver of soil profile with observed θ
         θobs 	    :: Array{Float64,2}
         ∑T  	    :: Vector{Float64}
      end # mutable struct

      function θDATA(clim, dateHypix, iScenario, Pathθobs)
      # Read data
         Data₀ = CSV.File(Pathθobs, header=true)

         Header = string.(Tables.columnnames(Data₀))

         Year   = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Year))
         Month  = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Month))
         Day    = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Day))
         Hour   = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Hour))
         Minute = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Minute))
         Second = convert(Vector{Int64}, Tables.getcolumn(Data₀,:Second))

         Data = Tables.matrix(Data₀)

         Nit    = length(Year)

         # READING THE DEPTH OF Θ MEASUREMENTS FROM HEADER: data having Z=
               Array_iHeader = Int64[]
               Ndepth = 0::Int64
               iCount = 0::Int64
               for iHeader in Header
                  iCount += 1
                  if occursin("Z=", iHeader) # Searching for 'Z=' in the header
                     Ndepth += 1
                     append!(Array_iHeader, iCount) 
                  end # occursin
               end # iHeader

            # Isolating data with Z= measurements
               # Nit,~ = size(Data)
               θobs = Data[:, minimum(Array_iHeader): maximum(Array_iHeader)]

            # The depths were we have θ measurements
               Z = fill(0.0::Float64, Ndepth)

               i = 0::Int64
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
            Date_Start_Calibr = DateTime(dateHypix.Year_Start_Sim[iScenario], dateHypix.Month_Start_Sim[iScenario], dateHypix.Day_Start_Sim[iScenario], dateHypix.Hour_Start_Sim[iScenario], dateHypix.Minute_Start_Sim[iScenario], dateHypix.Second_Start_Sim[iScenario])
            
            Date_End_Calibr = DateTime(dateHypix.Year_End[iScenario], dateHypix.Month_End[iScenario], dateHypix.Day_End[iScenario], dateHypix.Hour_End[iScenario], dateHypix.Minute_End[iScenario], dateHypix.Second_End[iScenario])

         # Compared to observed climate data
            Date_Start_Calibr = max(Date_Start_Calibr, clim.Date[2])
            Date_End_Calibr   = min(Date_End_Calibr, clim.Date[end])

         # SELECTING THE DATA WITHING FEASIBLE RANGE
            True = falses(Nit) # Initiating with false
            Date = fill(Date_Start_Calibr::DateTime, Nit)
            iCount = 0 ::Int64
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
            # if optionHypix.θobs_Reduced && ΔTimeStep < 86400
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
            # end # θobs_Reduced

         # This will be computed at PrioProcess
            ∑T        = fill(0.0::Float64, Nit)
            ithetaObs = fill(0::Int64, Ndepth)

         # SAVING SPACE 
            # Data = nothing
            # True = nothing

         # STRUCTURE
      return θOBSERVATION(Date, Z, ithetaObs, Nit, Ndepth, θobs, ∑T)
      end  # function: TIME_SERIES


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : LOOKUPTABLE
   #		Parameters as a function of time
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function LOOKUPTABLE_LAI(clim, optionHypix, PathLai::String, veg)	
         if optionHypix.LookupTable_Lai == true
            LookUpTable_Lai, ~   = tool.readWrite.READ_HEADER(PathLai, "Lai")
         end
         
         Laiᵀ_Norm = fill(0.0::Float64, clim.N_Climate) 
         for (i, Date) in enumerate(clim.Date)
            Month = month(Date)
            if optionHypix.LookupTable_Lai == true
               Laiᵀ_Norm[i] = LookUpTable_Lai[Month]
            else
               Laiᵀ_Norm[i] = veg.Lai
            end
         end
      
      return Laiᵀ_Norm
      end  # function: LOOKUPTABLE_LAI

      function LOOKUPTABLE_CROPCOEFICIENT(clim, optionHypix, PathLai::String, veg)
         if optionHypix.LookUpTable_CropCoeficient == true
            LookUpTable_CropCoeficient, ~   = tool.readWrite.READ_HEADER(PathLai, "CropCoeficient")
         end
         
         CropCoeficientᵀ_Norm = fill(0.0::Float64, clim.N_Climate) 
         for (i, Date) in enumerate(clim.Date)
            Month = month(Date)
            if optionHypix.LookUpTable_CropCoeficient == true
               CropCoeficientᵀ_Norm[i] = LookUpTable_CropCoeficient[Month]
            else
               CropCoeficientᵀ_Norm[i] = veg.CropCoeficient
            end
         end
      return CropCoeficientᵀ_Norm
      end  # function: LOOKUPTABLE_LAI
   
end  # module: readHypix
# ............................................................