# =============================================================
#		module: readLinkinFile
# =============================================================
module readLinkingFile
   import DelimitedFiles: readdlm
   import ..tool.readWrite: READ_HEADER_FAST

   Base.@kwdef struct PATHINPUT
      Climate          ::Vector{String}
      Discretisation   ::Vector{String}
      HydroInput       ::Vector{String}
      HydroRange       ::Vector{String}
      KsModel          ::Vector{String}
      LookUpTable_Crop ::Vector{String}
      LookUpTable_Lai  ::Vector{String}
      MultistepOpt     ::Vector{String}
      OptionHypix      ::Vector{String}
      ParamHypix       ::Vector{String}
      PathHypix        ::Vector{String}
      SoilLayer        ::Vector{String}
      Vegetation       ::Vector{String}
      θdata            ::Vector{String}
   end

   Base.@kwdef struct DATES
      Year_Start      ::Vector{UInt64}
      Year_Start_Sim  ::Vector{UInt64}
      Year_End        ::Vector{UInt64}

      Month_Start     ::Vector{UInt64}
      Month_Start_Sim ::Vector{UInt64}
      Month_End       ::Vector{UInt64}

      Day_Start       ::Vector{UInt64}
      Day_Start_Sim   ::Vector{UInt64}
      Day_End         ::Vector{UInt64}
      
      Hour_Start      ::Vector{UInt64}
      Hour_Start_Sim  ::Vector{UInt64}
      Hour_End        ::Vector{UInt64}

      Minute_Start      ::Vector{UInt64}
      Minute_Start_Sim  ::Vector{UInt64}
      Minute_End        ::Vector{UInt64}

      Second_Start      ::Vector{UInt64}
      Second_Start_Sim  ::Vector{UInt64}
      Second_End        ::Vector{UInt64}
   end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : LINKING_FILES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function LINKING_FILE(Path_Hypix::String, ProjectHypix::String)

      # PATHS
         Path_LinkingFile = Path_Hypix * "\\data\\INPUT\\Data_Hypix\\" * ProjectHypix * "\\" * ProjectHypix * "_LinkingFile.csv"

         @assert isfile(Path_LinkingFile)
   
         Path_Input = Path_Hypix * "\\data\\INPUT\\Data_Hypix\\" * ProjectHypix * "\\"
   
      # READING DATA
         Data, Header = readdlm(Path_LinkingFile, ',', header = true, use_mmap = true)
   
      # SELECTING DATA
         Id_True, ~ = READ_HEADER_FAST(Data, Header, "SELECT")
   
         IdSelect_True = convert.(Bool, Id_True)
   
         Data = @view Data[IdSelect_True, :]
   
      # NUMBER OF SCENARIOS
         N_Scenario = sum(Id_True)
   
      # Sorting the data with increasing Id
         Data = sortslices(Data, dims = 1)
   
      # READING
         Id, ~ = READ_HEADER_FAST(Data, Header, "Id")
   
         SiteName, ~ = READ_HEADER_FAST(Data, Header, "SiteName")
   
      ## READ DATES ===
         Year_Start, ~      = READ_HEADER_FAST(Data, Header, "Year_Start")
         Year_Start_Sim, ~  = READ_HEADER_FAST(Data, Header, "Year_Start_Sim")
         Year_End, ~        = READ_HEADER_FAST(Data, Header, "Year_End")
      
         Month_Start, ~     = READ_HEADER_FAST(Data, Header, "Month_Start")
         Month_Start_Sim, ~ = READ_HEADER_FAST(Data, Header, "Month_Start_Sim")
         Month_End, ~       = READ_HEADER_FAST(Data, Header, "Month_End")
      
         Day_Start, ~       = READ_HEADER_FAST(Data, Header, "Day_Start")
         Day_Start_Sim, ~   = READ_HEADER_FAST(Data, Header, "Day_Start_Sim")
         Day_End, ~         = READ_HEADER_FAST(Data, Header, "Day_End")
      
         Hour_Start, ~      = READ_HEADER_FAST(Data, Header, "Hour_Start")
         Hour_Start_Sim, ~  = READ_HEADER_FAST(Data, Header, "Hour_Start_Sim")
         Hour_End, ~        = READ_HEADER_FAST(Data, Header, "Hour_End")

         Minute_Start     = fill(0.0::Float64, N_Scenario)
         Minute_Start_Sim = fill(0.0::Float64, N_Scenario)
         Minute_End       = fill(0.0::Float64, N_Scenario)
   
         Second_Start     = fill(0.0::Float64, N_Scenario)
         Second_Start_Sim = fill(0.0::Float64, N_Scenario)
         Second_End       = fill(0.0::Float64, N_Scenario)
   
         dateHypix = DATES(Year_Start, Year_Start_Sim, Year_End, Month_Start, Month_Start_Sim, Month_End,Day_Start, Day_Start_Sim, Day_End, Hour_Start, Hour_Start_Sim, Hour_End, Minute_Start, Minute_Start_Sim, Minute_End, Second_Start, Second_Start_Sim , Second_End)
   
      ## READING PATH ===
         Climate = READ_PATH(Data, Header, Path_Input, "CLIMATE")
         Discretisation = READ_PATH(Data, Header, Path_Input, "DISCRETISATION")
         HydroInput = READ_PATH(Data, Header, Path_Input, "HYDRO_INPUT", AllowMissing = true)
         HydroRange = READ_PATH(Data, Header, Path_Input, "HYDRO_RANGE"; PathAdd = "ParamOptionPath\\HYDRO_RANGE", AllowMissing = true)
         KsModel = READ_PATH(Data, Header, Path_Input, "KSMODEL"; PathAdd = "ParamOptionPath\\KSMODEL")
         LookUpTable_Crop = READ_PATH(Data, Header, Path_Input, "LookUpTable_Crop"; PathAdd = "LookUpTable\\LookUpTable_Crop")
         LookUpTable_Lai = READ_PATH(Data, Header, Path_Input, "LookUpTable_Lai"; PathAdd = "LookUpTable\\LookUpTable_Lai")
         MultistepOpt = READ_PATH(Data, Header, Path_Input, "MULTISTEP_OPT", AllowMissing = true)
         OptionHypix = READ_PATH(Data, Header, Path_Input, "OPTION"; PathAdd = "ParamOptionPath\\OPTION")
         ParamHypix = READ_PATH(Data, Header, Path_Input, "PARAM"; PathAdd = "ParamOptionPath\\PARAM")
         PathHypix = READ_PATH(Data, Header, Path_Input, "PATH"; PathAdd = "ParamOptionPath\\PATH")
         SoilLayer = READ_PATH(Data, Header, Path_Input, "SOILLAYER")
         Vegetation = READ_PATH(Data, Header, Path_Input, "VEGETATION")
         θdata = READ_PATH(Data, Header, Path_Input, "SOILMOISTURE")
      
         pathInputHypix = PATHINPUT(Climate, Discretisation, HydroInput, HydroRange, KsModel, LookUpTable_Crop, LookUpTable_Lai, MultistepOpt, OptionHypix, ParamHypix, PathHypix, SoilLayer, Vegetation, θdata)
   
   return dateHypix, Id, N_Scenario, pathInputHypix, SiteName
   end  # function: LINKING_FILES
# ------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : READ_PATH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function READ_PATH(Data::Matrix{Any}, Header::Matrix{AbstractString}, Path_Input::String, PathName::String; PathAdd=""::String, AllowMissing=false::Bool)

      Output, ~ = READ_HEADER_FAST(Data, Header, PathName)	
      for i = eachindex(Output)
         if Output[i] ≠  "missing"
            if PathAdd == ""
               Output[i] = Path_Input * PathName * "//" * Output[i]     
            else
               Output[i] = Path_Input * PathAdd * "//" * Output[i]
            end
            if PathName ≠ "DISCRETISATION"
               @assert isfile(Output[i])
            end
         else
            if AllowMissing==false
               error("$PathName not allowed missing data")
            end
         end
      end # for
      
   return Output
   end  # function: READ_PATH
# ------------------------------------------------------------------
   
end  # module: readLinkinFile
# ............................................................