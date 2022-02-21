# =============================================================
#		module: readLinkinFile
# =============================================================
module readLinkingFile
   import DelimitedFiles: readdlm
   import ..tool.readWrite: READ_HEADER_FAST

   Base.@kwdef mutable struct PATHINPUT
      Climate::Vector{String}
		Discretization::Vector{String}
		DiscretizationAuto::Vector{String}
      θdata::Vector{String}
      Vegetation::Vector{String}
      Hydro::Vector{String}
      SoilLayer::Vector{String}
      LookUpTable_Lai::Vector{String}
      LookUpTable_Crop::Vector{String}

		# Dates::String
		# HyPix_HydroParam::String
		# HyPixParamOpt::String

		# HyPix_VegParam::String
		# IdSelect::String 
		# Input_OfStep::String
		# JulesMetadata::String
		# IdName_Hypix::String
		# obsTheta::String 
   end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : LINKING_FILES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function LINKING_FILE(Path_Hypix, SiteName_Hypix)

   
      Path_LinkingFile  = Path_Hypix * "\\data\\INPUT\\Data_Hypix\\" * SiteName_Hypix * "\\" * SiteName_Hypix * "_LinkingFile.csv"

      Path_Input = Path_Hypix * "\\data\\INPUT\\Data_Hypix\\" * SiteName_Hypix * "\\" 

      @assert isfile(Path_LinkingFile)

      Data, Header = readdlm(Path_LinkingFile, ',', header=true, use_mmap=true)

      # Selecting the data
         Id_True, ~ = READ_HEADER_FAST(Data, Header, "SELECT")
      
         IdSelect_True = convert.(Bool, Id_True)

         Data = @view Data[IdSelect_True, :]

         N_iZ_All = sum(Id_True)

      # Sorting the data with increasing Id
         Data = sortslices(Data, dims=1)

      # Reading data
         Id, ~  = READ_HEADER_FAST(Data, Header, "Id")

         B, ~ = READ_HEADER_FAST(Data, Header, "Soilname")

         a, ~ = READ_HEADER_FAST(Data, Header, "Year_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Year_Sim_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Year_End")

         a, ~ = READ_HEADER_FAST(Data, Header, "Month_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Month_Sim_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Month_End")

         a, ~ = READ_HEADER_FAST(Data, Header, "Day_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Day_Sim_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Day_End")

         a, ~ = READ_HEADER_FAST(Data, Header, "Hour_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Hour_Sim_Start")
         a, ~ = READ_HEADER_FAST(Data, Header, "Hour_End")

         a, ~ = READ_HEADER_FAST(Data, Header, "TopBoundary")
         a, ~ = READ_HEADER_FAST(Data, Header, "BottomBoundary")

         Climate = READ_PATH(Data, Header, Path_Input, "CLIMATE")
         Discretisation = READ_PATH(Data, Header, Path_Input, "DISCRETISATION")
         θdata = READ_PATH(Data, Header, Path_Input, "SOILMOISTURE")
         Vegetation = READ_PATH(Data, Header, Path_Input, "VEGETATION")
         Hydro = READ_PATH(Data, Header, Path_Input, "HYDRO")
         Optimisation = READ_PATH(Data, Header, Path_Input, "OPTIMISATION")
         SoilLayer = READ_PATH(Data, Header, Path_Input, "SOILLAYER")
         LookUpTable_Lai = READ_PATH(Data, Header, Path_Input, "LookUpTable_Lai"; PathAdd="LookUpTable\\LookUpTable_Lai")
         LookUpTable_Crop = READ_PATH(Data, Header, Path_Input, "LookUpTable_Crop"; PathAdd="LookUpTable\\LookUpTable_Crop")
         OptionHypix = READ_PATH(Data, Header, Path_Input, "OPTION"; PathAdd="ParamOptionPath\\OPTION")
         ParamHypix = READ_PATH(Data, Header, Path_Input, "PARAM"; PathAdd="ParamOptionPath\\PARAM")
  

         Option, ~ = READ_HEADER_FAST(Data, Header, "OPTION")
         Param, ~ = READ_HEADER_FAST(Data, Header, "PARAM")
         Path, ~ = READ_HEADER_FAST(Data, Header, "PATH")
         KsModel, ~ = READ_HEADER_FAST(Data, Header, "KSMODEL")

      @show Climate
      @show Discretisation
      @show θdata
      @show Vegetation

      # pathInputHypix = PATHINPUT(Climate, Dates, Discretization, DiscretizationAuto, HyPix_HydroParam, HyPixParamOpt, HyPix_VegParam, IdSelect, Input_OfStep, JulesMetadata, IdName_Hypix, obsTheta)

   return
   end  # function: LINKING_FILES
# ------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : READ_PATH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function READ_PATH(Data, Header, Path_Input, PathName; PathAdd="")
      Output, ~ = READ_HEADER_FAST(Data, Header, PathName)	
      for i = eachindex(Output)
         if PathAdd == ""
            Output[i] = Path_Input * PathName * "//" * Output[i]     
         else
            Output[i] = Path_Input * PathAdd * "//" * Output[i]
         end
         @show Output[i]
         @assert isfile(Output[i]) 
      end
      
   return Output
   end  # function: READ_PATH
# ------------------------------------------------------------------
   
end  # module: readLinkinFile
# ............................................................