# DEFINING STRUCTURE
   mutable struct EVAPOTRANSPIRATION
      Evaporation::Bool
      Transpiration::Bool
   end

   mutable struct SOIL
      Topsoil :: Bool
      Macropore :: Bool
   end

   mutable struct OPTION
      evapotranspiration::EVAPOTRANSPIRATION
      soil::SOIL 
   end

using TOML
function TOML_TEST()
   # PARSING TOML FILE
   Home2 = @__DIR__

   # perform cs..
      Home = dirname(Home2)
      Home = dirname(Home)

   # Change path name to /data/Private/
      Home = Home * "/data/" * "Private" * "/"


      # Path = PathHome *  "/Toml.toml"
      SiteName_Soilhyro = "NewFormat"

      FileDataSoilhydro_Input = Home * "INPUT/Data_SoilWater/" * SiteName_Soilhyro * "/" * SiteName_Soilhyro * "_"

      Path =FileDataSoilhydro_Input * "Path.toml"


      TomlParse = TOML.parsefile(Path)
      # println(TomlParse)

   # INITIAL VALUES OF STRUCTURE
      evapotranspiration=EVAPOTRANSPIRATION(false, false)
      evapotranspiration = TOML_2_STRUCT(evapotranspiration, TomlParse)
   
      soil = SOIL(false, false)
      soil = TOML_2_STRUCT(soil, TomlParse)

   # option = OPTION(evapotranspiration, soil)
   
   # println(option.evapotranspiration.Evaporation)
   # println(option.evapotranspiration.Transpiration)
   # println(option.soil.Topsoil)
   # println(option.soil.Macropore)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : TOML_2_STRUCT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function TOML_2_STRUCT(Structure, TomlParse)
   # LOOPING THROUGH THE DICT
   println(keys(TomlParse))

   for (iKey, iValue₀) in TomlParse
      for iValue in (keys(iValue₀))
         println(iValue)
         if iKey == string(Structure)
            setfield!(Structure, Symbol(iValue), TomlParse[iKey][iValue])
         end 
      end
   end
return Structure
end  # function: TOML_2_STRUCT

TOML_TEST()