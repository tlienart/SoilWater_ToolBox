# DEFINING STRUCTURE
using Base
   Base.@kwdef mutable struct EVAPOTRANSPIRATION
      Evaporation
      Transpiration
   end

   Base.@kwdef mutable struct SOIL
      Topsoil :: Bool
      Macropore :: Bool
   end

   Base.@kwdef mutable struct OPTION
      evapotranspiration::EVAPOTRANSPIRATION
      soil::SOIL 
   end

using TOML
function TOML_TEST()

   # PARSING TOML FILE
      Path= "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\Temporary\\Toml.toml"
      TomlParse = TOML.parsefile(Path)

   # STRUCTURE

   println(keys(TomlParse)) 
      evapotranspiration = TOML_2_STRUCT(EVAPOTRANSPIRATION, TomlParse; MyType_LowerCase=true)


      soil = TOML_2_STRUCT(SOIL,TomlParse; MyType_LowerCase=true)

   option = OPTION(evapotranspiration, soil)
   
   # TESTING
      println(option.evapotranspiration.Evaporation)
      println(option.evapotranspiration.Transpiration)
      println(option.soil.Topsoil)
      println(option.soil.Macropore)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : TOML_2_STRUCT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function TOML_2_STRUCT2(Structure, TomlParse)
      # LOOPING THROUGH THE DICT
      for (iKey, iValue₀) in TomlParse
      for iValue in (keys(iValue₀))
         if uppercase.(iKey) == (string(typeof(Structure)))
            setfield!(Structure, Symbol(iValue), TomlParse[iKey][iValue])
         end 
      end
   end
   return Structure
   end  # function: TOML_2_STRUCT

   function TOML_2_STRUCT(Structure, TomlParse; MyType_LowerCase=true, MyType=:MyType)
      if MyType_LowerCase == false
         MyType = string(MyType)
      else
         MyType = lowercase.(string(Structure))
      end

      Output = NamedTuple{Tuple(Symbol.(keys(TomlParse[MyType])))}(values(TomlParse[MyType]))
     return Structure(Output...)
   end



TOML_TEST()