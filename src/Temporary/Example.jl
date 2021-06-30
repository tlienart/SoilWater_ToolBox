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

using TOML, StructTypes
function TOML_TEST()
   # PARSING TOML FILE
      PathHome = @__DIR__
      Path = PathHome *  "/Toml.toml"
      # Dict{String, Any}("evapotranspiration" => Dict{String, Any}("Evaporation" => true, "Transpiration" => false), "soil" => Dict{String, Any}("Topsoil" => true, "Macropore" => false))
      TomlParse = TOML.tryparsefile(Path)

   # INITIAL VALUES OF STRUCTURE
      Evaporation::Bool =false
      Transpiration::Bool=false
      evapotranspiration=EVAPOTRANSPIRATION(Evaporation, Transpiration)
   
      Topsoil::Bool=false
      Macropore::Bool=false
      soil = SOIL(Topsoil, Macropore)

   # LOOPING THROUGH THE DICT
   for (iKey, iValue₀) in TomlParse

      for iValue in (keys(iValue₀))
         if iKey == "evapotranspiration"
            setfield!(evapotranspiration, Symbol(iValue), TomlParse[iKey][iValue])

         elseif iKey == "soil"
            setfield!(soil, Symbol(iValue), TomlParse[iKey][iValue])
         end 
      end
   end

   option = OPTION(evapotranspiration, soil)
   
   println(option.evapotranspiration.Evaporation)
   println(option.evapotranspiration.Transpiration)
   println(option.soil.Topsoil)
   println(option.soil.Macropore)
end

TOML_TEST()