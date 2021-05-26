mutable struct OPTION
   transpiration
   soil 
end

   mutable struct TRANSPIRATION
      Evaporation :: Bool
      Transpiration :: Bool
   end

   mutable struct SOIL
      Topsoil :: Bool
      Macropore :: Bool
   end

Evaporation = true
Transpiration = false
   transpiration = TRANSPIRATION(Evaporation, Transpiration)

   Topsoil = true
   Macropore = false
soil = SOIL(Topsoil, Macropore)

option = OPTION(transpiration, soil)