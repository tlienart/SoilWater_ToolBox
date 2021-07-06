# DEFINING STRUCTURE
using Configurations, YAML
   @option  mutable struct EVAPOTRANSPIRATION
      α::Symbol
      Evaporation::Int64
   end

   @option  mutable struct SOIL
      Topsoil::Float64 
      Macropore::String 
   end

   @option  mutable struct OPTION
      evapotranspiration::EVAPOTRANSPIRATION
      soil::SOIL 
   end

function TOML_TEST()
   # PARSING TOML FILE
      Path= "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\Temporary\\Toml.yaml"
      # TomlParse = TOML.parsefile(Path)

      YamlDict = YAML.load_file(Path)

      println(YamlDict)

      option = Configurations.from_dict(OPTION, YamlDict)
      # option = Configurations.from_toml(OPTION, Path)
  
   # TESTING
      println(option.evapotranspiration.Evaporation)
      println(option.evapotranspiration.α)
      println(option.soil.Macropore)
      println(option.soil.Topsoil)
end





TOML_TEST()