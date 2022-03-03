using Suppressor

Path_SoilWater = dirname(dirname(@__DIR__)) * "/src/" # moving down the path twice

@suppress begin	
   include(Path_SoilWater * "Cst.jl")
   include(Path_SoilWater * "Hydro/Wrc.jl")
   include(Path_SoilWater * "Hydro/Kunsat.jl")
   include(Path_SoilWater * "Tool.jl")
   include(Path_SoilWater * "Hydro/HydroStruct.jl")
   include(Path_SoilWater * "Hydro/HydroRelation.jl")
   include(Path_SoilWater * "Stats.jl")
   include("HorizonLayer.jl")

   include("ReadWriteHypix/ReadLinkingFile.jl")
   include("ReadWriteHypix/OptionsHypix.jl")
   include("ReadWriteHypix/ParamsHypix.jl")
   include("ReadWriteHypix/PathsHypix.jl")
   include("Opt/ThetaObs.jl")
   include("Climate.jl")
   include("Discretisation.jl")
   include("Memory.jl")
   include("Veg/VegStruct.jl")
   include("OutputHypix/TableHypix.jl")

   include("ReadWriteHypix/ReadHypix.jl")
   include("WaterBalance.jl")
   include("Veg/Interception.jl")
   include("Veg/RootWaterUptake.jl")
   include(Path_SoilWater * "Hydro/ΨminΨmax.jl")
   include("TimeStep.jl")
   include("Evaporation.jl")
   include("Interpolate.jl")
   include("Opt/OfHypix.jl")
   include("Veg/Pet.jl")
   include(Path_SoilWater * "Sorptivity/Sorptivity.jl")
   include("Ponding.jl")
   include("Flux.jl")
   include("Residual.jl")
   include("Richard.jl")
   include("HypixModel.jl")
   include("Opt/HypixOpt.jl")
   include("Other/θaver.jl")
   include("ΔΔtchange.jl")

   # include("θini.jl")

   include("OutputHypix/PlotHypix.jl")
   # include("OutputHypix/PlotOther.jl")

end # @suppress begin	
