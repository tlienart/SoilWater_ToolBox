# =============================================================
#		module: including
# =============================================================
# module including

using Suppressor

@suppress begin
    include("Option.jl")
    include("Path.jl")

    # include("Packages.jl")
    include("Tool.jl")
    include("Cst.jl")
    include("Param.jl")
    include("Hydro/ΨminΨmax.jl")
    include("Hydro/HydroStruct.jl")
    include("Hypix/HorizonLayer.jl")
    include("Hydro/HydroRelation.jl")
    include("Hydro/Wrc.jl")
    include("Optim/Optimize.jl")
    include("reading.jl")
    include("Hydro/Kunsat.jl")
    include("Table.jl")
    include("hydro/θψ2Ks.jl")
    include("Checking.jl")
    include("Stats.jl")
    include("RockFragment/RockFragment.jl")
    
    include("Psd/PsdThetar.jl")

    include("Hypix/VegStruct.jl")
    include("Hypix/Discretization.jl")
    include("NoCore/Smap/ReadSmap.jl")
    include("NoCore/Smap/TableSmap.jl")
    include("NoCore/Smap/Smap2Hypix.jl")

    include("NoCore/Smap/PlotSmap.jl")

    include("HydroLab/OfHydrolab.jl")
    include("HydroLab/HydrolabOpt.jl")

    include("Sorptivity/Sorptivity.jl")            
    include("Infilt/BestFunc.jl")
    include("Infilt/OfBest.jl")
    include("Infilt/QuasiExact.jl")
    include("Infilt/TimeTransSteady.jl")
    include("Infilt/InfiltStruct.jl")
    include("Infilt/InfiltInitialize.jl")
    include("Infilt/Infilt_START.jl")

    include("Psd/PsdStruct.jl")
    include("Psd/PsdInitialize.jl")
    include("Psd/PsdFunc.jl")
    include("Psd/PsdOpt.jl")
    include("Psd/Psd_START.jl")

    include("Plot.jl")

    include("Hypix/Interpolate.jl")
    include("Hypix/Opt/ThetaObs.jl")
    include("HyPix/θini.jl")
    include("Hypix/Opt/OfHypix.jl")
    include("Hypix/Interception.jl")
    include("Hypix/Flux.jl")
    include("Hypix/Ponding.jl")
    include("Hypix/Residual.jl")
    include("Hypix/Δtchange.jl")
    include("Hypix/TimeStep.jl")
    include("Hypix/Richard.jl")
    include("Hypix/WaterBalance.jl")
    include("Hypix/Evaporation.jl")
    include("Hypix/RootWaterUptake.jl")
    include("Hypix/CheckError.jl")
    include("Hypix/Pet.jl")
    include("HyPix/Other/θaver.jl")
    include("Hypix/Memory.jl")
    include("Hypix/Climate.jl")

    include("Hypix/PlotHypix.jl")

    include("Hypix/HypixModel.jl")
    include("Hypix/Opt/HypixOpt.jl")
    include("Hypix/HypixStart.jl")

    # include("NoCore/Jules/Jules.jl")

    include("Temporary/KS_SMAP.jl")

end # Suppressor 
