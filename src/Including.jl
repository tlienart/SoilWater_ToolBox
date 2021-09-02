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
    include("HyPix/HorizonLayer.jl")
    include("Hydro/HydroRelation.jl")
    include("Hydro/Wrc.jl")
    include("Hydro/Kunsat.jl")
    include("Stats.jl")
    include("Psd/PsdThetar.jl")
    include("Optim/Optimize.jl")
    include("Table.jl")
    include("Reading.jl")
    
    include("Distribution.jl")
    include("Ksmodel/θψ_2_KsModel.jl")
    include("Ksmodel/Opt_KsModel.jl")
    include("Ksmodel/Start_KsModel.jl")
    include("Ksmodel/Struct_Ksmodel.jl")
    
    include("Checking.jl")

    include("RockFragment/RockFragment.jl")
    
    include("HyPix/VegStruct.jl")
    include("HyPix/Discretization.jl")
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

    # include("HyPix/Interpolate.jl")
    # include("HyPix/Opt/ThetaObs.jl")
    # include("HyPix/θini.jl")
    # include("HyPix/Opt/OfHypix.jl")
    # include("HyPix/Interception.jl")
    # include("HyPix/Flux.jl")
    # include("HyPix/Ponding.jl")
    # include("HyPix/Residual.jl")
    # include("HyPix/ΔΔtchange.jl")
    # include("HyPix/TimeStep.jl")
    # include("HyPix/Richard.jl")
    # include("HyPix/WaterBalance.jl")
    # include("HyPix/Evaporation.jl")
    # include("HyPix/RootWaterUptake.jl")
    # include("HyPix/CheckError.jl")
    # include("HyPix/Pet.jl")
    # include("HyPix/Other/θaver.jl")
    # include("HyPix/Memory.jl")
    # include("HyPix/Climate.jl")
    # include("HyPix/PlotHypix.jl")
    # include("HyPix/HypixModel.jl")
    # include("HyPix/Opt/HypixOpt.jl")
    include("HyPix/HypixStart.jl")

    # include("NoCore/Jules/Jules.jl")
    # include("Temporary/Ks_Smap.jl")

    # include("NoCore/NSDR/ReadNsdr.jl")
end # Suppressor 