# =============================================================
#		module: including
# =============================================================
# module including

using Suppressor

@suppress begin
    include("Option.jl")
    include("Path.jl")

    # include("Packages.jl")
    # include("NaNTrap.jl")
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

    include("Plot.jl")
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

    include("Infilt/Infilt_START.jl")

    include("Psd/Psd_START.jl")

    include("HyPix/HypixStart.jl")

    # include("NoCore/Jules/Jules.jl")
    # include("Temporary/Ks_Smap.jl")

    # include("NoCore/NSDR/ReadNsdr.jl")
end # Suppressor 