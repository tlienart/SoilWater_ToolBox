# =============================================================
#		module: including
# =============================================================
# module including

    using Suppressor
    using Revise

    @suppress begin
        include("Option.jl")
            # Reading the options
                option = options.OPTIONS()

        # Install packages to run program
        if option.other.DownloadPackage
            include("Packages.jl")
        end # option.other.DownloadPackage
        include("Tool.jl")
        include("Path.jl")
        include("Cst.jl")
        include("Param.jl")
        include("Hydro/ΨminΨmax.jl")
        include("Hydro/HydroStruct.jl")
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

        if !(option.run.Hypix)
            include("Psd/PsdThetar.jl")
        end

        if option.run.Smap
            include("NoCore/Smap/ReadSmap.jl")
            include("NoCore/Smap/TableSmap.jl")
            if option.other.Ploting
                include("NoCore/Smap/PlotSmap.jl")
            end
        end

        if option.run.HydroLabθΨ⍰ ≠ :No && option.run.HydroLabθΨ⍰ ≠ :File &&  !(option.run.Hypix)
            include("HydroLab/OfHydrolab.jl")
            include("HydroLab/HydrolabOpt.jl")
        end
        
        if option.run.InfiltBest
            include("Sorptivity/Sorptivity.jl")            
            include("Infilt/BestFunc.jl")
            include("Infilt/OfBest.jl")
            include("Infilt/QuasiExact.jl")
            include("Infilt/TimeTransSteady.jl")
            include("Infilt/InfiltStruct.jl")
            include("Infilt/InfiltInitialize.jl")
            include("Infilt/InfiltStart.jl")
        end # option.run.InfiltBest
        
        if option.run.IntergranularMixingPsd
            include("Psd/PsdStruct.jl")
            include("Psd/PsdInitialize.jl")
            include("Psd/PsdFunc.jl")
            include("Psd/PsdOpt.jl")
            include("Psd/PsdStart.jl")
        end # option.run.IntergranularMixingPsd

        if option.other.Ploting && !(option.run.Hypix)
            include("Plot.jl")
        end # option.other.Ploting

        if option.run.Hypix
            include("Sorptivity/Sorptivity.jl")
            include("Hypix/Interpolate.jl")
            include("Hypix/Opt/ThetaObs.jl")
            include("HyPix/θini.jl")
            include("Hypix/Opt/OfHypix.jl")
            include("Hypix/VegStruct.jl")
            include("Hypix/HorizonLayer.jl")
            include("Hypix/Interception.jl")
            include("Hypix/Flux.jl")
            include("Hypix/Discretization.jl")
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
            if option.other.Ploting
                include("Hypix/PlotHypix.jl")
            end
            include("Hypix/HypixModel.jl")
            include("Hypix/Opt/HypixOpt.jl")
            include("Hypix/HypixStart.jl")
        end  # if: option.run.Hypix

        if option.dataFrom.Jules
            include("Hypix/VegStruct.jl")
            include("Hypix/Discretization.jl")
            include("NoCore/Jules/Jules.jl")
            include("HyPix/θini.jl")
            include("NoCore/Smap/Smap2hypix.jl")		
        end  # if: option.Temporay
    end # Suppressor 
