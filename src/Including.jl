# =============================================================
#		module: including
# =============================================================
# module including

    using Suppressor
    using Revise

    # @suppress begin
        include("Option.jl")
        # Install packages to run program
        if option.globalopt.DownloadPackage
            include("Packages.jl")
        end # option.globalopt.DownloadPackage
        include("Tool.jl")
        include("Hypix\\Other\\Sitename.jl")
        include("Path.jl")
        include("Cst.jl")
        include("Param.jl")
        include("Hydro\\ΨminΨmax.jl")
        include("Hydro\\HydroStruct.jl")
        include("Hydro\\HydroRelation.jl")
        include("Hydro\\Wrc.jl")
        include("Optim\\Optimize.jl")
        include("Reading.jl")

        if !(option.globalopt.Hypix)
            include("Hydro\\Φ.jl")
        end
        include("Checking.jl")
        include("Hydro\\Kunsat.jl")
        include("Stats.jl")
        if !(option.globalopt.Hypix)
            include("Table.jl")
            include("Psd\\PsdThetar.jl")
        end

        if option.globalopt.Smap
            include("Smap\\StoneSmap.jl")
            include("Smap\\ReadSmap.jl")
            include("Smap\\PlotSmap.jl")
            include("Smap\\TableSmap.jl")
        end

        if option.globalopt.θΨ ≠ :No && option.globalopt.θΨ ≠ :File &&  !(option.globalopt.Hypix)
            # include("HydroLab\\HydrolabInitialize.jl")	
            include("HydroLab\\OfHydrolab.jl")
            # include("HydroLab\\START_Lab.jl")
            include("HydroLab\\HydrolabOpt.jl")
            include("Hypix\\Other\\PlotOther.jl")
        end
        
        if option.globalopt.Infilt
            include("Sorptivity\\Sorptivity.jl")            
            include("Infilt\\BestFunc.jl")
            include("Infilt\\OfBest.jl")
            include("Infilt\\QuasiExact.jl")
            include("Infilt\\TimeTransSteady.jl")
            include("Infilt\\InfiltStruct.jl")
            include("Infilt\\InfiltInitialize.jl")
            include("Infilt\\InfiltStart.jl")
        end # option.globalopt.Infilt
        
        if option.globalopt.Psd
            include("Psd\\PsdStruct.jl")
            include("Psd\\PsdInitialize.jl")
            include("Psd\\PsdFunc.jl")
            include("Psd\\PsdOpt.jl")
            include("Psd\\PsdStart.jl")
        end # option.globalopt.Psd

        if option.globalopt.Ploting && !(option.globalopt.Hypix)
            include("Plot.jl")
        end # option.globalopt.Ploting

        println(option.globalopt.Hypix)
        if option.globalopt.Hypix
            include("Hypix\\PathHypix.jl")
            include("Sorptivity\\Sorptivity.jl")
            include("Hypix\\Interpolate.jl")
            include("Hypix\\Opt\\ThetaObs.jl")
            include("HyPix\\θini.jl")
            # include("Hypix\\Opt\\Signature.jl")
            include("Hypix\\Opt\\OfHypix.jl")
            include("Hypix\\TableHypix.jl")
            include("Hypix\\VegStruct.jl")
            include("Hypix\\HorizonLayer.jl")
            include("Hypix\\ReadHypix.jl")
            include("Hypix\\Interception.jl")
            include("Hypix\\Flux.jl")
            include("Hypix\\Discretization.jl")
            include("Hypix\\Ponding.jl")
            include("Hypix\\Residual.jl")
            include("Hypix\\Δtchange.jl")
            include("Hypix\\TimeStep.jl")
            include("Hypix\\Richard.jl")
            include("Hypix\\WaterBalance.jl")
            include("Hypix\\Evaporation.jl")
            include("Hypix\\RootWaterUptake.jl")
            include("Hypix\\CheckError.jl")
            include("Hypix\\Pet.jl")
            include("HyPix\\Other\\θaver.jl")
            include("Hypix\\Memory.jl")
            include("Hypix\\Climate.jl")
            if option.globalopt.Ploting
                include("Hypix\\Other\\PlotOther.jl")
                include("Hypix\\PlotHypix.jl")
            end
            include("Hypix\\HypixModel.jl")
            include("Hypix\\Opt\\HypixOpt.jl")
            include("Hypix\\HypixStart.jl")
        end  # if: option.globalopt.Hypix

        if option.globalopt.Jules
            include("Hypix\\PathHypix.jl")
            include("Hypix\\VegStruct.jl")
            include("Hypix\\Discretisation.jl")
            include("Jules\\Jules.jl")
            include("HyPix\\ThetaIni.jl")
            include("Smap\\Smap2hypix.jl")		
        end  # if: option.Temporay
    # end # Suppressor 
