# =============================================================
#		module: tableHyix
# =============================================================
module tableHypix

   import ..cst, ..tool, ..wrc, ..kunsat
   import DelimitedFiles
   import Dates: value, DateTime, year, month, day, hour, minute, second
   export TABLE_HYPIX

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : TABLE_HYPIX
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function TABLE_HYPIX(∑∑ΔSink, ∑Pr, ∑T, ∑T_Climate, ∑T_Reduced, ∑WaterBalance_η, ∑WaterBalanceη_Reduced, ∑ΔQ_Bot, CccBest, clim, Date_Reduced, discret, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, Hpond, hydroHorizon, iMultistep, iNonConverge_iOpt, iScenario, N_Layer, Nit, Nit_Reduced, NiZ, NseBest, optionHypix, paramHypix, pathOutputHypix, Q, SwcRoots, veg, WilmotBest, WofBest, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr, ΔPr_Reduced, ΔQ_Reduced, ΔRunTimeHypix, ΔSink_Reduced, ΔT, ΔT_Average, θ_Reduced, θobs_Reduced, θsim_Aver, Ψ_Reduced)
     
			println("		=== === START: Table === ===")

         # Writing values of hydraulic parameters
         tableHypix.HYDRO(hydroHorizon, iMultistep, N_Layer, pathOutputHypix)

         # Writing values of veg parameters
         tableHypix.VEG(veg, iMultistep, pathOutputHypix)

         tableHypix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iMultistep, iScenario, NseBest, paramHypix,  pathOutputHypix, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average)	
            

         if optionHypix.Table_Discretization
            tableHypix.DISCRETISATION_RRE(discret, NiZ, Z[1:NiZ], pathOutputHypix)
         end
         if optionHypix.Table_TimeSeriesDaily
            tableHypix.TIME_SERIES_DAILY(∑T_Reduced[1:Nit_Reduced], ∑WaterBalanceη_Reduced[1:Nit_Reduced], Date_Reduced[1:Nit_Reduced], iMultistep, Nit_Reduced, ΔEvaporation_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced, NiZ+1], ΔPet_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], ΔPr_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced], pathOutputHypix)
         end
         if optionHypix.Table_θ
            tableHypix.θ(Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,1:NiZ], discret.Znode[1:NiZ], iMultistep, pathOutputHypix)
         end
         if optionHypix.Table_Ψ
            tableHypix.Ψ(Date_Reduced[1:Nit_Reduced], Ψ_Reduced[1:Nit_Reduced,1:NiZ], discret.Znode[1:NiZ], iMultistep, pathOutputHypix)
         end
         if optionHypix.Table_Q
            tableHypix.Q(Date_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced,1:NiZ+1], Z[NiZ], discret.Znode[1:NiZ], iMultistep, pathOutputHypix)
         end
         if optionHypix.Tabule_θΨ
            tableHypix.θΨ(hydroHorizon, iMultistep, N_Layer, optionHypix, paramHypix, pathOutputHypix)
            tableHypix.KΨ(hydroHorizon, iMultistep, N_Layer, optionHypix, paramHypix, pathOutputHypix)
         end
         if optionHypix.Table_Climate
            tableHypix.DAILY_CLIMATE(∑T_Climate, clim, iMultistep, pathOutputHypix)
         end
         if optionHypix.θavr_RootZone && optionHypix.θobs
            tableHypix.θAVERAGE(Date_Reduced[1:Nit_Reduced], iMultistep, θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], pathOutputHypix)
         end

      println("		=== === END: Table === === \n")
      return nothing
      end  # function: TABLE_HYPIX
   # ------------------------------------------------------------------

   # ===================================================
   #          DISCRETISATION AUTO
   # ===================================================
      function DISCRETISATION_AUTO(Flag_θΨini::Symbol, Layer::Vector{Int64}, PathDiscretisation::String, Z::Vector{Float64}, θini_or_Ψini_Cell::Vector{Float64})

         # println("			~  $(PathDiscretisation) ~")

         if Flag_θΨini == :Ψini
            Header = ["iZ";"Z"; "Layer"; "Ψini"]

         elseif Flag_θΨini == :θini
            Header = ["iZ";"Z"; "Layer"; "θini"]
         end

         iZ = collect(1:1:length(Z))

         open(PathDiscretisation, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",") # Header
            DelimitedFiles.writedlm(io, [iZ Z Layer θini_or_Ψini_Cell], ",")
         end # open

      return nothing
      end # Table DISCRETISATION_AUTO
   #------------------------------------------------------


   # ===================================================
   #          Discretization
   # ===================================================
      function DISCRETISATION_RRE(discret, NiZ, Z, pathHyPix)
         println("			~  $(pathHyPix.Table_Discretisation) ~")

         Header =  ["Z" "ΔZ" "ΔZ_⬓" "Znode" "ΔZ_Aver" "ΔZ_W" "Z_CellUp"]

         open(pathHyPix.Table_Discretisation, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [Z[1:NiZ] discret.ΔZ[1:NiZ] discret.ΔZ_⬓[1:NiZ] discret.Znode[1:NiZ] discret.ΔZ_Aver[1:NiZ] discret.ΔZ_W[1:NiZ] discret.Z_CellUp[1:NiZ]], ",")
         end
      return nothing
      end # Table DISCRETISATION
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function HYDRO(hydroHorizon, iScenario, N_Layer, pathOutputHypix)
         Path = pathOutputHypix.Table_Hydro  * "_" * string(iScenario) * ".csv"
         println("			~ $(Path) ~")

         Id = 1:1:N_Layer

         Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_Layer, hydroHorizon)
               
         pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

         open(Path, "w") do io
            DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
            DelimitedFiles.writedlm(io, [Int64.(Id) Matrix], ",")
         end
      return nothing			
      end  # function: HYDRO
   #------------------------------------------------------
      

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : veg
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function VEG(veg, iScenario, pathHyPix)
         Path = pathHyPix.Table_Veg * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(1, veg)

         open(Path, "w") do io
            DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
            DelimitedFiles.writedlm(io, [Matrix], ",")
         end
      return nothing
      end  # function: VEG
   #------------------------------------------------------


   # ===================================================
   #          TimeStep at ΔT
   # ===================================================
      function TIME_SERIES(∑T, ΔT, ∑Pr, ΔPr, Hpond, Recharge, ∑WaterBalance_η, iScenario, pathHyPix)		
         Path = pathHyPix.Table_TimeSerie * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")
         
         Header =  ["∑T[mm]" "ΔT[mm]" "∑Pr[mm/ΔT]" "ΔPr[mm/ΔT]" "Hpond[mm]" "Recharge[mm/ΔT]" "∑WaterBalance_η[mm]"]

         open(Path, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [∑T ΔT ∑Pr ΔPr Hpond Recharge ∑WaterBalance_η], ",")
         end
      return nothing
      end # Table DISCRETISATION
   #------------------------------------------------------


   # ===================================================
   #          TimeStep daily
   # ===================================================
      function TIME_SERIES_DAILY(∑T_Reduced, ∑WaterBalanceη_Reduced, Date_Reduced, iScenario, Nit_Reduced, ΔEvaporation_Reduced, ΔRecharge_Plot, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔSink_Reduced, pathHyPix)
         Header =  ["iD" "Year" "Month" "Day" "Hour" "Minute" "Second" "∑T[Hour]" "ΔPr_Through[mm/day]" "ΔPet[mm/day]" "ΔSink[mm/day]" "ΔEvaporation[mm/day]" "Hpond[mm]" "Recharge[mm/day]" "∑WaterBalance_η_Profile[mm/day]"]

         Path = pathHyPix.Table_TimeSerie_Daily * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         Id = 1:1:Nit_Reduced

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         open(Path, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [Id Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ∑T_Reduced ΔPr_Reduced ΔPet_Reduced ΔSink_Reduced ΔEvaporation_Reduced ΔPond_Reduced ΔRecharge_Plot ∑WaterBalanceη_Reduced], ",")
         end
      return nothing
      end # Table  TIME_SERIES_DAILY
   #------------------------------------------------------


   # ===================================================
   #          θ
   # ===================================================
      function θ(Date_Reduced, θ_Reduced, Znode, iScenario, pathHyPix)
         Path = pathHyPix.Table_θ * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         # Adding an other column
            Header = ["Year" "Month" "Day" "Hour" "Minute" "Second"]

         DelimitedFiles.writedlm(Path, [Header -transpose(Znode); Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ θ_Reduced], ",")
      return nothing
      end  # Table θ
   #------------------------------------------------------


   # ===================================================
   #          Q
   # ===================================================
      function Q(Date_Reduced, ΔQ_Reduced, Z_Bottom, Znode, iScenario, pathHyPix)	
         Path = pathHyPix.Table_Q * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end
         
         # Adding an other column
         append!(Znode, Z_Bottom)

         Header = ["Year" "Month" "Day" "Hour" "Minute" "Second"]

         DelimitedFiles.writedlm(Path, [Header transpose(Znode); Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ΔQ_Reduced], ",")
      return nothing
      end  # function Q
   #------------------------------------------------------


   # ===================================================
   #          Ψ
   # ===================================================
      function Ψ(Date_Reduced, Ψ_Reduced, Znode, iScenario, pathHyPix)
         Path = pathHyPix.Table_Ψ * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         Header = ["Year" "Month" "Day" "Hour" "Minute" "Second"]

         DelimitedFiles.writedlm(Path, [Header -transpose(Znode); Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ -Ψ_Reduced], ",")
      return nothing
      end  # function Ψ
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θΨ
   # 		Tabular values of the hydroParam model
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θΨ(hydroHorizon, iScenario, N_Layer, optionₘ, paramHypix, pathHyPix)		
         Path = pathHyPix.Table_θΨ * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         N_θΨobs = Int64(length(paramHypix.ploting.θΨ_Table))

         # Writting the Header
            FieldName_String = fill(""::String, N_θΨobs)

            for i =1:N_θΨobs
               FieldName_String[i] = string(paramHypix.ploting.θΨ_Table[i] * cst.Mm_2_Cm) * "cm"
            end
            pushfirst!(FieldName_String, string("Layer")) # Write the "Id" at the very begenning
         
         # Computing θ at required θ
            θ_Mod = fill(0.0::Float64, (N_Layer, N_θΨobs))
            for iZ=1:N_Layer, iΨ =1:N_θΨobs
                  Ψ_Mod =paramHypix.ploting.θΨ_Table[iΨ]
                  θ_Mod[iZ, iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψ_Mod, iZ, hydroHorizon)
            end # iZ

         # Concatenating the 2 matrices
         Id = 1:1:N_Layer

         θ_Mod = hcat(Id, θ_Mod)

         # Writting the table
            open(Path, "w") do io
               DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
               for iZ = 1:N_Layer
                  # DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
                  DelimitedFiles.writedlm(io, [θ_Mod[iZ, 1:N_θΨobs+1]], ",")
               end # i
            end # Path
      return nothing	
      end  # function:  θΨK_PSD
   #------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : KΨ
   # 		Tabular values of the hydroParam model
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function KΨ(hydroHorizon, iScenario, N_Layer, optionₘ, paramHypix, pathHyPix)				
         Path = pathHyPix.Table_KΨ * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         N_θΨobs = Int64(length(paramHypix.ploting.θΨ_Table))

         # Writting the Header
            FieldName_String = fill(""::String, N_θΨobs)

            for i =1:N_θΨobs
               FieldName_String[i] = string(paramHypix.ploting.θΨ_Table[i] * cst.Mm_2_Cm) * "cm"
            end
            pushfirst!(FieldName_String, string("Layer Cm/H")) # Write the "Id" at the very begenning
         
         # Computing θ at required θ
            K_Mod = fill(0.0::Float64, (N_Layer, N_θΨobs))
            for iZ=1:N_Layer, iΨ =1:N_θΨobs
                  Ψ_Mod =paramHypix.ploting.θΨ_Table[iΨ]
                  K_Mod[iZ, iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_Mod, iZ, hydroHorizon) .* cst.MmS_2_CmH
            end # iZ

         # Concatenating the 2 matrices
         Id = 1:1:N_Layer

         K_Mod = hcat(Id, K_Mod)

         # Writting the table
            open(Path, "w") do io
               DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
               for iZ = 1:N_Layer
                  DelimitedFiles.writedlm(io, [K_Mod[iZ,1:N_θΨobs+1]], ",")
               end # i
            end # Path
      return nothing	
      end  # function:  θΨK_PSD
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : PERFORMACE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iOpt, iScenario, NseBest, paramHypix, pathHyPix, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average)	
         
         iSim₀ = paramHypix.iOptMultiStep_Start + iOpt	

         Path = pathHyPix.Table_Performance * "_" * string(iSim₀) * ".csv"
         println("			~  $(Path) ~")

         Header = ["Id" "WofBest" "NseBest" "CccBest" "WilmotBest" "Efficiency" "Global_WaterBalance" "Global_WaterBalance_NormPr" "ΔT_Average" "∑∑ΔSink" "∑ΔQ_Bot" "SwcRoots" "iNonConverge" "ΔRunTimeHypix"]

         Id = 1:1:length(WofBest)

         open(Path, "w") do io
            # DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
            DelimitedFiles.write(io, [0xef,0xbb,0xbf]) 
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [Id  WofBest NseBest CccBest WilmotBest Efficiency Global_WaterBalance Global_WaterBalance_NormPr ΔT_Average ∑∑ΔSink ∑ΔQ_Bot SwcRoots iNonConverge_iOpt ΔRunTimeHypix], ",")
         end
      return nothing
      end # function PERFORMACE
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : DAILY_CLIMATE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function DAILY_CLIMATE(∑T_Climate, clim, iScenario, pathHyPix)
         Path = pathHyPix.Table_DailyClimate * "_" * string(iScenario) * ".csv"
         println("			~  $(Path) ~")

         local ∑T_Int = ceil.(Int, ∑T_Climate[1:clim.N_Climate] .* cst.Second_2_Day)

         Header = ["Year" "Month" "Day" "Pr" "Pr_Ground"]

         open(Path, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [year.(clim.Date[1:clim.N_Climate]) month.(clim.Date[1:clim.N_Climate]) day.(clim.Date[1:clim.N_Climate]) clim.Pr[1:clim.N_Climate] clim.Pr_Through] , ",")
         end
      return nothing	
      end  # function: DAILY_CLIMATE
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θAVERAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θAVERAGE(Date_Reduced, iScenario, θobs_Reduced, θsim_Aver, pathHyPix)
         Path = pathHyPix.Table_θaverage * ".csv"
         println("			~  $(Path) ~")

         Header = ["Id", "Year","Month","Day" ,"θobs_Aver", "θsim_Aver"]

         Id = 1:1:length(θsim_Aver)

         Year = year.(Date_Reduced)
         Month = month.(Date_Reduced)
         Day = day.(Date_Reduced)

         open(Path, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [Id Year Month Day θobs_Reduced θsim_Aver] , ",")
         end # open
      return nothing			
      end # function: θAVERAGE
   #------------------------------------------------------

end  # module: tableHyix
# ............................................................