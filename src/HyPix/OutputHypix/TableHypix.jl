# =============================================================
#		module: tableHyix
# =============================================================
module tableHypix

   import ..cst, ..tool, ..wrc, ..kunsat
   import Dates: value, DateTime, year, month, day, hour, minute, second

   import CSV, Tables
   export TABLE_HYPIX

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : TABLE_HYPIX
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function TABLE_HYPIX(∑∑ΔSink,  ∑WaterBalanceη_Reduced, ∑ΔQ_Bot, CccBest,  Date_Reduced, discret, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, hydroHorizon, iMultistep, iNonConverge_iOpt, N_Layer,  Nit_Reduced, NiZ, NseBest, optionHypix, paramHypix, pathOutputHypix, SwcRoots, veg, WilmotBest, WofBest, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔRunTimeHypix, ΔSink_Reduced, ΔT_Average, θ_Reduced, θobs_Reduced, θsim_Aver, Ψ_Reduced)
     
			println("		=== === START: Table === ===")

         # Writing values of hydraulic parameters
         tableHypix.HYDRO(hydroHorizon, iMultistep, N_Layer, pathOutputHypix)

         # Writing values of veg parameters
         tableHypix.VEG(veg, iMultistep, pathOutputHypix)

         tableHypix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iMultistep, iMultistep, NseBest, paramHypix,  pathOutputHypix, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average)	
            

         if optionHypix.Table_Discretization
            tableHypix.DISCRETISATION_RRE(discret, NiZ, Z[1:NiZ], pathOutputHypix)
         end
         if optionHypix.Table_TimeSerie
            tableHypix.TIME_SERIES_DAILY(  ∑WaterBalanceη_Reduced[1:Nit_Reduced], Date_Reduced[1:Nit_Reduced], iMultistep, Nit_Reduced, pathOutputHypix, ΔEvaporation_Reduced[1:Nit_Reduced], ΔPet_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], ΔPr_Reduced[1:Nit_Reduced], ΔPrGross_Reduced[1:Nit_Reduced],  ΔQ_Reduced[1:Nit_Reduced, NiZ+1], ΔRunoff_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced])
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
         if optionHypix.θavr_RootZone && optionHypix.θobs
            tableHypix.θAVERAGE(Date_Reduced[1:Nit_Reduced], iMultistep, θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], pathOutputHypix)
         end

      println("		=== === END: Table === === \n")
      return nothing
      end  # function: TABLE_HYPIX
   # ------------------------------------------------------------------


   # ===================================================
   #          TimeStep daily
   # ===================================================
      function TIME_SERIES_DAILY(∑WaterBalanceη_Reduced, Date_Reduced, iMultistep, Nit_Reduced, pathOutputHypix, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Reduced, ΔPrGross_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced)

         Header =  ["iD" ,"Year" ,"Month" ,"Day" ,"Hour" ,"Minute" ,"Second" ,"ΔPrGross[mm]", "ΔPrSoil[mm]" , "ΔPet[mm]" ,"ΔSink[mm]","ΔTranspiration[mm]" ,"ΔEvaporation[mm]" ,"ΔRecharge[mm]" ,"Hpond[mm]" ,"ΔRunoff[mm]","∑WaterBalance_η[mm]"]

         Path = pathOutputHypix.Table_TimeSerie_Daily * "_" * string(iMultistep) * ".csv"

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

         CSV.write(Path, Tables.table([Id Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ΔPrGross_Reduced ΔPr_Reduced ΔPet_Reduced ΔSink_Reduced ΔSink_Reduced.-ΔEvaporation_Reduced ΔEvaporation_Reduced ΔQ_Reduced ΔPond_Reduced ΔRunoff_Reduced ∑WaterBalanceη_Reduced]), writeheader=true, header=Header, bom=true)
      return nothing
      end # Table  TIME_SERIES_DAILY
   #------------------------------------------------------


   # ===================================================
   #          θ
   # ===================================================
      function θ(Date_Reduced, θ_Reduced, Znode, iMultistep, pathHyPix)
         Path = pathHyPix.Table_θ * "_" * string(iMultistep) * ".csv"

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
            Header = ["θ[mm² mm⁻²]_Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm]"]
            Header = vcat(Header, string.(-Znode))
      
         CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ θ_Reduced]), writeheader=true, header=Header, bom=true)
      return nothing
      end  # Table θ
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θAVERAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θAVERAGE(Date_Reduced, iMultistep, θobs_Reduced, θsim_Aver, pathHyPix)
         Path = pathHyPix.Table_θaverage * string(iMultistep) * ".csv"

         Header = ["Id", "Year","Month","Day" ,"θobs_Aver", "θsim_Aver"]

         Id = 1:1:length(θsim_Aver)

         Year = year.(Date_Reduced)
         Month = month.(Date_Reduced)
         Day = day.(Date_Reduced)

         CSV.write(Path, Tables.table([Id Year Month Day θobs_Reduced θsim_Aver]), writeheader=true, header=Header, bom=true)
      return nothing			
      end # function: θAVERAGE
   #------------------------------------------------------


   # ===================================================
   #          Q
   # ===================================================
      function Q(Date_Reduced, ΔQ_Reduced, Z_Bottom, Znode, iMultistep, pathHyPix)	
         Path = pathHyPix.Table_Q * "_" * string(iMultistep) * ".csv"

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

         Header = ["ΔFlux[mm]_Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm]"]
         Header = vcat(Header, string.(-Znode))

         CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ΔQ_Reduced]), writeheader=true, header=Header, bom=true) 
      return nothing
      end  # function Q
   #------------------------------------------------------


   # ===================================================
   #          Ψ
   # ===================================================
      function Ψ(Date_Reduced, Ψ_Reduced, Znode, iMultistep, pathHyPix)
         Path = pathHyPix.Table_Ψ * "_" * string(iMultistep) * ".csv"

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

         Header = ["Ψ[mm] Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm"]
         Header = vcat(Header, string.(-Znode))

         CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ -Ψ_Reduced]), writeheader=true, header=Header, bom=true)
      return nothing
      end  # function Ψ
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : PERFORMACE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, iOpt, iMultistep, NseBest, paramHypix, pathHyPix, SwcRoots, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average)	
            
         iSim₀ = paramHypix.opt.iOptMultiStep_Start + iOpt - 1	

         Path = pathHyPix.Table_Performance * "_" * string(iSim₀) * ".csv"

         Header = ["Multisteps", "WofBest[mm]", "NseBest[-]" ,"CccBest[-]", "WilmotBest[-]" ,"∑∑ΔSink[mm]","∑ΔQ_Bot[mm]" , "SwcRoots[mm]" ,"Global_WaterBalance[mm]" ,"Global_WaterBalance_NormPr[-]", "ΔT_Average[second]","iNonConverge[count]" , "Efficiency[count]", "ΔRunTimeHypix[second]"]

         Id = 1:1:length(WofBest)

         CSV.write(Path, Tables.table([Id  WofBest NseBest CccBest WilmotBest ∑∑ΔSink ∑ΔQ_Bot SwcRoots Global_WaterBalance Global_WaterBalance_NormPr ΔT_Average iNonConverge_iOpt Efficiency ΔRunTimeHypix]), writeheader=true, header=Header, bom=true)
      return nothing
      end # function PERFORMACE
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θΨ
   # 		Tabular values of the hydroParam model
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θΨ(hydroHorizon, iMultistep, N_Layer, optionₘ, paramHypix, pathHyPix)		
         Path = pathHyPix.Table_θΨ * "_" * string(iMultistep) * ".csv"

         N_θΨobs = Int64(length(paramHypix.ploting.θΨ_Table))

         # Writting the Header
            FieldName_String = fill(""::String, N_θΨobs)

            for i =1:N_θΨobs
               FieldName_String[i] = string(-paramHypix.ploting.θΨ_Table[i] )
            end
            pushfirst!(FieldName_String, string("Ψ[mm] / θ(Ψ)[mm² mm⁻²]")) # Write the "Id" at the very begenning
         
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
            CSV.write(Path, Tables.table(θ_Mod[1:N_Layer, 1:N_θΨobs+1]), writeheader=true, header=FieldName_String, bom=true)

      return nothing	
      end  # function:  θΨK_PSD
   #------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : KΨ
   # 		Tabular values of the hydroParam model
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function KΨ(hydroHorizon, iMultistep, N_Layer, optionₘ, paramHypix, pathHyPix)				
         Path = pathHyPix.Table_KΨ * "_" * string(iMultistep) * ".csv"

         N_θΨobs = Int64(length(paramHypix.ploting.θΨ_Table))

         # Writting the Header
            FieldName_String = fill(""::String, N_θΨobs)

            for i =1:N_θΨobs
               FieldName_String[i] = string(paramHypix.ploting.θΨ_Table[i])
            end
            pushfirst!(FieldName_String, string("Ψ[mm] / K(Ψ)[mm/hour]")) # Write the "Id" at the very begenning
         
         # Computing θ at required θ
            K_Mod = fill(0.0::Float64, (N_Layer, N_θΨobs))
            for iZ=1:N_Layer, iΨ =1:N_θΨobs
               Ψ_Mod =paramHypix.ploting.θΨ_Table[iΨ]
               K_Mod[iZ, iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_Mod, iZ, hydroHorizon) .* cst.MmS_2_MmH
            end # iZ

         # Concatenating the 2 matrices
         Id = 1:1:N_Layer

         K_Mod = hcat(Id, K_Mod)

         # Writting the table
            CSV.write(Path, Tables.table(K_Mod[1:N_Layer,1:N_θΨobs+1]), writeheader=true, header=FieldName_String, bom=true)
      return nothing	
      end  # function:  θΨK_PSD
   #------------------------------------------------------


   # ===================================================
   #          DISCRETISATION AUTO
   # ===================================================
      function DISCRETISATION_AUTO(Flag_θΨini::Symbol, Layer::Vector{Int64}, PathDiscretisation::String, Z::Vector{Float64}, θini_or_Ψini_Cell::Vector{Float64})

         if Flag_θΨini == :Ψini
            Header = ["iZ","Z", "Layer", "Ψini"]

         elseif Flag_θΨini == :θini
            Header = ["iZ","Z", "Layer", "θini"]
         end

         iZ = collect(1:1:length(Z))

         CSV.write(PathDiscretisation, Tables.table([iZ Z Layer θini_or_Ψini_Cell]), writeheader=true, header=Header, bom=true)
      return nothing
      end # Table DISCRETISATION_AUTO
   #------------------------------------------------------


   # ===================================================
   #          Discretization
   # ===================================================
      function DISCRETISATION_RRE(discret, NiZ, Z, pathHyPix)
         Header =  ["Z", "ΔZ", "ΔZ_⬓", "Znode", "ΔZ_Aver", "ΔZ_W", "Z_CellUp"]

         CSV.write(pathHyPix.Table_Discretisation, Tables.table( [Z[1:NiZ] discret.ΔZ[1:NiZ] discret.ΔZ_⬓[1:NiZ] discret.Znode[1:NiZ] discret.ΔZ_Aver[1:NiZ] discret.ΔZ_W[1:NiZ] discret.Z_CellUp[1:NiZ]]), writeheader=true, header=Header, bom=true)
      return nothing
      end # Table DISCRETISATION
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function HYDRO(hydroHorizon, iMultistep, N_Layer, pathOutputHypix)
         Path = pathOutputHypix.Table_Hydro  * "_" * string(iMultistep) * ".csv"

         Id = 1:1:N_Layer

         Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_Layer, hydroHorizon)
               
         pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

         CSV.write(Path, Tables.table([Int64.(Id) Matrix]), writeheader=true, header=FieldName_String, bom=true)
      return nothing			
      end  # function: HYDRO
   #------------------------------------------------------
   

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : veg
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function VEG(veg, iMultistep, pathOutputHypix)
         Path = pathOutputHypix.Table_Veg * "_" * string(iMultistep) * ".csv"

         Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(1, veg)
       
         CSV.write(Path, Tables.table(Matrix), writeheader=true, header=FieldName_String, bom=true)
      return nothing
      end  # function: VEG
   #------------------------------------------------------

end  # module: tableHyix
# ............................................................