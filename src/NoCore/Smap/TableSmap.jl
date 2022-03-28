# =============================================================
#		module: tableSmap
# =============================================================
module tableSmap
   import ..tool, ..wrc, ..kunsat
   import DelimitedFiles, CSV, Tables

   	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function θΨK(hydro, hydroOther, IdSelect, KₛModel, NiZ, Path, Smap_Depth, Soilname)
            println("    ~  $(Path) ~")

            Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydro)

            Matrix2, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydroOther)

            # Concatenating matrices
               Matrix = hcat(Matrix, Matrix2)

               Matrix = hcat(Matrix, KₛModel)

               FieldName_String = vcat(FieldName_String, FieldName_String2)

               FieldName_String = vcat("Id", "SoilName", "Depth", FieldName_String, "Ks_Model[mm/s]")

                open(Path, "w") do io
                  DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
                  DelimitedFiles.writedlm(io, [string.(IdSelect[1:NiZ]) Soilname[1:NiZ] Smap_Depth[1:NiZ] Matrix], ",")
               end
         return nothing
         end  # function:  θΨK

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : Smap
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   """
   TopnetModel = ["ThetaS_Ch[mm3 mm-3]";"ThetaR_Ch[mm3 mm-3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]";" Hga_Ch[mm]";"Ks_Vg[mm s-1]; "0mm"; "500mm"; "1000mm"; "2000mm"; "4000mm"; "10000mm"; "150000mm"]

   JulesModel_CH =  ["ThetaS_Ch[mm3 mm-3]";"ThetaR_Ch[mm3 mm-3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]"; "Ks_Ch[mm s-1]";" Hga_Ch[mm]";"3300mm";"10000mm" ;"150000mm"]

   JulesModel_VangenuchtenJules = ["ThetaS_VgJules[mm3 mm-3]";"ThetaR_VgJules[mm3 mm-3]";"n_VgJules[-]";"Hvg_VgJules[mm]"; "Ks_VgJules[mm s-1]";"3300mm";"10000mm"]

   """
      function SMAP(hydro, IdSelect, IsTopsoil, NiZ, optionₘ, param, path, RockFragment, Smap_Depth, Smap_MaxRootingDepth, Smap_PermeabilityClass, Smap_RockDepth, Smap_SmapFH, Soilname)

         println("    ~  $(path.tableSmap.Table_Smap) ~")

         HeaderSmap = true # <true> the greek characters are replaced by alphabet; <false> original parameter names with no units usefull to use values in SoilWater-ToolBox

         # User input
            Option_BrooksCorey       = false
            Option_ClappHornberger   = false
            Option_VanGenuchten      = false
            Option_VanGenuchtenJules = false
            Option_Kosugi            = true
            Option_Kosugi_Table_θψ   = true
            Option_Kosugi_Table_Kψ   = true

            

         Header = ["Id"; "SoilName"; "Depth_mm"; "IsTopsoil"; "RockFragment_%";"RockDepth_mm"; "MaxRootingDepth_mm"; "PermeabilityClass"; "SmapFH"]
         Data = []
      
      # Select data
         # HydroModel_θΨ == "BrooksCorey"
         if Option_BrooksCorey # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "BrooksCorey"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"
            
            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"λbc";"Ψbc"; "Ks"; "Ψga"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
               
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Bc[mm3 mm-3]";"ThetaR_Bc[mm3 mm-3]";"LambdaBc_Bc[-]";"Hbc_Bc[mm]";"Ks_Bc[mm s-1]";"Hga_Bc[mm]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_BrooksCorey


         # HydroModel_θΨ == "ClappHornberger"
         if  Option_ClappHornberger # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "ClappHornberger"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"λch";"Ψch";"Ks";"Ψga"]

               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end

               if HeaderSmap
                  Header_θΨ = ["ThetaS_Ch[mm3 mm-3]";"ThetaR_Ch[mm3 mm-3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]"; "Ks_Ch[mm s-1]";" Hga_Ch[mm]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end #  Option_ClappHornberger

         
         # HydroModel_θΨ == "Vangenuchten"
         if Option_VanGenuchten # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "Vangenuchten"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Vg[mm3 mm-3]";"ThetaR_Vg[mm3 mm-3]";"N_Vg[-]";"Hvg_Vg[mm]"; "Ks_Vg[mm s-1]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_VanGenuchten


         # HydroModel_θΨ == "VangenuchtenJules"
         if Option_VanGenuchtenJules # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "VangenuchtenJules"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_VgJules[mm3 mm-3]";"ThetaR_VgJules[mm3 mm-3]";"n_VgJules[-]";"Hvg_VgJules[mm]"; "Ks_VgJules[mm s-1]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_VanGenuchten


         # HydroModel_θΨ == "Kosugi"
         if Option_Kosugi # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "Kosugi"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ =["θs";"θr";"Ks";"Ψm";"σ";"θsMacMat_ƞ";"σMac";"ΨmMac"; "θsMacMat";"Φ"]
                   
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Kg[mm3 mm-3]";"ThetaR_Kg[mm3 mm-3]";"Ks_Kg[mm s-1]";"Hm_Kg[mm]";"Sigma_Kg";"ThetaSMacMatNorm_Kg[mm3 mm-3]";"SigmaMac_Kg";"HmMac_Kg[mm]";"θsMacMat_Kg[mm3 mm-3]";"TotalPorosity[mm3 mm-3]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_Kosugi


      # HydroModel_θΨ == "Option_Kosugi_Table_θψ"
      if Option_Kosugi_Table_θψ # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
         HydroModel_θΨ = "Kosugi"

         N_Ψ = length(param.smap.Ψ_Table[:])
         θ₂ = fill(0.0::Float64, (NiZ, N_Ψ))

         for iZ=1:NiZ
            for iΨ =1:N_Ψ
               Ψ₂ = param.smap.Ψ_Table[iΨ]
               θ₂[iZ, iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψ₂, iZ, hydro)
            end # iΨ
         end # iZ

         if isfile(Path_θΨ)
            Select_θΨ = "ThetaH_" .* string.(Int64.(param.smap.Ψ_Table)) .* "_mm"            

            Data = hcat(Data[1:NiZ, :], θ₂[1:NiZ, :])
      
            Header_θΨ = Select_θΨ

            Header =  append!(Header, Header_θΨ)
         end
      end # Option_Kosugi



      # HydroModel_θΨ == "Option_Kosugi_Table_Kψ"
      if Option_Kosugi_Table_Kψ # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
         HydroModel_θΨ = "Kosugi"

         N_Ψ = length(param.hydro.K_Table[:])
         K₂ = fill(0.0::Float64, (NiZ, N_Ψ))

         for iZ=1:NiZ
            for iΨ =1:N_Ψ
               Ψ₂ = param.hydro.K_Table[iΨ]
               K₂[iZ, iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ₂, iZ, hydro)
            end # iΨ
         end # iZ

         if isfile(Path_θΨ)
            Select_θΨ = "KunsatH_" .* string.(Int64.(param.hydro.K_Table)) .* "_mm s-1"            

            Data = hcat(Data[1:NiZ, :], K₂[1:NiZ, :])
      
            Header_θΨ = Select_θΨ

            Header =  append!(Header, Header_θΨ)
         end
      end # Option_Kosugi

         
      # COMBINING OUTPUTS   
         open(path.tableSmap.Table_Smap, "w") do io
            DelimitedFiles.writedlm(io,[Header] , ",",) # Header
            DelimitedFiles.writedlm(io, [string.(IdSelect[1:NiZ]) Soilname[1:NiZ] Smap_Depth[1:NiZ] IsTopsoil[1:NiZ] RockFragment[1:NiZ] Smap_RockDepth[1:NiZ] Smap_MaxRootingDepth[1:NiZ] Smap_PermeabilityClass[1:NiZ] Smap_SmapFH[1:NiZ] Data[1:NiZ,:]], ",")
         end	
      return nothing
      end  # function:  smap
	# ............................................................
  
end  # module: tableSmap
# ............................................................