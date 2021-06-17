# =============================================================
#		module: tableSmap
# =============================================================
module tableSmap
   import ..tool, ..wrc
   import CSV, Tables, DataFrames, DelimitedFiles

   	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function θΨK(hydro, hydroOther, IdSelect, Kₛ_Model, N_iZ, smap, Path)
            println("    ~  $(Path) ~")

            Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ, hydro)

            Matrix2, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ, hydroOther)

            # Concatenating matrices
               Matrix = hcat(Matrix, Matrix2)

               Matrix = hcat(Matrix, Kₛ_Model)

               FieldName_String = vcat(FieldName_String, FieldName_String2)

               pushfirst!(FieldName_String, string("Depth"))
               pushfirst!(FieldName_String, string("SoilName"))
               pushfirst!(FieldName_String, string("Id")) 
               push!(FieldName_String, string("KunsatModel"))
               
            CSV.write(Path, Tables.table([IdSelect smap.Soilname[1:N_iZ] smap.Depth[1:N_iZ] Matrix]), header=FieldName_String )
      
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
      function SMAP(optionₘ, IdSelect, N_iZ, smap, param, path)
         println("    ~  $(path.Table_Smap) ~")

         # Header
            HeaderSmap = false # <true> the greek characters are replaced by alphabet; <false> original parameter names with no units usefull to use values in SoilWater-ToolBox

         # User input
            Option_BrooksCorey       = true
            Option_ClappHornberger   = true
            Option_VanGenuchten      = true
            Option_VanGenuchtenJules = true
            Option_Kosugi            = true
            Option_Kosugi_Table      = true

         Header = ["Id"; "SoilName"; "Depth_mm"; "IsTopsoil"; "RockFragment_%";"RockDepth_mm"; "MaxRootingDepth_mm"]

         Data = []
      
      # Select data
         # HydroModel_θΨ == "BrooksCorey"
         if Option_BrooksCorey # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "BrooksCorey"

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(optionₘ.σ_2_Ψm) *  "_" * path.option.ModelName * "_" * path.tableSmap.Table_θΨK₀
            
            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"λbc";"Ψbc"; "Ks"; "Ψga"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:N_iZ, :], Data_θΨ[1:N_iZ, :])
               catch
                  Data = Data_θΨ[1:N_iZ, :]
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

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(optionₘ.σ_2_Ψm) *  "_" * path.option.ModelName * "_" * path.tableSmap.Table_θΨK₀ 

            if isfile(Path_θΨ)

               Select_θΨ = ["θs";"θr";"λch";"Ψch";"Ks";"Ψga"]

               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:N_iZ, :], Data_θΨ[1:N_iZ, :])
               catch
                  Data = Data_θΨ[1:N_iZ, :]
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

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(optionₘ.σ_2_Ψm) *  "_" * path.option.ModelName * "_" * path.tableSmap.Table_θΨK₀ 

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:N_iZ, :], Data_θΨ[1:N_iZ, :])
               catch
                  Data = Data_θΨ[1:N_iZ, :]
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

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(optionₘ.σ_2_Ψm) *  "_" * path.option.ModelName * "_" * path.tableSmap.Table_θΨK₀ 

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:N_iZ, :], Data_θΨ[1:N_iZ, :])
               catch
                  Data = Data_θΨ[1:N_iZ, :]
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

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(optionₘ.σ_2_Ψm) *  "_" * path.option.ModelName * "_" * path.tableSmap.Table_θΨK₀ 

            if isfile(Path_θΨ)
               Select_θΨ =["θs";"θr";"Ks";"Ψm";"σ";"σMac";"ΨmMac";"θsMacMat"; "θsMacMat_ƞ"]
                   
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:N_iZ, :], Data_θΨ[1:N_iZ, :])
               catch
                  Data = Data_θΨ[1:N_iZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Kg[mm3 mm-3]";"ThetaR_Kg[mm3 mm-3]";"Ks_Kg[mm s-1]";"Hm_Kg[mm]";"Sigma_Kg";"SigmaMac_Kg";"HmMac_Kg[mm]";"ThetaSMacMat_Kg[mm3 mm-3]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_Kosugi


      # HydroModel_θΨ == "Option_Kosugi_Table"
      if Option_Kosugi_Table # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
         # HydroModel_θΨ = "Kosugi"

         Path_θΨ =  path.tableSoilwater.Table_KosugiθΨ

         if isfile(Path_θΨ)
            Select_θΨ = string.(Int64.(param.hydro.smap.Ψ_Table)) .* "mm"
            
            Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
               
            try
               Data = hcat(Data[1:N_iZ, :], Data_θΨ[1:N_iZ, :])
            catch
               Data = Data_θΨ[1:N_iZ, :]
            end
      
            Header_θΨ = Select_θΨ

            Header =  append!(Header, Header_θΨ)
         else
            @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
         end # if isfile(Path_θΨ)
      end # Option_Kosugi

         
      # COMBINING OUTPUTS   
         Output = Tables.table([string.(IdSelect[1:N_iZ]) smap.Soilname[1:N_iZ] smap.Depth[1:N_iZ] smap.IsTopsoil[1:N_iZ] smap.RockFragment[1:N_iZ] smap.RockDepth[1:N_iZ] smap.MaxRootingDepth[1:N_iZ] Data[1:N_iZ,:]])
      
         CSV.write(path.Table_Smap, Output, header=Header)	

      return nothing
      end  # function:  smap
	# ............................................................

  
end  # module: tableSmap
# ............................................................