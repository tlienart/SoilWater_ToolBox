# =============================================================
#		module: smap2hypix
# =============================================================
module smap2hypix
   import ..tool, ..cst, ..discretisation, ..hydroStruct, ..reading, ..vegStruct, ..wrc
   import DelimitedFiles, Tables, CSV

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SMAP_2_HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SMAP_2_HYPIX(NiZ, optionₘ, param, path, Smap_Depth, Smap_MaxRootingDepth, Soilname)
         # Reducing Ks
            Ks_FactorReduce = 0.1

         # Hydraulic structures
            hydroSmap = hydroStruct.HYDROSTRUCT(optionₘ, NiZ) # Making a structure

         # Index of each soil profile in the spreadsheet
            iSoilProfile_End, iSoilProfile_Start, N_SoilProfile, Soilname_SoilProfile = SMAP_SOILPROFILE(NiZ, Soilname)

         # Input path of the hydraulic parameters
            Path_Input =  path.tableSoilwater.Path_Soilwater_Table *  "_" * "Kosugi" *  "_" * "Table_SmapThetaHK.csv"

         # For every soil profile
         for iSoilProfile=1:N_SoilProfile

            # Making directory
               Path_Output = path.smap2Hypix.Path_Smap2Hypix * "/" * Soilname[iSoilProfile] 
               mkpath(Path_Output)

            # Reading the hydraulic parameters of each soil profile
               hydroSmap, N_iSoilProfile = tool.readWrite.READ_STRUCT(hydroSmap, Path_Input; iStart=iSoilProfile_Start[iSoilProfile], iEnd=iSoilProfile_End[iSoilProfile])

               Path_Output₂ = Path_Output * "/HypixHydro" * "_" * Soilname[iSoilProfile] * ".csv"
               TABLE_HYDRO_VEG(hydroSmap, N_iSoilProfile, Path_Output₂)
            
            # Getting the depth of soils
                Smap_Depth_Soil = Smap_Depth[iSoilProfile_Start[iSoilProfile]:iSoilProfile_End[iSoilProfile]]

            # Discretisation
               LAYER_DISCRETISATION(param, Path_Output, Smap_Depth_Soil, Soilname[iSoilProfile])
            
            # Vegetation parameter
               vegSmap = vegStruct.VEGSTRUCT()

               # Abstracting data
               Path_Vegetaion ="D:\\DATAraw\\JULESdata\\Vegetation\\Vegetation.csv"
               
               vegSmap, NiZ = tool.readWrite.READ_STRUCT(vegSmap, Path_Vegetaion; iStart=1, iEnd=1)

               vegSmap.Zroot = min(vegSmap.Zroot, Smap_MaxRootingDepth[1])

               Path_Output₂ = Path_Output * "/Vegetation" * "_" * Soilname[iSoilProfile] * ".csv"

               TABLE_HYDRO_VEG(vegSmap, 1, Path_Output₂)
         end #  iSoilProfile=1:N_SoilProfile      
   
      return nothing
      end  # function: SMAP_2_HYDRO

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #		FUNCTION : DISCRETISATION
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function LAYER_DISCRETISATION(param, Path_Output, Smap_Depth_Soil, Soilname₀)
            iLayers = collect(1:1:length(Smap_Depth_Soil))

            Path_Output₂ = Path_Output * "/Layer" * "_" * string(Soilname₀) * ".csv"
            TABLE_DISCRETISATION(Path_Output₂, Smap_Depth_Soil, iLayers)
                  
            # Automatic Disscretizing of SoilProfiles per soil =====
               SoilProfile, Z = discretisation.DISCRETISATION_AUTO(param, N_Layer=length(Smap_Depth_Soil), Zlayer=Smap_Depth_Soil, Zroot=800.0)

               Path_Output₃ = Path_Output * "/Discretization_2" * "_" * string(Soilname₀) * ".csv"
               TABLE_DISCRETISATION(Path_Output₃, SoilProfile, Z)

            return nothing
         end  # function: DISCRETISATION


      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #		FUNCTION : SMAP_SOIL_iSoilProfile
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function SMAP_SOILPROFILE(NiZ, Soilname)
            iSoilProfile_End = []
            iSoilProfile_Start = [1]
            Soilname_Initial = Soilname[1]
            Soilname_SoilProfile = [Soilname[1]]
            i = 1
            N_SoilProfile = 1
            for iSoilname in Soilname
               # if soil changes
               if iSoilname ≠ Soilname_Initial
                  append!(iSoilProfile_Start, i)
                  append!(iSoilProfile_End, i-1)
                  push!(Soilname_SoilProfile, iSoilname)

                  Soilname_Initial = Soilname[i] # New soil
                  N_SoilProfile += 1
               elseif  i == NiZ
                  append!(iSoilProfile_End, i)  
               end  # if: name
               i += 1
            end # iSoilname
            
            return iSoilProfile_End, iSoilProfile_Start, N_SoilProfile, Soilname_SoilProfile
         end  # function: SMAP_SOIL_iSoilProfile


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : COMPUTE_θINI
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function COMPUTE_θINI(hydroSmap, iSiteName, Path_Output_θini, Path_SmapHydro, θᵢₙᵢ)
         # READING HYDRAULIC PARAMETERS
            # Deriving the number of soil SoilProfiles
            println(Path_SmapHydro)
               Data = DelimitedFiles.readdlm(Path_SmapHydro, ',')

               NiZ = size(Data)[1] - 1

         println(iSiteName," ====" ,θᵢₙᵢ)

         # COMPUTING θini
            θ₁ = fill(0.0::Float64, NiZ)

            θ₁[1] = max( min(θᵢₙᵢ, hydroSmap.θs[1]), hydroSmap.θr[1])

            Se = wrc.θ_2_Se(θ₁[1], 1, hydroSmap)

            # Assuming that all SoilProfiles have the same Se
            for iZ=2:NiZ
               θ₁[iZ] = wrc.Se_2_θ(Se, iZ, hydroSmap)
            end

         # Computing 1..NiZ for output file
            iZ = collect(1:1:NiZ)

         # Writing to file
            Header = ["iZ";"θini"; "SoilProfile"]

            Output = Tables.table( [iZ θ₁ iZ])
      
            CSV.write(Path_Output_θini, Output, header=Header)	   

      return θ₁
      end  # function: COMPUTE_θINI

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : TABLE_DISCRETISATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function TABLE_DISCRETISATION(Path, Z, iLayers)
         Header = ["iZ";"Z"; "Layer"]

         iZ = collect(1:1:length(Z))

         open(Path, "w") do io
				DelimitedFiles.writedlm(io,[Header] , ",") # Header
				DelimitedFiles.writedlm(io, [iZ Z iLayers], ",")
			end # open
      return nothing
      end  # function: TABLE


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_VEG
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TABLE_HYDRO_VEG(hydroSmap, N_iSoilProfiles, Path)
			println("			~ $(Path) ~")

			Id = 1:1:N_iSoilProfiles

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iSoilProfiles, hydroSmap)
					
			pushfirst!(FieldName_String, string("iSite")) # Write the "Id" at the very begenning

			open(Path, "w") do io
				DelimitedFiles.writedlm(io,[FieldName_String] , ",") # Header
				DelimitedFiles.writedlm(io, [Int64.(Id) Matrix], ",")
			end # open
      return nothing
		end  # function: TABLE_HYDRO

end # module smap2hypix
