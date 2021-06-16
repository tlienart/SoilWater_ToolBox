# =============================================================
#		module: readSmap
# =============================================================
module readSmap
   import ..tool, ..cst

   using DelimitedFiles
   export DATA2D, SMAP, ROCKFRAGMENT_WETTABLE_STRUCT

	
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #		FUNCTION : DATA2D
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function DATA2D(Path)
            println("    ~  $(Path) ~")

            # Read data
				   Data = DelimitedFiles.readdlm(Path, ',')
            # Read header
               Header = Data[1,1:end]
            # Remove first READ_ROW_SELECT
               Data = Data[2:end,begin:end]
            # Sort data
               Data = sortslices(Data, dims=1)

            # Read data of interest
               Id₂, N_SoilSelect = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")

               Soilname₂, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")

               Ψdata = []
               θData = []
               for iHeader in Header
                  if occursin("wrc", iHeader)
                     θ₀, N_SoilSelect = tool.readWrite.READ_HEADER_FAST(Data, Header, iHeader)

                     iHeader = replace(iHeader, "wrc" => "")
                     iHeader = replace(iHeader, "kpa" => "")
                     iHeader = replace(iHeader, " " => "")
                     iHeader_Float=  parse(Float64, iHeader)

                     iHeader_Float = iHeader_Float * cst.kPa_2_Mm

                     append!(Ψdata, iHeader_Float)

                     try
                        θData = hcat(θData[1:N_SoilSelect, :], θ₀[1:N_SoilSelect])
                     catch
                        θData = θ₀[1:N_SoilSelect]
                     end
                  end # occursin("wrc", iHeader)
               end # for iHeader in Header

               θ_θΨ₂ = zeros(Float64, N_SoilSelect, length(Ψdata))
               Ψ_θΨ₂ = zeros(Float64, N_SoilSelect, length(Ψdata))
               N_θΨ₂ = zeros(Int64, N_SoilSelect)
    
               for iZ=1:N_SoilSelect
                  iΨ_Count = 1
                  for iΨ=1:length(Ψdata)
                     if !isnan(θData[iZ, iΨ])
                        Ψ_θΨ₂[iZ, iΨ_Count] = Ψdata[iΨ]
                        θ_θΨ₂[iZ, iΨ_Count] = θData[iZ, iΨ]
                        N_θΨ₂[iZ] += 1
                        iΨ_Count += 1
                     end #  !isnan(θData[iZ, iΨ])
                  end # iΨ
               end # iZ

         return Id₂, N_θΨ₂, Soilname₂, θ_θΨ₂, Ψ_θΨ₂
         end  # function: DATA2D
   
end  # module: readSmap
# ............................................................