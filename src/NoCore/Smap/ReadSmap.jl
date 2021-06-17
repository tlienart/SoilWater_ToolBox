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
               Id₂, N_iZ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")

               Soilname₂, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")

               Ψdata = []
               θData = []
               for iHeader in Header
                  if occursin("wrc", iHeader)
                     θ₀, N_iZ = tool.readWrite.READ_HEADER_FAST(Data, Header, iHeader)

                     iHeader = replace(iHeader, "wrc" => "")
                     iHeader = replace(iHeader, "kpa" => "")
                     iHeader = replace(iHeader, " " => "")
                     iHeader_Float=  parse(Float64, iHeader)

                     iHeader_Float = iHeader_Float * cst.kPa_2_Mm

                     append!(Ψdata, iHeader_Float)

                     try
                        θData = hcat(θData[1:N_iZ, :], θ₀[1:N_iZ])
                     catch
                        θData = θ₀[1:N_iZ]
                     end
                  end # occursin("wrc", iHeader)
               end # for iHeader in Header

               θ_θΨobs₂ = zeros(Float64, N_iZ, length(Ψdata))
               Ψ_θΨobs₂ = zeros(Float64, N_iZ, length(Ψdata))
               N_θΨobs₂ = zeros(Int64, N_iZ)
    
               for iZ=1:N_iZ
                  iΨ_Count = 1
                  for iΨ=1:length(Ψdata)
                     if !isnan(θData[iZ, iΨ])
                        Ψ_θΨobs₂[iZ, iΨ_Count] = Ψdata[iΨ]
                        θ_θΨobs₂[iZ, iΨ_Count] = θData[iZ, iΨ]
                        N_θΨobs₂[iZ] += 1
                        iΨ_Count += 1
                     end #  !isnan(θData[iZ, iΨ])
                  end # iΨ
               end # iZ

         return Id₂, N_θΨobs₂, Soilname₂, θ_θΨobs₂, Ψ_θΨobs₂
         end  # function: DATA2D
   
end  # module: readSmap
# ............................................................