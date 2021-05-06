# =============================================================
#		module: jules
# =============================================================
module jules
   import ..path, ..option, ..tool
   import DelimitedFiles, Dates

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : START_JULES
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function START_JULES()
      println("    ~  $(path.JulesMetadata) ~")

      Path_Climate = "D:\\DATAraw\\JULESdata\\Climate\\VCSN_Obs\\"

      # Read data
         Data = DelimitedFiles.readdlm(path.JulesMetadata, ',')
      # Read header
         Header = Data[1,1:end]
      # Remove first READ_ROW_SELECT
         Data = Data[2:end,begin:end]
      # Reading
         SiteName, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SiteName")
         SiteNumber, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "VCSNgridnumber")

      # Dictionary
      # SiteName2SiteNumber::Dict{String, Int64} 

         i = 1
         SiteName2SiteNumber = Dict("a"=>9999) # Initializing
         for iSiteName in SiteName
            Path_Output =  path.Home * "//INPUT//DataHyPix//JULES//JulesInput//" * iSiteName
            mkpath(Path_Output) 
            SiteName2SiteNumber[iSiteName] = SiteNumber[i]
            i += 1
            Path_Climate_Input = Path_Climate * "VCSN_obsSM_" * string(SiteName2SiteNumber[iSiteName]) * ".csv"

            READ_CLIMATE(Path_Climate_Input)

            Path_Climate_Output = Path_Output * "//VCSN_obsSM_" * string(SiteName2SiteNumber[iSiteName]) * ".csv"

            # Copy paste
            # cp(Path_Climate_Input, Path_Climate_Output, force=true )
         
            # println(Path_Climate_Input)
         end

      
      return
   end  # function: START_JULES

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : READ_CLIMATE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function READ_CLIMATE(Path)

       # Read data
          Data = DelimitedFiles.readdlm(Path, ',')
       # Read header
          Header = Data[1,1:end]
       # Remove first READ_ROW_SELECT
          Data = Data[2:end,begin:end]
       # Reading
         Date, N = tool.readWrite.READ_HEADER_FAST(Data, Header, "OBS_DATE")
         Rain, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "rain_fill")
         Pet, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "pet_fao56")

         Year  = Dates.year.(Dates.DateTime.(Date))
         Month = Dates.month.(Dates.DateTime.(Date))
         Day   = Dates.day.(Dates.DateTime.(Date))

         # for iDate in Date
         #    Year = Dates.year(Dates.DateTime(iDate))
         #    Month = Dates.month(Dates.DateTime(iDate))
         #    Day = Dates.day(Dates.DateTime(iDate))

         #    println(iDate,",",Year, ",", Month,",", Day)
         # end
        
        @show Day
      
      return
   end  # function: READ_CLIMATE
   
end  # module: jules
# ............................................................