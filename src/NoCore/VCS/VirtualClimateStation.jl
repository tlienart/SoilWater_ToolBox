module virtualClimateStation

   PathInputVcsAgentList = r"D:\DATAraw\CLIMATE\VCS\AGENTlist\Agent_Lat_Longt.csv"
   PathInputCliflo =r"D:\DATAraw\CLIMATE\CliFlo\CliFlo precipitation and PET sites.csv"
   PathOutput = r"D:\Main\MODELS\SoilWater_ToolBox\data\OUTPUT\VirtialClimateStation"

   import CSV, Tables

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : READ_VCS
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function READ_VCS( PathInputVcsAgentList, PathInputCliflo)

      # ~~~~~~ READING PathInputCliflo ~~~~~~~
         Data = CSV.File(PathInputCliflo, header=true)
         Header = string.(Tables.columnnames(Data))

         Agent_Cliflo = convert(Vector{Float64}, Tables.getcolumn(Data, :Agent))
         Lat_Cliflo = convert(Vector{Float64}, Tables.getcolumn(Data, :Lat))
         Long_Cliflo =convert(Vector{Float64}, Tables.getcolumn(Data, :Long))
      #-------------------------------------------------------------------------


      # ~~~~~~ READING PathInputVcsAgentList ~~~~~~~
      Data = CSV.File(PathInputVcsAgentList, header=true)
      Header = string.(Tables.columnnames(Data))

      Agent_Vcs = convert(Vector{Float64}, Tables.getcolumn(Data, :AGENT))
      Lat_Vcs = convert(Vector{Float64}, Tables.getcolumn(Data, :LATITUDE))
      Long_Vcs =convert(Vector{Float64}, Tables.getcolumn(Data, :LONGTITUDE))
   #-------------------------------------------------------------------------


   return 
   end
   #-------------------------------------------------------------------------

end