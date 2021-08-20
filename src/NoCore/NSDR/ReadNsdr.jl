# =============================================================
#		module: readNsd
# =============================================================


module readNsdr
	import DelimitedFiles


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : READ_NSDR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function READ_NSDR()

		Path_Input = "D:\\Main\\MODELS\\SoilWater_ToolBox\\data\\INPUT\\Data_SoilWater\\Nsdr\\Nsdr_ThetaH_Raw.csv"

		# Read data
			Data = DelimitedFiles.readdlm(Path_Input, ',')
		# Read header
			Header = Data[1,1:end]
		# Remove first READ_ROW_SELECT
			Data = Data[2:end,begin:end]
		# Reading
			Date, N = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")


		
		return
	end  # function: READ_NSDR
	
end  # module: readNsd
# ............................................................

