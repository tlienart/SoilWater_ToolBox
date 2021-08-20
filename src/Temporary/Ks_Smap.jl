

module temporary
	import ..tool
	import DelimitedFiles 
	using CairoMakie, ColorSchemes, Colors

	CairoMakie.activate!()

	function KS_SMAP()
		PathInput = "D:\\DATAraw\\SoilWater_ToolBox\\Ksat\\Smap Hydrological_Ksat20210624_Cleaned.csv"
		PathOutput = "D:\\Main\\MODELS\\SoilWater_ToolBox\\data\\OUTPUT\\Temporary\\PlotSmap.svg"

		println(" === KS_SMAP ====")

		# READING
			Data = DelimitedFiles.readdlm(PathInput, ',')
			Header = Data[1,1:end]
			Data = Data[2:end,begin:end]

			FHPermeabilityCode, N = tool.readWrite.READ_HEADER_FAST(Data, Header,"FHPermeabilityCode")
			KsModel, ~            = tool.readWrite.READ_HEADER_FAST(Data, Header, "Ks_Kg")
			KsTrevor, ~           = tool.readWrite.READ_HEADER_FAST(Data, Header, "FHPermeabilityBestValue_mm")
			SmapFH, ~             = tool.readWrite.READ_HEADER_FAST(Data, Header, "SmapFH")

		# SMALL CONVERSION
			KsModel_Log10 = fill(0.0::Float64, N)
			KsTrevor_Log10 = fill(0.0::Float64, N)
			for i=1:N
				KsModel_Log10[i] = (KsModel[i] * 60.0)  # Convert tmm/s -> mm/h
				KsTrevor_Log10[i] = log10(KsTrevor[i]) # mm/h
			end

		# UNIQUE FHPermeabilityCode
			Categorie_Unique = unique(FHPermeabilityCode)
			N_Categories = length(Categorie_Unique)

		# DERIVING CLASSES
			KsModel_Class = fill(0.0::Float64, (10000, N_Categories))
			KsTrevor_Class = fill(0.0::Float64, (10000, N_Categories))
			iN_Class = fill(1, length(Categorie_Unique))
			i = 1		
			for i = 1:N		
				iPos = findfirst(x -> x==FHPermeabilityCode[i], Categorie_Unique)

				KsModel_Class[iN_Class[iPos], iPos] = KsModel_Log10[i]
				KsTrevor_Class[iN_Class[iPos], iPos] = KsTrevor_Log10[i] 
				iN_Class[iPos] = iN_Class[iPos] + 1 
			end # for i = 1:N	

			# Small correction
			for iPos = 1:N_Categories
				iN_Class[iPos] = max(iN_Class[iPos] - 1, 1)
			end

		# PLOTTING	
			colors = ColorScheme(range(colorant"black", colorant"red", length=N_Categories))
		
			Fig = CairoMakie.Figure(resolution = (700, 800))
			Axis1 = CairoMakie.Axis(Fig, xlabel = "variable", yscale=log10, ylabel = "Ks [mm / h]", xticks = (1:N_Categories, Categorie_Unique))
			Axis2 = CairoMakie.Axis(Fig, xlabel = "variable", yscale=log10, ylabel = "Ks [mm / h]", xticks = (1:N_Categories, Categorie_Unique))

			for iPos=1:N_Categories
				X = fill(iPos, iN_Class[iPos])
				Y = KsModel_Class[1:iN_Class[iPos], iPos]

				boxplot!(Axis1, X, Y; whiskerwidth = 1, width = 0.35, color=(colors[iPos], 0.45), whiskercolor = (colors[iPos], 1), mediancolor = :black)
				 
				# violin!(Axis2, X, Y; width = 0.35, color=(colors[iPos], 0.45), strokecolor = colors[iPos], show_median = true, mediancolor = :black)
			end

		Fig[1,1] = Axis1
		Fig[2,1] = Axis2
		CairoMakie.save(PathOutput, Fig)
	end # KS_SMAP


	# function CLASSES(Categorie_Code, Categorie_Unique, N)
	# 	iN_Class = fill(1, length(Categorie_Unique))
	# 	for i = 1:N		
	# 		iPos = findfirst(x -> x==Categorie_Code[i], Categorie_Unique)

	# 		KsModel_Class[iN_Class[iPos], iPos] = KsModel_Log10[i]

	# 		iN_Class[iPos] = iN_Class[iPos] + 1 
	# 	end # for i = 1:N	

	# 	# Small correction
	# 	for iPos = 1:length(iN_Class)
	# 		iN_Class[iPos] = max(iN_Class[iPos] - 1, 1)
	# 	end
		
	# end

end # module
