# =============================================================
#		MODULE: table
# =============================================================
module table
	import Tables, CSV
	export TABLE_ID

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TABLE_EXTRAPOINTS_K
	# 		Tabular values of the PSD model
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TABLE_ID(N_iZ::Int64, Path::String)
			println("    ~  $(Path) ~")

			IdSelect = collect(1:1:N_iZ)

			Select = fill(1::Int64, N_iZ)

			FieldName_String = ["Id", path.option.Select]

			Output = Tables.table( [IdSelect[1:N_iZ] Select[1:N_iZ]] )
			
			CSV.write(Path, Output, header=FieldName_String, delim=',')
		return nothing
		end  # function:  TABLE_ID


	# =============================================================
	#		MODULE: hydroLab
	# =============================================================
	module hydroLab
		import  ...tool, ...wrc, ...kunsat
		import DelimitedFiles, Tables, CSV
		export θΨK

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK(hydro, hydroOther, IdSelect, Kₛ_Model, N_iZ::Int64, Path)
				println("    ~  $(Path) ~")

				Matrix₁, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ, hydro)

				Matrix₂, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ, hydroOther)

				Header = vcat("Id", FieldName_String, FieldName_String2, "KsModel")

				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[Header] , ",",) # Header
					DelimitedFiles.writedlm(io, [IdSelect Matrix₁ Matrix₂ Kₛ_Model], ",")
				end
			return nothing
			end  # function:  θΨK


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TABLE_θΨK
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TABLE_EXTRAPOINTS_θΨ(optionₘ, hydro, IdSelect, N_iZ::Int64, Path, Ψ_Table; Orientation="Horizontal")
				println("    ~  $(Path) ~")

				N_Ψ = Int64(length(Ψ_Table))

				if Orientation == "Horizontal" # <>=<>=<>=<>=<>=<>
					# Writting the Header
						FieldName_String = fill(""::String, N_Ψ)
						for i =1:N_Ψ
							FieldName_String[i] = string(Int64(Ψ_Table[i]) ) * "mm"
						end
						pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
							
					# Computing θ at required θ
						θ₂ = fill(0.0::Float64, (N_iZ, N_Ψ))

						for iZ=1:N_iZ
							for iΨ =1:N_Ψ
								Ψ₂ = Ψ_Table[iΨ]
								θ₂[iZ, iΨ] = wrc. Ψ_2_θDual(optionₘ, Ψ₂, iZ, hydro)
							end # iΨ
						end # iZ

						open(Path, "w") do io
							DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
							DelimitedFiles.writedlm(io, [Int64.(IdSelect) θ₂], ",")
						end

				elseif Orientation == "Vertical" # <>=<>=<>=<>=<>=<>
					FieldName_String = ["Id","H[mm]","Theta[-]"]
					N = N_Ψ * N_iZ
					Id₂ = Vector{Int64}(undef, N)
					Ψ₂  = Vector{Float64}(undef,  N)
					θ₂  = Vector{Float64}(undef, N)
					iCount = 1

					for iZ=1:N_iZ
						for iΨ =1:N_Ψ
							Id₂[iCount] = IdSelect[iZ]
							Ψ₂[iCount] = Ψ_Table[iΨ]
							θ₂[iCount] = wrc. Ψ_2_θDual(optionₘ,Ψ₂[iCount], iZ, hydro)
	
							iCount+=1
						end # iΨ
					end # iZ

					open(Path, "w") do io
						DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
						DelimitedFiles.writedlm(io, [Id₂[1:N] Ψ₂[1:N] θ₂[1:N]], ",")
					end
				else
					error("SoilWaterToolBox Error in TABLE_EXTRAPOINTS_θΨ $Orientation not understood")
				end
		return nothing
		end  # function:  θΨ


		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TABLE_EXTRAPOINTS_K
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TABLE_EXTRAPOINTS_Kθ(optionₘ, hydroParam, IdSelect, K_Table, Kₛ_Model, N_iZ::Int64, Path::String)
				println("    ~  $(Path) ~")

				N_K = Int64(length(K_Table))

			# Writting the Header
				FieldName_String =["Id", "H[mm]" ,"Kunsat[mm_s]"]
							
			# Computing K at required Ψ
				N = N_K *N_iZ
				Id₂     = Vector{Int64}(undef, N)
				Ψ₂      = Vector{Float64}(undef,  N)
				Kunsat₂ = Vector{Float64}(undef, N)
				iCount  = 1
				hydroParam₂ = deepcopy(hydroParam)
				for iZ=1:N_iZ
					for iK =1:N_K
						hydroParam₂.Ks[iZ] = Kₛ_Model[iZ]

						Id₂[iCount] = IdSelect[iZ]
						Ψ₂[iCount] = K_Table[iK]
						Kunsat₂[iCount] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ₂[iCount], iZ, hydroParam₂)

						iCount+=1
					end # iΨ
				end # iZ

				Output = Tables.table([string.(Id₂[1:N]) Ψ₂[1:N] Kunsat₂[1:N]])
				CSV.write(Path, Output, header=FieldName_String, delim=',')
				
		return nothing
		end  # function:  θΨ
			
	end  # module hydro
	# ............................................................


	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		import...tool, ...wrc, ...cst
		import DelimitedFiles
		export PSD, θΨK_PSD

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD(IdSelect, N_iZ::Int64, paramPsd, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ,  paramPsd)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Int64.(IdSelect) round.(Matrix,digits=5)], ",")
				end
			return nothing
			end


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK_PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK_PSD(hydroPsd, IdSelect, KunsatModel_Psd, N_iZ::Int64, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ, hydroPsd)

				Matrix = hcat(Matrix, KunsatModel_Psd)

				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				push!(FieldName_String, "Kunsat_Model")

				Matrix =  round.(Matrix, digits=5)
				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [IdSelect Matrix], ",")
				end
			return nothing
			end  # function:  θΨK_PSD



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD_θΨ_θ
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD_θΨ_θ(hydroPsd, IdSelect, N_iZ::Int64, param, Path::String)
				println("    ~  $Path ~")

				N_Ψ = Int64(length(param.psd.Ψ_Table))

				# Writting the Header
					FieldName_String = Array{String}(undef, (N_Ψ))

					for i =1:N_Ψ
						FieldName_String[i] = string(param.psd.Ψ_Table[i] * cst.Mm_2_Cm) * "cm"
					end
					pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				
				# Computing θ at required θ
					θ = Array{Float64}(undef, (N_iZ, N_Ψ))
					for iZ=1:N_iZ
						for iΨ =1:N_Ψ
							Ψ = param.psd.Ψ_Table[iΨ]
							θ[iZ, iΨ] = wrc. Ψ_2_θDual(optionₘ,Ψ, iZ, hydroPsd)
						end # iΨ
					end # iZ

				# Concatenating the 2 matrices
					θ = hcat(IdSelect, θ)
					θ = round.(θ, digits=5)

				# Writting the table
					open(Path, "w") do io
						DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
						for i = 1:length(IdSelect)
							DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
							DelimitedFiles.writedlm(io, [θ[i, 1:N_Ψ+1]], ",")
						end # i
					end # Path
			return nothing
			end  # function:  θΨK_PSD

	end  # module psd
	# ............................................................
	
	

	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		import ...tool
		import DelimitedFiles
		export HYDRO_INFILT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO_INFILT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO_INFILT(hydroInfilt, IdSelect, KunsatModel_Infilt, N_iZ::Int64, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ, hydroInfilt)

				Matrix = hcat(Matrix, KunsatModel_Infilt)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				push!(FieldName_String, string("Kunsat_θΨ"))

				Matrix =  round.(Matrix, digits=10)
				open(Path, "w") do io
					DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [IdSelect Matrix], ",")
				end
			return nothing
			end  # function: HYDRO_INFILT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : infilt
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function INFILT(IdSelect, N_iZ, infiltOutput, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iZ::Int64, infiltOutput)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				Matrix =  round.(Matrix, digits=5)
				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [IdSelect Matrix], ",")
				end
			return nothing
			end  # function: HYDRO_INFILT
		
	end  # module: infilt

	# .........................................................................................


	# =============================================================
	#		module: tableHypix
	# =============================================================
	module hyPix
		import ...cst, ...tool, ...wrc, ...kunsat
		import DelimitedFiles
		import Dates: value, DateTime, year, month, day, hour, minute, second
		export DAILY_CLIMATE, DISCRETIZATION, HYDRO, PERFORMANCE, Q, TIME_SERIES, TIME_SERIES_DAILY, VEG, θ, θΨ, Ψ, θAVERAGE

		# ===================================================
		#          Discretization
		# ===================================================
			function DISCRETIZATION(discret, N_iZ, Z, pathHyPix)
				println("			~  $(pathHyPix.Table_Discretisation) ~")

				Header =  ["Z" "ΔZ" "ΔZ_⬓" "Znode" "ΔZ_Aver" "ΔZ_W" "Z_CellUp"]

				open(pathHyPix.Table_Discretisation, "w") do io
					DelimitedFiles.writedlm(io,[Header] , ",",) # Header
					DelimitedFiles.writedlm(io, [Z[1:N_iZ] discret.ΔZ[1:N_iZ] discret.ΔZ_⬓[1:N_iZ] discret.Znode[1:N_iZ] discret.ΔZ_Aver[1:N_iZ] discret.ΔZ_W[1:N_iZ] discret.Z_CellUp[1:N_iZ]], ",")
				end
			return nothing
			end # Table DISCRETIZATION


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO(hydroHorizon, iSim, N_iHorizon, pathHyPix)
				Path = pathHyPix.Table_Hydro * "_" * string(iSim) * ".csv"
				println("			~ $(Path) ~")

				Id = 1:1:N_iHorizon

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iHorizon, hydroHorizon)
						
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Int64.(Id) Matrix], ",")
				end
			return nothing			
			end  # function: HYDRO

			
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : veg
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function VEG(veg, iSim, pathHyPix)
				Path = pathHyPix.Table_Veg * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(1, veg)

				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Matrix], ",")
				end
			return nothing
			end  # function: VEG


		# ===================================================
		#          TimeStep at ΔT
		# ===================================================
			function TIME_SERIES(∑T, ΔT, ∑Pr, ΔPr, Hpond, Recharge, ∑WaterBalance_η, iSim, pathHyPix)		
				Path = pathHyPix.Table_TimeSerie * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")
				
				Header =  ["∑T[mm]" "ΔT[mm]" "∑Pr[mm/ΔT]" "ΔPr[mm/ΔT]" "Hpond[mm]" "Recharge[mm/ΔT]" "∑WaterBalance_η[mm]"]

				open(Path, "w") do io
					DelimitedFiles.writedlm(io,[Header] , ",",) # Header
					DelimitedFiles.writedlm(io, [∑T ΔT ∑Pr ΔPr Hpond Recharge ∑WaterBalance_η], ",")
				end
			return nothing
			end # Table DISCRETIZATION


		# ===================================================
		#          TimeStep daily
		# ===================================================
			function TIME_SERIES_DAILY(∑T_Plot, ∑WaterBalance_η_Plot, Date_Plot, iSim, N_∑T_Plot, ΔEvaporation_Plot, ΔRecharge_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, pathHyPix)
				Header =  ["iD" "Year" "Month" "Day" "Hour" "Minute" "Second" "∑T[Hour]" "ΔPr_Through[mm/day]" "ΔPet[mm/day]" "ΔSink[mm/day]" "ΔEvaporation[mm/day]" "Hpond≈[mm]" "Recharge[mm/day]" "∑WaterBalance_η_Profile[mm/day]"]

				Path = pathHyPix.Table_TimeSerie_Daily * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")

				Id = 1:1:N_∑T_Plot

				Year₁   =fill(0::Int64, N_∑T_Plot)
				Month₁  =fill(0::Int64, N_∑T_Plot)
				Day₁    =fill(0::Int64, N_∑T_Plot)
				Hour₁   =fill(0::Int64, N_∑T_Plot)
				Minute₁ =fill(0::Int64, N_∑T_Plot)
				Second₁ =fill(0::Int64, N_∑T_Plot)

				for iT=1:N_∑T_Plot
					Year₁[iT] = year(Date_Plot[iT])
					Month₁[iT] = month(Date_Plot[iT])
					Day₁[iT] = day(Date_Plot[iT])
					Hour₁[iT] = hour(Date_Plot[iT])
					Minute₁[iT] = minute(Date_Plot[iT]) 
					Second₁[iT] = second(Date_Plot[iT])
				end

				open(Path, "w") do io
					DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[Header] , ",",) # Header
					DelimitedFiles.writedlm(io, [Id Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ∑T_Plot 10.0.*ΔPr_Plot ΔPet_Plot ΔSink_Plot ΔEvaporation_Plot 10.0.*ΔPond_Plot ΔRecharge_Plot ∑WaterBalance_η_Plot], ",")
				end
			return nothing
			end # Table  TIME_SERIES_DAILY


		# ===================================================
		#          θ
		# ===================================================
			function θ(∑T, θ, Znode, iSim, pathHyPix)
				Path = pathHyPix.Table_θ * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")

				# Adding an other column
				prepend!(Znode, -999)

				DelimitedFiles.writedlm(Path, [transpose(Znode); ∑T θ], ",")
			return nothing
			end  # Table θ


		# ===================================================
		#          Q
		# ===================================================
			function Q(∑T, Q, Z_Bottom, Znode, iSim, pathHyPix)	
				Path = pathHyPix.Table_Q * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")
				
				# Adding an other column
				prepend!(Znode, -999)
				append!(Znode, Z_Bottom)

				DelimitedFiles.writedlm(Path, [transpose(Znode); ∑T Q], ",")
			return nothing
			end  # function Q


		# ===================================================
		#          Ψ
		# ===================================================
			function Ψ(∑T, Ψ, Znode, iSim, pathHyPix)
				Path = pathHyPix.Table_Ψ * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")

				# Adding an other column
				prepend!(Znode, -999)

				DelimitedFiles.writedlm(Path, [transpose(Znode); ∑T Ψ], ",")
			return nothing
			end  # function Ψ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨ
		# 		Tabular values of the hydroParam model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨ(hydroHorizon, iSim, N_iHorizon, optionₘ, param, pathHyPix)		
				Path = pathHyPix.Table_θΨ * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")

				N_θΨobs = Int64(length(param.hyPix.ploting.θΨ_Table))

				# Writting the Header
					FieldName_String = fill(""::String, N_θΨobs)

					for i =1:N_θΨobs
						FieldName_String[i] = string(param.hyPix.ploting.θΨ_Table[i] * cst.Mm_2_Cm) * "cm"
					end
					pushfirst!(FieldName_String, string("Layer")) # Write the "Id" at the very begenning
				
				# Computing θ at required θ
					θ_Mod = fill(0.0::Float64, (N_iHorizon, N_θΨobs))
					for iZ=1:N_iHorizon, iΨ =1:N_θΨobs
							Ψ_Mod =param.hyPix.ploting.θΨ_Table[iΨ]
							θ_Mod[iZ, iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψ_Mod, iZ, hydroHorizon)
					end # iZ

				# Concatenating the 2 matrices
				Id = 1:1:N_iHorizon

				θ_Mod = hcat(Id, θ_Mod)

				# Writting the table
					open(Path, "w") do io
						DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
						for iZ = 1:N_iHorizon
							# DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
							DelimitedFiles.writedlm(io, [θ_Mod[iZ, 1:N_θΨobs+1]], ",")
						end # i
					end # Path
			return nothing	
			end  # function:  θΨK_PSD


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KΨ
		# 		Tabular values of the hydroParam model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KΨ(hydroHorizon, iSim, N_iHorizon, optionₘ, param, pathHyPix)				
				Path = pathHyPix.Table_KΨ * "_" * string(iSim) * ".csv"
				println("			~  $(Path) ~")

				N_θΨobs = Int64(length(param.hyPix.ploting.θΨ_Table))

				# Writting the Header
					FieldName_String = fill(""::String, N_θΨobs)

					for i =1:N_θΨobs
						FieldName_String[i] = string(param.hyPix.ploting.θΨ_Table[i] * cst.Mm_2_Cm) * "cm"
					end
					pushfirst!(FieldName_String, string("Layer Cm/H")) # Write the "Id" at the very begenning
				
				# Computing θ at required θ
					K_Mod = fill(0.0::Float64, (N_iHorizon, N_θΨobs))
					for iZ=1:N_iHorizon, iΨ =1:N_θΨobs
							Ψ_Mod =param.hyPix.ploting.θΨ_Table[iΨ]
							K_Mod[iZ, iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_Mod, iZ, hydroHorizon) .* cst.MmS_2_CmH
					end # iZ

				# Concatenating the 2 matrices
				Id = 1:1:N_iHorizon

				K_Mod = hcat(Id, K_Mod)

				# Writting the table
					open(Path, "w") do io
						DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
						for iZ = 1:N_iHorizon
							DelimitedFiles.writedlm(io, [K_Mod[iZ,1:N_θΨobs+1]], ",")
						end # i
					end # Path
			return nothing	
			end  # function:  θΨK_PSD


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PERFORMACE
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iSim, iSim, RmseBest, SwcRoots, WofBest, ΔRunTimeHypix, ΔT_Average, SiteName_Hypix, pathHyPix)
				Path = pathHyPix.Table_Performance * ".csv"
				println("			~  $(Path) ~")

				Header = ["Id" "SiteName" "WofBest" "NseBest" "Efficiency" "Global_WaterBalance" "Global_WaterBalance_NormPr" "ΔT_Average" "∑∑ΔSink" "∑ΔQ_Bot" "SwcRoots" "iNonConverge" "ΔRunTimeHypix"]

				Id = 1:1:length(WofBest)

				open(Path, "w") do io
					# DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[Header] , ",",) # Header
					DelimitedFiles.writedlm(io, [Id SiteName_Hypix WofBest RmseBest Efficiency Global_WaterBalance Global_WaterBalance_NormPr ΔT_Average ∑∑ΔSink ∑ΔQ_Bot SwcRoots iNonConverge_iSim ΔRunTimeHypix], ",")
				end
			return nothing
			end # function PERFORMACE


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : DAILY_CLIMATE
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function DAILY_CLIMATE(∑T_Climate, clim, iSim, pathHyPix)
					Path = pathHyPix.Table_DailyClimate * "_" * string(iSim) * ".csv"
					println("			~  $(Path) ~")

					local ∑T_Int = ceil.(Int, ∑T_Climate[1:clim.N_Climate] .* cst.Second_2_Day)

					Header = ["Year" "Month" "Day" "Pr" "Pr_Ground"]

					open(Path, "w") do io
						DelimitedFiles.writedlm(io,[Header] , ",",) # Header
						DelimitedFiles.writedlm(io, [year.(clim.Date[1:clim.N_Climate]) month.(clim.Date[1:clim.N_Climate]) day.(clim.Date[1:clim.N_Climate]) clim.Pr[1:clim.N_Climate] clim.Pr_Through] , ",")
					end
				return nothing	
				end  # function: DAILY_CLIMATE


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : θAVERAGE
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function θAVERAGE(Date_Plot, iSim, θobs_Plot, θsim_Aver, pathHyPix)
					Path = pathHyPix.Table_θaverage * ".csv"
					println("			~  $(Path) ~")

					Header = ["Id", "Year","Month","Day" ,"θobs_Aver", "θsim_Aver"]

					Id = 1:1:length(θsim_Aver)

					Year = year.(Date_Plot)
					Month = month.(Date_Plot)
					Day = day.(Date_Plot)

					open(Path, "w") do io
						DelimitedFiles.writedlm(io,[Header] , ",",) # Header
						DelimitedFiles.writedlm(io, [Id Year Month Day θobs_Plot θsim_Aver] , ",")
					end # open
				return nothing			
				end # function: θAVERAGE

	end  # module hyPix ............................................................

end  # module table
# ............................................................