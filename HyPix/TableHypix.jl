# =============================================================
#		module: tableHypix
# =============================================================
module tableHypix

	import ..cst, ..param, ..path, ..tool, ..wrc, ..kunsat
	import DelimitedFiles
	import Dates: value, DateTime, year, month, day, hour, minute, second
	
	export DAILY_CLIMATE, DISCRETIZATION, HYDRO, PERFORMANCE, Q, SIGNATURE, TIME_SERIES, TIME_SERIES_DAILY, VEG, θ, θΨ, Ψ

	# ===================================================
	#          Discretization
	# ===================================================
		function DISCRETIZATION(discret, N_iZ, Z)
			println("			~  $(path.Table_Discretisation) ~")

			Header =  ["Z" "ΔZ" "ΔZ_⬓" "Znode" "ΔZ_Aver" "ΔZ_W" "Z_CellUp"]

			open(path.Table_Discretisation, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				DelimitedFiles.writedlm(io, [Z[1:N_iZ] discret.ΔZ[1:N_iZ] discret.ΔZ_⬓[1:N_iZ] discret.Znode[1:N_iZ] discret.ΔZ_Aver[1:N_iZ] discret.ΔZ_W[1:N_iZ] discret.Z_CellUp[1:N_iZ]], ",")
			end
		end # Table DISCRETIZATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDRO(hydroHorizon, iSim, N_iHorizon)
			Path = path.Table_Hydro * "_" * string(iSim) * ".csv"
			println("			~ $(Path) ~")

			Id = 1:1:N_iHorizon

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iHorizon, hydroHorizon)
					
			pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				DelimitedFiles.writedlm(io, [Int64.(Id) Matrix], ",")
			end
		end  # function: HYDRO

		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : veg
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function VEG(veg, iSim)

			Path = path.Table_Veg * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME_PARAM(veg)

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				DelimitedFiles.writedlm(io, [Matrix], ",")
			end
			
		end  # function: VEG


	# ===================================================
	#          TimeStep at ΔT
	# ===================================================
		function TIME_SERIES(∑T, ΔT, ∑Pr, ΔPr, Hpond, Recharge, ∑WaterBalance_η, iSim)
			Header =  ["∑T[mm]" "ΔT[mm]" "∑Pr[mm/ΔT]" "ΔPr[mm/ΔT]" "Hpond[mm]" "Recharge[mm/ΔT]" "∑WaterBalance_η[mm]"]

			Path = path.Table_TimeSerie * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				DelimitedFiles.writedlm(io, [∑T ΔT ∑Pr ΔPr Hpond Recharge ∑WaterBalance_η], ",")
			end
		end # Table DISCRETIZATION


	# ===================================================
	#          TimeStep daily
	# ===================================================
		function TIME_SERIES_DAILY(∑T_Plot, ∑WaterBalance_η_Plot, Date_Plot, iSim, N_∑T_Plot, ΔEvaporation_Plot, ΔRecharge_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot)
			Header =  ["iD" "Year" "Month" "Day" "Hour" "Minute" "Second" "∑T[Hour]" "ΔPr_Through[mm/day]" "ΔPet[mm/day]" "ΔSink[mm/day]" "ΔEvaporation[mm/day]" "Hpond≈[mm]" "Recharge[mm/day]" "∑WaterBalance_η_Profile[mm/day]"]

			Path = path.Table_TimeSerie_Daily * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			Id = 1:1:N_∑T_Plot

			Year₁   = Vector{Int64}(undef, N_∑T_Plot)
			Month₁  = Vector{Int64}(undef, N_∑T_Plot)
			Day₁    = Vector{Int64}(undef, N_∑T_Plot)
			Hour₁   = Vector{Int64}(undef, N_∑T_Plot)
			Minute₁ = Vector{Int64}(undef, N_∑T_Plot)
			Second₁ = Vector{Int64}(undef, N_∑T_Plot)

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
		end # Table DISCRETIZATION


	# ===================================================
	#          θ
	# ===================================================
		function θ(∑T, θ, Znode, iSim)
			Path = path.Table_θ * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			# Adding an other column
			prepend!(Znode, -999)

			DelimitedFiles.writedlm(Path, [transpose(Znode); ∑T θ], ",")
		end  # function θ


	# ===================================================
	#          Q
	# ===================================================
		function Q(∑T, Q, Z_Bottom, Znode, iSim)	
			Path = path.Table_Q * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")
			
			# Adding an other column
			prepend!(Znode, -999)
			append!(Znode, Z_Bottom)

			DelimitedFiles.writedlm(Path, [transpose(Znode); ∑T Q], ",")
		end  # function Q


	# ===================================================
	#          Ψ
	# ===================================================
		function Ψ(∑T, Ψ, Znode, iSim)

			Path = path.Table_Ψ * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			# Adding an other column
			prepend!(Znode, -999)

			DelimitedFiles.writedlm(Path, [transpose(Znode); ∑T Ψ], ",")
		end  # function Ψ



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θΨ
	# 		Tabular values of the hydroParam model
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨ(hydroHorizon, iSim, N_iHorizon)
			
			Path = path.Table_θΨ * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			N_θΨ = Int64(length(param.hypix.plot.θΨ_Table))

			# Writting the Header
				FieldName_String = Array{String}(undef, N_θΨ)

				for i =1:N_θΨ
					FieldName_String[i] = string(param.hypix.plot.θΨ_Table[i] * cst.Mm_2_Cm) * "cm"
				end
				pushfirst!(FieldName_String, string("Horizon")) # Write the "Id" at the very begenning
			
			# Computing θ at required θ
				θ_Mod = Array{Float64}(undef, (N_iHorizon, N_θΨ))
				for iZ=1:N_iHorizon, iΨ =1:N_θΨ
						Ψ_Mod = param.hypix.plot.θΨ_Table[iΨ]
						θ_Mod[iZ, iΨ] = wrc.Ψ_2_θDual(Ψ_Mod, iZ, hydroHorizon)
				end # iZ

			# Concatenating the 2 matrices
			Id = 1:1:N_iHorizon

			θ_Mod = hcat(Id, θ_Mod)

			# Writting the table
				open(Path, "w") do io
					DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
					for iZ = 1:N_iHorizon
						# DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
						DelimitedFiles.writedlm(io, [θ_Mod[iZ, 1:N_θΨ+1]], ",")
					end # i
				end # Path
			
		end  # function:  θΨK_PSD


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KΨ
	# 		Tabular values of the hydroParam model
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KΨ(hydroHorizon, iSim, N_iHorizon)
				
			Path = path.Table_KΨ * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			N_θΨ = Int64(length(param.hypix.plot.θΨ_Table))

			# Writting the Header
				FieldName_String = Array{String}(undef, N_θΨ)

				for i =1:N_θΨ
					FieldName_String[i] = string(param.hypix.plot.θΨ_Table[i] * cst.Mm_2_Cm) * "cm"
				end
				pushfirst!(FieldName_String, string("Horizon Cm/H")) # Write the "Id" at the very begenning
			
			# Computing θ at required θ
				K_Mod = Array{Float64}(undef, (N_iHorizon, N_θΨ))
				for iZ=1:N_iHorizon, iΨ =1:N_θΨ
						Ψ_Mod = param.hypix.plot.θΨ_Table[iΨ]
						K_Mod[iZ, iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Mod, iZ, hydroHorizon) .* cst.MmS_2_CmH
				end # iZ

			# Concatenating the 2 matrices
			Id = 1:1:N_iHorizon

			K_Mod = hcat(Id, K_Mod)

			# Writting the table
				open(Path, "w") do io
					DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
					for iZ = 1:N_iHorizon
						# DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
						DelimitedFiles.writedlm(io, [K_Mod[iZ,1:N_θΨ+1]], ",")
					end # i
				end # Path
			
		end  # function:  θΨK_PSD


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PERFORMACE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iSim, iSim, RmseBest, SwcRoots, WofBest, ΔRunTimeHypix, ΔT_Average)
			Path = path.Table_Performance * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			Header = ["Id" "WofBest" "RmseBest" "Efficiency" "Global_WaterBalance" "Global_WaterBalance_NormPr" "ΔT_Average" "∑∑ΔSink" "∑ΔQ_Bot" "SwcRoots" "iNonConverge" "ΔRunTimeHypix"]

			Id = 1:1:length(WofBest)

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				DelimitedFiles.writedlm(io, [Id WofBest RmseBest Efficiency Global_WaterBalance Global_WaterBalance_NormPr ΔT_Average ∑∑ΔSink ∑ΔQ_Bot SwcRoots iNonConverge_iSim ΔRunTimeHypix], ",")
			end
		end # function PERFORMACE


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SIGNATURE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SIGNATURE(iSim, Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim)
			
			Path = path.Table_Signature * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			Header = ["Month" "Sign_Deficit_Obs" "Sign_Max_Obs" "Sign_Saturated_Obs" "Sign_Senescence_Obs" "Sign_Deficit_Sim" "Sign_Max_Sim" "Sign_Saturated_Sim" "Sign_Senescence_Sim"]

			Month = 1:1:12

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				DelimitedFiles.writedlm(io, [Month Signature_Deficit_Obs Signature_Max_Obs Signature_Saturated_Obs Signature_Senescence_Obs Signature_Deficit_Sim Signature_Max_Sim Signature_Saturated_Sim Signature_Senescence_Sim], ",")
			end
		end  # function: SIGNITURE

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : DAILY_CLIMATE
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DAILY_CLIMATE(∑T_Climate, clim, iSim)
			Path = path.Table_DailyClimate * "_" * string(iSim) * ".csv"
			println("			~  $(Path) ~")

			local ∑T_Int = ceil.(Int, ∑T_Climate[1:clim.N_Climate] .* cst.Second_2_Day)

			Header = ["Year" "Month" "Day" "Pr" "Pr_Ground"]

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				DelimitedFiles.writedlm(io, [year.(clim.Date[1:clim.N_Climate]) month.(clim.Date[1:clim.N_Climate]) day.(clim.Date[1:clim.N_Climate]) clim.Pr[1:clim.N_Climate] clim.Pr_Through] , ",")
			end
			
		end  # function: DAILY_CLIMATE

	
end  # module tableHypix
# ............................................................