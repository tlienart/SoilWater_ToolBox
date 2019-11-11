# =============================================================
#		MODULE: table
# =============================================================
module table


	# =============================================================
	#		MODULE: name
	# =============================================================
	module hydroParam #TODO need to finish
		import ...path
		import DelimitedFiles, CSV, DataFrames
		export θΨK


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨK(Id_Select, Of, Of_θΨ, Of_Kunsat, N_SoilSelect, hydro)
			# Header =  ["Id" "Thetas" "Thetar" "Ks" "Sigma" "Hm" "ThetasMat" "SigmaMac" "HmMac" "Nse" "Nse_ThetaH" "Nse_Kunsat"]

			Matrix = Array{Float64}(undef, ( N_SoilSelect, length(hydro.FieldName)-1))
			for i=1:length(hydro.FieldName)-1
				Struct_Array = getfield(hydro, hydro.FieldName[i])
				Matrix[1:N_SoilSelect,i] = Struct_Array[1:N_SoilSelect]
			end

			FieldName_String = Array{String}(undef, length(hydro.FieldName))
			i=1

			for FieldNames in hydro.FieldName
				FieldName_String[i] =  String(FieldNames)
				i += 1
			end
			


		println(FieldName_String)


			Nse = 1 .- Of
			Nse_θΨ = 1 .- Of_θΨ
			Nse_Kunsat = 1 .- Of_Kunsat

			# DelimitedFiles.writedlm(path.Table_θΨK, [Header; Id_Select Matrix Nse Nse_θΨ Nse_Kunsat] , ",")

			open(path.Table_θΨK, "w") do io
				DelimitedFiles.writedlm(io,[ FieldName_String] , ",",)
				DelimitedFiles.writedlm(io, [Id_Select Matrix], ",")
			end;
	 
			return
		end  # function:  θΨK
		
	end  # module hydroParam
	# ............................................................


	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		
	end  # module psd
	# ............................................................






	

		




		





	
end  # module table
# ............................................................