# =============================================================
#		MODULE: table
# =============================================================
module table
	using ..path
	import DelimitedFiles
	export θΨK

	# ===================================================
    #          Discretization
    # ===================================================
	function θΨK(Id_Select, Of, Of_θΨ, Of_Kunsat, N_SoilSelect, hydro)
		Header =  ["Id" "Thetas" "Thetar" "Ks" "Sigma" "Hm" "ThetasMat" "SigmaMac" "HmMac" "Nse" "Nse_ThetaH" "Nse_Kunsat"]

		Nse = 1 .- Of
		Nse_θΨ = 1 .- Of_θΨ
		Nse_Kunsat = 1 .- Of_Kunsat

		DelimitedFiles.writedlm(path.Table_θΨK, [Header; Id_Select hydro.θs hydro.θr hydro.Ks hydro.σ hydro.Ψm hydro.θsMat hydro.σMac hydro.ΨmMac Nse Nse_θΨ Nse_Kunsat] , ",")

		println("Nse = $(sum(Nse))")
	 end # Table DISCRETIZATION

	
	
end  # module table
# ............................................................