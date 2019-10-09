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
		Header =  ["Id" "Thetas" "Thetar" "Ks" "Sigma" "Hm" "ThetasMat" "SigmaMac" "HmMac" "Of" "Of_ThetaH" "Of_Kunsat"]

		DelimitedFiles.writedlm(path.Table_θΨK, [Header; Id_Select hydro.θs hydro.θr hydro.Ks hydro.σ hydro.Ψm hydro.θsMat hydro.σMac hydro.ΨmMac Of Of_θΨ Of_Kunsat] , ",")
	 end # Table DISCRETIZATION

	
	
end  # module table
# ............................................................