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
	function θΨK(Of, Of_θΨ, Of_Kunsat, N_SoilSelect, hydro)
		Header =  ["Thetas" "Thetar" "Ks" "Sigma" "Hm" "ThetasMat" "SigmaMac" "HmMac" "Of" "Of_ThetaH" "Of_Kunsat"]

		DelimitedFiles.writedlm(path.Table_θΨK, [Header; hydro.θs hydro.θr hydro.Ks hydro.σ hydro.Ψm hydro.θsMat hydro.σMac hydro.ΨmMac Of Of_θΨ Of_Kunsat] , ",")
	 end # Table DISCRETIZATION

	
	
end  # module table
# ............................................................