# =============================================================
#		MODULE: cst
# =============================================================
module cst
	const bd 		= 0.916 		# bulk density [g cm-3] 0.95  Sam's dataset 0.916 # in
	const pd 		= 2.65 			# particle density [g cm-3] 2.59 # Sam's dataset 2.65  # in (g cm-3)
	const γ 		= 0.75 			# Shape param
	const β 		= 0.6  			# 0.60 Ratio defined using the ratio between 2 estimators of the Sorptivity
	const ϵ 		= 0.00000001
	const Y 		= 0.149 * 10.0 ^2. # [mm^2]
	const Kconst 	= (10.0 / (24.0 * 60.0 * 60.)) *( 1.03663 * 10.0 ^9.) #  convert from cm/day to mm s−1.
	const Pσ_1 		= 0.5920
	const Pσ_2 		= 0.7679
end  # module cst
# ............................................................


