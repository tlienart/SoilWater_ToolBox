# =============================================================
#		module: ofKsModel
# =============================================================
module optKsModel
	import ..stats

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KSMODEL(hydro, hydroParam, KsModel, N_iZ)
			return Of_Ks = stats.NASH_SUTCLIFE_MINIMIZE(log1p(hydro.Ks[1:N_iZ]) , log1p(KsModel[1:N_iZ]))	
		end  # function: OF_KSMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
end  # module: ofKsModel
# ............................................................