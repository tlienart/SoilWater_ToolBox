# =============================================================
#		module: kunsatModel
# =============================================================
module kunsatModel

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSAT_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""
		function KUNSAT_MODEL(N_SoilSelect, option, param; IsTopsoil=false )

				KunsatModel_Lab = fill(0.0::Float64, N_SoilSelect)
				if option.hydro.HydroModel==:Kosugi 
				println("	=== START: Dering Ks from lab θ(Ψ) data ===")
					for iZ=1:N_SoilSelect
						if option.dataFrom.Smap
							if smap.IsTopsoil[iZ] == 1
								TopsoilSubsoil="Topsoil"
							else
								TopsoilSubsoil="Subsoil"
							end
						else
							TopsoilSubsoil="Topsoil"
						end
						KunsatModel_Lab[iZ] = kunsat.θΨ_2_KUNSAT(option.hydro, param, 0.9999, iZ, hydro, 0.0; TopsoilSubsoil=TopsoilSubsoil)

						if option.hydro.KunsatΨ == false
							hydro.Ks[iZ] = KunsatModel_Lab[iZ]
						end
					end # for
				println("		~~~ END: Dering Ks from lab θ(Ψ) data ~~~ \n")
				end
			
			return hydro
		end  # function: KUNSAT_MODEL
		
end  # module: kunsatModel
# ............................................................