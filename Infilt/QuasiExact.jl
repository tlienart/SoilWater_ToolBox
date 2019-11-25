# =============================================================
#		MODULE: quasiExact
# =============================================================
module quasiExact

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILT3D_2_INFILT1D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILT3D_2_INFILT1D(iSoil, Infilt_3D, Sorptivity, Tinfilt, infiltParam, hydro)
			Δθ = hydro.θs[iSoil] - hydro.θr[iSoil]
			return max(Infilt_3D - ( Tinfilt * infiltParam.γ[iSoil] * Sorptivity ^ 2.0) / (infiltParam.RingRadius[iSoil] * Δθ), eps())
		end  # function: INFILT3D_2_INFILT1D


    #= COMPUTE NORMALISED TIME: Tη =#
    function TIME_2_TIMEη(iT, Sorpt, ΔK)
        return Time_η = iT * 2.0* (ΔK/ Sorpt)^2.
    end




    function INFILTRATION3D_2_INFILTRATION1D(iT, Inf_3D_Obs, Sorpt, RingRadius, θs, θr)
		Δθ = θs - θr
        return Inf_1D = max(Inf_3D_Obs - (iT * cst.γ * Sorpt^2.) / (RingRadius * Δθ), cst.ϵ)
    end



    # TRANSFORMS NORMALIZED INFILTRATION TO INFILTRATION-3D
    function INFILTTRATIONη_2_INFILTRATION3D(iT, Time_η, Inf_η, Sorpt, ΔK,  K_θini, Δθ, RingRadius)
        # Function compute infiltration-1d from normalized infiltration
        function INFILTTRATIONη_2_INFILTRATION1D(iT, Inf_η, Sorpt, ΔK,  K_θini)
            return INFILTRATION_1D = K_θini*iT + Inf_η * (Sorpt^2.) / (2.0 *ΔK)
        end
        ΔI_η = Time_η * cst.γ
       Inf_3D_Obs = INFILTTRATIONη_2_INFILTRATION1D(iT, Inf_η, Sorpt, ΔK, K_θini) + ΔI_η * (Sorpt^4.) / (2.0*RingRadius*Δθ*(ΔK^2.))
        return Inf_3D_Obs
    end



    # TRANSFORMS INFILTRATION-3D TO NORMALIZED INFILTRATION
    function INFILTRATION3D_2_INFILTTRATIONη(iT,Inf_3D_Obs, Sorpt, ΔK, K_θini, Δθ, RingRadius)
          Inf_η = max((2.0 *ΔK / Sorpt^2.) * (Inf_3D_Obs  - K_θini*iT -  cst.γ * TIME_2_TIMEη(iT, Sorpt, ΔK) * (Sorpt^4.) /  (RingRadius * Δθ * 2.0* (ΔK^2.)) ), cst.ϵ)
         return Inf_η
    end

	
end  # module: quasiExact
# ............................................................