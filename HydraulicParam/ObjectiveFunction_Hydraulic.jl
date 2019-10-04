module of 
      include("Option.jl")
      include("Param.jl")
      include("Cst.jl")
      include("Stats.jl")
      include("Wrc.jl")
      include("Kunsat.jl")

	  export  WRC_KUNSAT
	  
      # OF BIMODAL and UNIMODAL
      function WRC_KUNSAT(Option_Data_Kθ, iSoil, Ψ_θΨ, θ_θΨ, N_θΨ, K_Kθ, Ψ_Kθ, N_Kθ, θsMac, θr, σMat, ΨkgMat, θsMat, σMac, ΨkgMac, KsMac) 
         # OF θΨ==
         θ_Obs = zeros(Float64, N_θΨ[iSoil])
         θ_Sim = zeros(Float64, N_θΨ[iSoil])
		 @simd for iH in 1:N_θΨ[iSoil]
            θ_Obs[iH]= θ_θΨ[iSoil,iH]
            H_Obs = Ψ_θΨ[iSoil,iH]
            θ_Sim[iH] = wrc.kg.Ψ_2_θdual(H_Obs, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)
         end

		 Of_θh = stats.NASH_SUTCLIFFE_ERRORmin(θ_Obs[1:N_θΨ[iSoil]], θ_Sim[1:N_θΨ[iSoil]], Power=2.0)

         # OF Kunsat==
         if Option_Data_Kθ
            Kunsat_Obs = zeros(N_Kθ[iSoil])
            Kunsat_Sim = zeros(N_Kθ[iSoil])
            @simd for iH in 1: N_Kθ[iSoil]
               Kunsat_Obs[iH] = K_Kθ[iSoil,iH]
               H_Obs =  Ψ_Kθ[iSoil,iH]
               θ_Sim = wrc.kg.Ψ_2_θdual(H_Obs, θsMac, θr, ΨkgMat, σMat, θsMat, ΨkgMac, σMac)
               Se = wrc.se.θ_2_Se(θ_Sim, θsMac, θr)
      
			   Kunsat_Sim[iH] = kunsat.kg.Se_2_KUNSAT(Se, θsMac, θr, σMat, KsMac, θsMat, σMac)
            end

			Of_Kunsat = stats.NASH_SUTCLIFFE_ERRORmin(log.(1.0 .+ Kunsat_Obs[1:N_Kθ[iSoil]]) , log.(1.0 .+ Kunsat_Sim[1:N_Kθ[iSoil]]))
         else
            Of_Kunsat = 0.0

         end

         return Of = param.Of_W * Of_θh + (1.0 - param.Of_W) * Of_Kunsat

      end # WRC_KUNSAT (OF BIMODAL and UNIMODAL) 

end # module objective function