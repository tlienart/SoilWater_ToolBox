# =============================================================
#		module: θaver
# =============================================================
module θaver
   import ..discretization
   export θAVER

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θ_AVERAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θAVER(discret; Z=Z, θ_Plot=θ_Plot, N_iZ=N_iZ, N_∑Treduced=N_∑Treduced, Zaver=250.0)
         # Make sure that Zaver is whitin physical bounds
            Zaver = min(Zaver, Z[N_iZ])

         # Memory
            θsim_Aver = fill(0.0::Float64, N_∑Treduced)

         # For every time step
         for iT = 1:N_∑Treduced
            iZ_Max = 1
            for iZ = 1:N_iZ
               if Z[iZ] ≤ Zaver
                  θsim_Aver[iT] += θ_Plot[iT, iZ] * discret.ΔZ[iZ]
               else
                  iZ_Max = iZ
                  break
               end
               iZ_Max = iZ
            end # iZ=1:N_iZ
            if Z[iZ_Max-1] + eps(100.0) < Zaver
               θsim_Aver[iT] += θ_Plot[iT, iZ_Max] * (Zaver - Z[iZ_Max-1])
            end
            
            θsim_Aver[iT] =  θsim_Aver[iT] / Zaver

         # println(iT," , ",  " , ",θsim_Aver[iT]," ; ",θ_Plot[iT,5] )
         end # for iT

      return θsim_Aver
      end  # function: θ_AVERAGE
end  # module: θaver
# ............................................................

