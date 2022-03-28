# =============================================================
#		module: θaver
# =============================================================
module θaver
   import ..discretisation
   export θAVER

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θ_AVERAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, NiZ=NiZ, Nit_Reduced=Nit_Reduced, Zaver=250.0)
         # Make sure that Zaver is whitin physical bounds
            Zaver = min(Zaver, Z[NiZ])

         # Memory
            θsim_Aver = fill(0.0::Float64, Nit_Reduced)

         # For every time step
         for iT = 1:Nit_Reduced
            iZ_Max = 1
            for iZ = 1:NiZ
               if Z[iZ] ≤ Zaver
                  θsim_Aver[iT] += θ_Reduced[iT, iZ] * discret.ΔZ[iZ]
               else
                  iZ_Max = iZ
                  break
               end
               iZ_Max = iZ
            end # iZ=1:NiZ
            if Z[iZ_Max-1] + eps(100.0) < Zaver
               θsim_Aver[iT] += θ_Reduced[iT, iZ_Max] * (Zaver - Z[iZ_Max-1])
            end
            
            θsim_Aver[iT] =  θsim_Aver[iT] / Zaver
         end # for iT

      return θsim_Aver
      end  # function: θ_AVERAGE
end  # module: θaver
# ............................................................

