# =============================================================
#		module: θINITIAL()
# =============================================================
module θini
   import ..wrc
   export θINI_TOP

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θINI_TOP
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   """ 
      θINI_TOP(hydro, NiZ, θ)
   Compute θ[t=1:1:NiZ] when only θ[1,1] is known
   It assumes that Se[t=1:1:NiZ] is constant through the soil profile"""
      function θINI_TOP(hydro, NiZ, θ)
         Se = wrc.θ_2_Se(θ[1,1], 1, hydro)

         for iZ=2:NiZ
            θ[1,iZ] = wrc.Se_2_θ(Se, 1, hydro) 
         end
      return θ
      end  # function: θINI_TOP

end  # module: θini
# ............................................................