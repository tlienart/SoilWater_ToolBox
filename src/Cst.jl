# =============================================================
#		MODULE: cst
# =============================================================
   module cst
      const Cm_2_Mm        = 10.0
      const Day_2_Second   = 86400.0
      const Hour_2_Second  = 3600.0
      # const KunsatModel    = 119980.32407407407 # [mm s ⁻¹]
      const KunsatModel    = 9595.085098 # [mm⁻¹ s ⁻¹]

      const MmS_2_CmH      = 3600. / 10.0
      const MmS_2_MmH      = 3600.0
      const Mm_2_Cm        = 0.1
      const Mm_2_kPa       = 0.01
      const Mm_2_kPa_Exact = 1.0 / 101.97162
      const Second_2_Day   = 1.0 / 86400.0
      const Second_2_Hour  = 1.0 / 3600.0
      const Y              = 14.9 # Young Laplace equation constant [mm2]
      const kPa_2_Mm       = 100.0
      const kPa_2_Mm_Exact = 101.97162
      const β              = 0.6 # Best infiltration parameter
   end  # module cst
# ............................................................