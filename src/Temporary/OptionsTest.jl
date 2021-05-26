# =============================================================
#		module: test
# =============================================================
module test

   mutable struct OPTION
      hydro
      infilt
      P5
   end
   
   mutable struct HYDRO
      θs
      θr
   end
   
   mutable struct INFILT
      P1
      P2
   end

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : test
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function TEST()
      P1=1
      P2=3
      P5=10

      θs =3
      θr=4

      infilt = INFILT(P1,P2)
      hydro = HYDRO(θs, θr)
      option = OPTION(hydro, infilt, P5)   
      
      return infilt, hydro, option
   end  # function: name

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : TEST2
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function TEST2()
      infilt, hydro, option = TEST()

      @show hydro.θs
      @show hydro.θr
      @show infilt.P1
      @show infilt.P2

      option.P5 = 10

      option.hydro.θs =5
      option.hydro.θr =6

      option.infilt.P1=5
      option.infilt.P2 =6

      @show option.hydro.θs
      @show option.hydro.θr

      @show  option.infilt.P1
      @show option.infilt.P2

      @show  option.P5

      return
   end  # function: TEST2
end  # module: test

test.TEST2()
# ............................................................