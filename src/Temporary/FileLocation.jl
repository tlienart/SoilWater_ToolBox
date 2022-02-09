function PATH_LOCATION()
   Home = @__DIR__
   println(Home)

   Home2 = dirname(Home)
   println(Home2)
end

PATH_LOCATION()