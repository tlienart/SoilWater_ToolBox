function PATH_LOCATION()
   Home = @__DIR__
   println(Home)

   Home2 = dirname(Home)
   println(Home2)
end

cd(raw"D:\Main\MODELS\SoilWater_ToolBox\data\INPUT\Data_Hypix\LYSIMETERS\OPTIMISATION")
A=readdir()
sort(A)

PATH_LOCATION()