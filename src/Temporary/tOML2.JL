Let's add a new (failing) test for getting the General registry path:

@testset "Interacting with the registry" begin 
  @testset "The General registry is accessible" begin 
    IssueReporter.generalregistrypath() |> Base.Filesystem.isdir 
  end 
end 
As for the implementation, we will want to loop over each entry in the DEPOT_PATH and check if it contains a registries/General path of directories. These should be in the user depot, but a more extensive lookup will make our code more robust:

function generalregistrypath() 
  for i in DEPOT_PATH 
    if isdir(joinpath(i, "registries", "General")) 
      return joinpath(i, "registries", "General") 
    end 
  end 
end 
Once we have the path to the General registry, we'll want to parse the Registry.toml file and extract the information corresponding to the package we'll be searching for. Once parsed, the Registry.toml file produces in a dictionary with five entries:

Dict{String,Any} with 5 entries: 
  "name"        => "General" 
  "repo"        => "https://github.com/JuliaRegistries/General.git" 
  "packages"    => Dict{String,Any}("c786d6c3-4fbc-59fc-968c-e848efb65d2d"=>Dict{String,Any}("name"=>"ScHoLP","path"=>"S/ScHoLP"),"88634af6-177f-5301-88b8-7819386cfa38"=>Dict{String,Any}("name"=>"SaferIntegers","path"=>"S/SaferIntegers")... 
  "uuid"        => "23338594-aafe-5451-b93e-139f81909106" 
  "description" => "Official general Julia package registry where people  
# output omitted #  
We're only interested in the packages data, which looks like this:

Dict{String,Any} with 2358 entries: 
  "c786d6c3-4fbc-59fc-968c-e848efb65d2d" => Dict{String,Any}("name"=>"ScHoLP","path"=>"S/ScHoLP") 
  "88634af6-177f-5301-88b8-7819386cfa38" => Dict{String,Any}("name"=>"SaferIntegers","path"=>"S/SaferIntegers") 
  "aa65fe97-06da-5843-b5b1-d5d13cad87d2" => Dict{String,Any}("name"=>"SnoopCompile","path"=>"S/SnoopCompile") 
# output truncated # 
And actually, we don't want all that, since we don't care about the UUID; only the name and the path. Let's add a new function to the IssueReporter module, to reflect this spec:

function generalregistry() 
    TOML.parsefile(joinpath(generalregistrypath(), "Registry.toml"))["packages"] |> values |> collect 
end 
The output of the function is similar to this, an array of Dict elements:

2358-element Array{Any,1}: 
 Dict{String,Any}("name"=>"ScHoLP","path"=>"S/ScHoLP") 
 Dict{String,Any}("name"=>"SaferIntegers","path"=>"S/SaferIntegers") 
 Dict{String,Any}("name"=>"SnoopCompile","path"=>"S/SnoopCompile") 
# output truncated # 
Once we have this, it is very easy to perform a package lookup by name. We simply iterate over each item and compare the "name" value against the search string:

function searchregistry(pkgname::String) 
  for item in generalregistry() 
    item["name"] == pkgname && return item 
  end 
end 
After we identify a package by name, we can use the path value to build the path to the folder that contains the package's metadata information. Remember that we're after the Package.toml file, as this contains the repo URI.

Putting it all together, we can finally write our IssueReporter.packageuri function:

function packageuri(pkgname::String) 
  TOML.parsefile(joinpath(generalregistrypath(), searchregistry(pkgname)["path"], "Package.toml"))["repo"] 
end 
Your IssueReporter.jl file should look like this:

module IssueReporter 
 
using Pkg, Pkg.TOML 
 
function generalregistrypath() 
  for i in DEPOT_PATH 
    if isdir(joinpath(i, "registries", "General")) 
      return joinpath(i, "registries", "General") 
    end 
  end 
end 
 
function generalregistry() 
    TOML.parsefile(joinpath(generalregistrypath(), "Registry.toml"))["packages"] |> values |> collect 
end 
 
function searchregistry(pkgname::String) 
  for item in generalregistry() 
    item["name"] == pkgname && return item 
  end 
end 
 
function packageuri(pkgname::String) 
  TOML.parsefile(joinpath(generalregistrypath(), searchregistry(pkgname)["path"], "Package.toml"))["repo"] 
end 
 
end # module 