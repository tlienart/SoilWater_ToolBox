import Pkg

function PACKAGES(Option_PackageUpdate)
    """This would add the packages automatically if not available"""

    function PACKAGE_MANAGER(Package)
    """ONLY ADD PACKAGE IF NOT AVAILABLE"""
        if  haskey(Pkg.dependencies() , Package) == false
            println("Adding $Package package because not available...")
            Pkg.add(Package)
            println("$Package package added")
        end
    end # PACKAGE_MANAGER

	# PACKAGES
		PACKAGE_MANAGER("BlackBoxOptim")
		PACKAGE_MANAGER("CSV")
		PACKAGE_MANAGER("CairoMakie")
		PACKAGE_MANAGER("DataFrames")
		PACKAGE_MANAGER("ForwardDiff")
		PACKAGE_MANAGER("GLMakie")
		PACKAGE_MANAGER("LaTeXStrings")
		PACKAGE_MANAGER("Makie")
		PACKAGE_MANAGER("Optim")
		PACKAGE_MANAGER("PGFPlots")
		PACKAGE_MANAGER("PGFPlotsX")
		PACKAGE_MANAGER("Plots")
		PACKAGE_MANAGER("Polynomials")
		PACKAGE_MANAGER("QuadGK")
		PACKAGE_MANAGER("Revise")
		PACKAGE_MANAGER("SpecialFunctions")
		PACKAGE_MANAGER("Statistics")
		PACKAGE_MANAGER("Tables")
 

	if Option_PackageUpdate
		println("Updating metdata...")
		Pkg.update()
	end
end # PACKAGES

PACKAGES(false)