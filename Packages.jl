import Pkg

function PACKAGES(Option_PackageUpdate)
    """This would add the packages automatically if not available"""

    function PACKAGE_MANAGER(Package)
    """ONLY ADD PACKAGE IF NOT AVAILABLE"""
        if  haskey(Pkg.installed(), Package) == false
            println("Adding $Package package because not available...")
            Pkg.add(Package)
            println("$Package package added")
        end
    end # PACKAGE_MANAGER

    # PACKAGE_MANAGER("JuliaDB")
    PACKAGE_MANAGER("SpecialFunctions")
    PACKAGE_MANAGER("BlackBoxOptim")
    PACKAGE_MANAGER("Optim")
    PACKAGE_MANAGER("DataFrames")
    PACKAGE_MANAGER("CSV")
	PACKAGE_MANAGER("Statistics")
	PACKAGE_MANAGER("QuadGK")
	PACKAGE_MANAGER("Suppressor")
	PACKAGE_MANAGER("Compat")
    PACKAGE_MANAGER("PGFPlots")
	# PACKAGE_MANAGER("Revise")
    # PACKAGE_MANAGER("ORCA")
    # PACKAGE_MANAGER("PlotlyJS")
    # PACKAGE_MANAGER("GR")
    # PACKAGE_MANAGER("LaTeXStrings")
    # PACKAGE_MANAGER("Polynomials")

    if Option_PackageUpdate
        println("Updating metdata...")
        Pkg.update()
    end

end # PACKAGES

PACKAGES(false)