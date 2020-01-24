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

	# PACKAGES
		PACKAGE_MANAGER("SpecialFunctions")
		PACKAGE_MANAGER("BlackBoxOptim")
		PACKAGE_MANAGER("Optim")
		PACKAGE_MANAGER("CSV")
		PACKAGE_MANAGER("Statistics")
		PACKAGE_MANAGER("QuadGK")
		PACKAGE_MANAGER("Suppressor")
		# PACKAGE_MANAGER("Distributed")
		PACKAGE_MANAGER("Plots")
		PACKAGE_MANAGER("LaTeXStrings")
		# PACKAGE_MANAGER("PGFPlots")

			if Option_PackageUpdate
				println("Updating metdata...")
				Pkg.update()
			end
end # PACKAGES

PACKAGES(false)