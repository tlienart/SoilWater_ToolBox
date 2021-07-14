<!-- MathJax -->

<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>


# POTENTIAL EVAPOTRANSPIRATION MODEL

The *potential evapotranspiration* depth, $\varDelta Pet$ [L], for a time-step $\varDelta T$ is computed using the Penman–Monteith equation. Input variables used the estimates derived from New Zealand Virtual Climate Stations network (VCS), which are based on the spatial interpolation of actual data observations made at climate stations located around the country ([Tait *et al*., 2006](#_ENREF_4)). We assume that:

 $$ \varDelta Pet=\varDelta Pet_{int} $$
 
where $\varDelta Pet_{int}$ [L] is the *interception potential evaporation* depth, which is the potential evaporation from a wet canopy. The remaining energy that is not used to evaporate water from a wet canopy (e.g., periods with no rainfall) allows the *potential evapotranspiration* depth, $\varDelta Pet_{et}$ [L], to be computed as:

$$\varDelta Pet_{et}=\varDelta Pet-\varDelta Evap_{int}$$

where $\varDelta Evap_{int}$ [L] is derived from the actual energy used to evaporate water from a wet canopy.

The *potential transpiration* depth of vegetation, $\varDelta Pet_{transp}$ [L], and *potential evaporation* depth of soil, $\varDelta Pet_{evap}$ [L], is partitioned from $\varDelta Pet_{et}$ by using the Beer–Lambert law, which uses as parameters the leaf area index, $LAI$ [-]. The Beer–Lambert law assumes that the net radiation inside the canopy decreases exponentially. Therefore, the partitioning of $\varDelta Pet_{et}$ is given by:

$$\varDelta Pet_{evap}=\varDelta Pet_{et}\,\,e^{-K_{\mathrm{g}}×LAI} $$
$$\varDelta Pet_{transp}=\varDelta Pet_{et}^{}-\varDelta Pet_{evap} $$

where the extinction coefficient for solar radiation, $K_{g}$ [-], is set to 0.5 (e.g. [Varado *et al*., 2006](#_ENREF_5)).
