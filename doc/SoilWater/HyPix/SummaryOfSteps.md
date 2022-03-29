<!-- MathJax -->
  <script type="text/x-mathjax-config">
    MathJax.Hub.Config({
		     TeX: {
      equationNumbers: {
        autoNumber: "AMS"
      }
    },
      tex2jax: {
        skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
        inlineMath: [['$','$']]
      }
    });
  </script>
<script id="MathJax-script" async src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML"></script>


# Summary of steps in HyPix 

## Initialisation of the vegetation and soil hydraulic parameters

-	Compute the initial soil water pressure, $ψ_{ini}$ [L], from initial observed $θ^{t=1}$.
-	Calculate $\varDelta ψ_{max}$ for every cell, which will be used for the variable time-step [[Eq.(16 FROM RICHARDSEQUATION!!!!)]]().
-	Estimate the percentage of roots, $\varDelta Rdf_i$, per cell, derived from the vegetation parameters [[Eq. (61 FROM VEGETATION PARAMETERS!!!!!)]]().

## LOOPING 1: Compute at a daily time-step:

As described in [Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp), HyPix first computes, for all data:

-	vegetation parameters $LAI$ and $K_C$, which vary monthly (section 3.1.4)
-	precipitation reaching the ground, $\varDelta Pr_{ground}$ [(section 7.4 !!!!!!!!!!]())
-	ΔPettransp and ΔPetevap [(section7.3 !!!!!!!!!!)]().

## LOOPING 2: Iteration of the NR, which solves $\psi$ of the RE for a given $\varDelta T$

### Compute at the beginning of every time-step: 

-	new time-step, $\varDelta T$, computed from $ψ^{t-1}$ [(section 2.2.2.4 !!!!!)](),
-	sorptivity from $ψ^{t-1}$ to derive the infiltration rate [(section 2.2.1.1 !!!!!)](),
-	interpolate $\varDelta Pr_{ground}$ and $\varDelta Pet_{evap}$ for a given $\varDelta T$,
-	evaporation from $ψ^{t-1}$ [(section 7.5.1 !!!!!)](),
-	root water uptake from $ψ^{t-1}$ [(section 7.5.2 !!!!!)]().

### Iteration:

-	The iteration loop of the NR (section 2.2.2) is continued until either **(a)** the convergence criterion, [Eq. (11) FROM RICHARDSEQUATION !!!!!](), is met, or **(b)** the maximum allowed iteration is reached.
### After iteration is computed, determine if recompute is needed using a smaller time-step

-	Rerun the model if $\varDelta T^t(\psi ^t)\ll \varDelta T^t(\psi ^{t-1})$ with a new derived time-step [(section 2.2.2.5 !!!!!!!)]().
### Results

-	Compute the soil water balance [(section 2.2.2.6.1 !!!!!!)]().
-	Compute the efficiency of simulation [(section 2.2.2.6.2 !!!!!!)]().
