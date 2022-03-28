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

# HYPIX OUTPUTS

The following files are table output in HyPix with the following Path for which in this case the project name is TESTCASE for name of site TC1: test

SoilWater_Toolbox/data/OUTPUT/Hypix/TESTCASE/Tc1/TableSoilWater_Toolbox/data/OUTPUT/Hypix/TESTCASE/Tc1/Table

## PARAMETERS

### TABLE_Hydro.csv

These are the optimal Kosugi hydraulic parameters for every layers:

| Parameters    | Units                  | Explantions                                                                                                               |
| ------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Id            | [-]                    | Layers                                                                                                                    |
| θs           | [mm$^3$ mm$^{-3}$] | Saturated volumetric soil water content.                                                                                  |
| θr           | [mm$^3$ mm$^{-3}$] | Residual volumetric soil water content.                                                                                   |
| Ks            | [mm s$^{-1}$]        | Saturated hydraulic conductivity.                                                                                         |
| Ψm           | [mm]                   | Mean of ln ψ in the matrix domaine.                                                                                      |
| σ            | [-]                    | Standard deviation of ln ψ in the matrix domaine.                                                                        |
| θsMacMat\_ƞ | [mm$^3$ mm$^{-3}$] | Normalized volumetric saturated water content that differentiates inter-aggregate pores and matrix domains.               |
| σMac         | [-]                    | Standard deviation of ln ψ in the macropore domaine.                                                                     |
| ΨmMac        | [mm]                   | Mean of ln ψ in the macropore domain.                                                                                    |
| So            | \[mm$^{-1}$\]        | Fluid compressibility.                                                                                                    |
| θsMacMat     | [mm$^3$ mm$^{-3}$] | Volumetric saturated water content that differentiates inter-aggregate pores and matrix domains with θsMacMat$\le$θs. |

### TABLE_Veg.csv

## HYDRAULIC RELATIONSHIPS

### TABLE\_$\theta \Psi$.csv

For every layer, the water retention curve,  $\theta(\Psi)$ , which is the relationship between the *soil water pressure*, $\Psi$ [mm] and the *volumetric soil water content* *θ* [mm$^3$ mm$^{-3}$] , for  values of $\Psi$ selected by the users.

### TABLE\_$\theta K \Psi$.csv

For every layer, the hydraulic conductivity  $K(\Psi)$ , which is the relationship between the *soil water pressure*, $\Psi$ [cm] and the *unsaturated hydraulic conductivity* *K* [mm h$^{-1}$] , for  values of $\Psi$ selected by the users.

## TIME SERIES

### TABLE_TimeSerie.csv

The meaning of Δ*X* means:

$$  \Delta X = \sum_{t=1}^T{X^t-}\sum_{t=1}^{T-1}{X^t} $$

| Id                 | Units    | Explanations                                                                         |
| ------------------ | -------- | ------------------------------------------------------------------------------------ |
| Year               | [Year]   | Year                                                                                 |
| Month              | Month    | Month                                                                                |
| Day                | [Day]    | Day                                                                                  |
| Hour               | [Hour]   | Hour                                                                                 |
| Minute             | [Minute] | Minute                                                                               |
| Second             | [Second] | Second                                                                               |
| ΔPrGross          | [mm]     | Precipitation reaching the top of vegetation                                         |
| ΔPrSoil           | [mm\]    | Precipitation infiltrating into the soil through throughfall and which is not runoff |
| ΔPet              | \[mm\]   | Potential evapotranspiration                                                         |
| ΔSink             | \[mm\]   | Sink term =Transpiration + Evaporation                                               |
| ΔTranspiration    | [mm\]    | Transpiration                                                                        |
| ΔEvaporation      | mm\]     | Evaporation                                                                          |
| ΔRecharge         | \[mm\]   | Recharge at the bottom cell                                                          |
| Hpond              | [mm\]    | Maximum ponding depth between 2 time steps                                           |
| ΔRunoff           | \[mm\]   | Runoff                                                                               |
| ∑WaterBalance\_η | [mm\]    | Normalized water balance                                                             |

## TIME SERIES WITH DEPTH

The following are output with time series at time step interval defined by user with depth. The header are the the depth of the center of the cell (Znode):

### TABLE\_$\theta$.csv

A matrice of time series with depth of volumetric soil water content $\theta$ [m$^3 $ m$^{-3}$].

### TABLE\_$\Psi$.csv

A matrice of time series with depth of volumetric soil water content $\Psi$ [mm].

### TABLE\_$\Delta$Flux.csv

A matrice of time series with depth of flux $\Delta Q$ [mm], where positive is water movong downwards and negative where water is moving upwards

## # HOW GOOD THE SIMULATIONS

### TABLE_Performance.csv
