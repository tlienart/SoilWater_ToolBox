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


# SINK TERM MODEL

The *sink term*, $\varDelta Sink$, described in the Richards equation is a function of the *soil evaporation* depth, $ \varDelta Evap $ [L], and the *root water uptake* depth of the vegetation, $\varDelta Rwu$ [L]:

$$\begin{equation}
\varDelta Sink=\varDelta Evap+\varDelta Rwu
\end{equation}$$

##	*Soil evaporation*: $\varDelta Evap$

The evaporation model based on [Romano and Giudici (2009)](#_ENREF_6) and adapted by [Pollacco and Mohanty (2012)](#_ENREF_7) is computed as:

$$\begin{equation}
\varDelta Evap\,\,=\,\,\varDelta Pet_{evap}^{}\,\,Se_{1}^{}
\end{equation}$$

where $Se_{1}$ [0-1] [-] refers to the *effective soil water content* defined for the top cell. $\varDelta Pet_{evap}$ [L] is the *potential evaporation* depth computed with the Beer–Lambert law.

##	*Root water uptake*: $\varDelta Rwu$

The *root water uptake* computes the volume of water removed per unit time from a unit volume of soil in the root zone and is computed for each cell $i$:

$$\begin{equation}
\varDelta Rwu_i=K_c\cdot \varDelta Pet_{transp}^{}\cdot \varDelta Rdf_i\cdot FwaterStress_i\cdot RootComp_i 
\end{equation}$$

where $\varDelta Rwu$ [L] is the *root water uptake* depth; $Kc$ [-] is the *crop coefficient*; $\varDelta Pet_{\mathrm{transp}}$ [L] corresponds to the *potential transpiration* depth; $\varDelta Rdf_{i}$ [-] refers to the percentage of roots per cell; $FwaterStress$ [-] is the *water stress function* per cell, which computes the reduction of transpiration based on $ψ$; and $RootComp$ [-] is the *root compensation* by enabling water uptake from deeper layers when the upper layers are depleted.

### *Root density function*: $\varDelta RootDensity$

The percentage of roots per cell is given by $Rdf_{i}$, which defines the general shape of the roots, as given by the example provided in [Figure 6](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure6.bmp). The root distribution is based on an empirical function of [Gale and Grigal (1987)](#_ENDREF_8), which was modified by [Pollacco et al. (2008)](#_ENDREF_9). The model requires, as parameters, the *maximum rooting depth*, $Z_{Nroot}$ [L] and the percentage of roots in the top 30 cm (other values can also be considered), $\varDelta Rdf_{Top}$ [%].

The fraction of roots, $Rdf_{i}$ [-], for each cell, $i$, is computed as:

$$\begin{equation}
\begin{cases}	\varDelta Rdf_i=\frac{R_d^{Z_{i+\frac{1}{2}}}-R_d^{Z_{i-\frac{1}{2}}}}{1-R_d^{Z_{Nroot}}}\\	with\,\,\sum_{i=1}^{i=Nroot}{\varDelta Rdf_i=1}\\\end{cases}
\end{equation}$$

where $Z_{i-\frac{1}{2}}$ and $Z_{i+\frac{1}{2}}$ [L] are respectively the depth of the top and the bottom of cell $i$; $Rd$ [-] is the routing distribution parameter (computed numerically from $Z_{Nroot}$ and $\varDelta Rdf_{Top}$). $R_{d}$ varies between 0.7000 and 0.9999, such that when $R_{d}$ is close to 0.7, all the roots are distributed at the top, and when $R_{d}$ is close to 1, the roots are evenly distributed within the root zone.

The value of $R_{d}$ is computed by solving the following equation:

$$\begin{equation}
\varDelta Rdf_{Top}=\frac{R_d^0-R_d^{30}}{1-R_d^{\left| Z_{root} \right|}}=\frac{1-R_d^{30}}{1-R_d^{\left| Z_{root} \right|}}
\end{equation}$$

An example of $\varDelta Rdf$ is provided in [Figure 6a](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure6.bmp) for a soil that has equal discretization.

![HyPix](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure6.bmp "Figure 6. Pasture grass models used for all lysimeters: (a) root density function plotted with depth and (b) schematic of the [Feddes et al. (1978)](#_ENDREF_10) plant water stress function.")


#### *Water Stress response function*: $F_{WaterStress}$

The *water stress response function*, $F_{WaterStress}$ [-] shown in [Figure 6b](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure6.bmp) is a prescribed dimensionless function of $ψ$, where (0 $\le$ $F_{WaterStress}$ $\le$ 1). Following [Feddes et al. (1978)](#_ENDREF_10), $F_{WaterStress}$ is defined by using four soil water pressure values leading to a trapezoidal curve, so that:

* $F_{WaterStress} = 0$ close to saturation at a pressure greater than $\varPsi _{\mathrm{Feddes}1}$ and also above the permanent wilting point soil water pressure, $\varPsi _{\mathrm{Feddes}4}$
*  $F_{WaterStress} = 1$ between soil water pressures $\varPsi _{\mathrm{Feddes}2}$ and $\varPsi _{\mathrm{Feddes}3}$
* for soil water pressure between $\varPsi _{\mathrm{Feddes}1}$ and $\varPsi _{\mathrm{Feddes}2}$, $F_{WaterStress}$ increases linearly
* for soil water pressure between $\varPsi _{\mathrm{Feddes}3}$ and $\varPsi _{\mathrm{Feddes}4}$, $F_{WaterStress}$ decreases linearly.

A schematic plot of this stress response function is depicted in [Figure 6b](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure6.bmp).


#### *Compensation mechanism*: $Root_{comp}$

A *root water uptake compensation* mechanism is introduced to improve the prediction of transpiration by enabling water uptake from deeper layers when the upper layers are depleted, although the percentage of roots at deeper depth is limited. The compensation mechanism of [Li et al. (2001)](#_ENDREF_11), validated by [Braud et al. (2005)](#_ENDREF_12), is introduced. The model requires the *compensation mechanism parameter* $C$ [-], which accounts for the general soil water content profile before computing the water uptake from individual cell i and is derived as:

$$\begin{equation}
\varDelta RootComp_{i}^{}=\frac{FwaterStress_i\,\,\varDelta Rdf_{i}^{C-1}}{\sum_{i=1}^{i=N_{root}}{FwaterStress_i\,\,\varDelta Rdf_i^C}}
\end{equation}$$

where $Rdf_{i}$ is the vertical fraction of the root density function for each $i$ cell [%]; $F_{WaterStress}$ is the reduction of root water uptake at pressure head $ψ$ for every cell; and $C$ is a parameter such that when $C$ = 1 the model is not compensated, and when $C$ = 0 the $Rdf_{i}$ becomes constant throughout the whole root-zone depth ($i=N_{root}$). In HyPix, $C$ = 0.5, as suggested by [Li et al. (2001)](#_ENDREF_11) and [Braud et al. (2005)](#_ENDREF_12).
