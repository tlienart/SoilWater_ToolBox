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


# HYPIX

Modelling unsaturated flow in highly heterogeneous soils can be accurately performed by solving the [Richards (1931)](#_ENDREF_13) equation (RE), which is commonly adopted by soil vegetation atmosphere transfer models. However, RE is highly nonlinear, and despite numerous efforts over the last decade its solution using numerical methods is demanding, and problems finding techniques to achieve fast and accurate solutions are unresolved (e.g. [Zha *et al*., 2019](#_ENDREF_14)). Here, we propose improvements to RE and implement them in HyPi, using the mixed form of RE, as recommended by [Celia *et al*. (1990)](#_ENDREF_14). The solution of RE is based on [Hassane Maina and Ackerer (2017)](#_ENDREF_15), for which the RE partial differential equation is solved using a *cell-centered finite-volume (implicit finite differences)* scheme for the spatial discretization, with an implicit Euler scheme for the temporal discretization by using the weighted average inter-cell hydraulic conductivity.

##	*Richards equation of HyPix model*

Assuming a rigid solid matrix ([Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp)), the mixed form of RE is written as:

$$\begin{equation}
\frac{\theta _i\left( \psi _{i}^{t} \right) -\theta _i\left( \psi _{i}^{t-1} \right)}{\varDelta T^t}-S_o\frac{\theta _i\left( \psi _{i}^{t} \right)}{\theta _{s_i}}\frac{\left| \psi _{i}^{t} \right|-\left| \psi _{i}^{t-1} \right|}{\varDelta T^t}=\frac{Q_{i-\frac{1}{2}}^{t}-Q_{i+\frac{1}{2}}^{t}}{\varDelta Z_i}-Sink_i\left( \psi _{i}^{t-1} \right) 
\end{equation}$$

where $\varDelta T^t$ [T] is the time-step at time $t$; $\varDelta Z_{i}$ [L] is the mesh size of the cell $i$, with the vertical coordinate positive downwards; $θ_{i}$ [L<sup>3</sup> L<sup>-3</sup>] is the volumetric soil water content of the cell $i$; $θ_{s}$ [L<sup>3</sup> L<sup>-3</sup>]  is the saturated volumetric soil water content; $S_{0}$ [L<sup>-1</sup>] is a parameter that accounts for fluid compressibility, which is assumed to be constant with depth; $ψ_{i}$ [L] is the soil water pressure of cell $i$, considering $ψ > 0$ for unsaturated soils; $Q$ [L T<sup>-1</sup>] is the soil water flux based on the extended Darcy’s law, positive downward and negative when water moves upwards; and $Sink_{i}$ [L<sup>3</sup> L<sup>-3</sup>], taken as positive, is the sink term defined as the volume of water per unit time removed from cell $i$ by *evaporation* and *root water uptake*.

![HyPix](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp "Figure 1. Diagram describing the 1D vertical discretization of the Richards equation, where i  is the cell number (cell 1 is the top cell and cell Ni is the bottom cell, therefore cell = i  is below cell = i¬–1); ΔPr  [L] is the precipitation reaching the top of the canopy; ΔPrground [L]  is the precipitation reaching the soil surface (cell = 1); ΔHpond  [L] is the ponding water; and ΔQi+1/2 = Qi+1/2 ΔT [L] is the inter-cell water volume (positive downwards). Water is removed from the soil profile by transpiration, ΔTransp [L], and evaporation, ΔEvapo [L], depending on θ and potential evapotranspiration, ΔPet [L] (partitioned between potential evaporation, ΔPetevap [L], and potential transpiration, ΔPettransp [L]).")

### *Computing water infiltration into the soil profile*

The amount of water infiltrating into the top cell, $i = 1$, for a period of time, $\varDelta T$, is computed by $\varDelta Q_{1/2}$ [L] ([Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/SoilWater-ToolBox-FlowChart.bmp)). As RE cannot compute the top air–soil boundary, we compute $\varDelta Q_{1/2}$ using the two-term approximation of [Haverkamp *et al*. (1994)](#_ENDREF_16), as suggested by [Fernández-Gálvez *et al*. (2019)](#_ENDREF_17). We do not include the 3D radial flux considered in the two-term expansions because HyPix computes 1D water infiltration. The maximum infiltration depth for a given $\varDelta T$ is $\varDelta Qmax_{\frac{1}{2}}^{t}$ [L], and is computed as:

$$\begin{equation}
\begin{cases}	B=K_s\left( \frac{2-\beta}{3}+\frac{1+\beta}{3}\frac{K\left( \theta _{1}^{t-1} \right)}{K_s} \right)\\	\varDelta Qmax_{\frac{1}{2}}^{t}=Qmax_{\frac{1}{2}}^{t}\,\,\varDelta T^t=\cos  \alpha \left[ Sorpt\left( \theta _{1}^{t-1} \right) \,\,\sqrt{\varDelta T^t}\,\, +B\,\,Ks_1\,\,\varDelta T^t \right]\\\end{cases}
\end{equation}$$

where $Qmax_{\frac{1}{2}}^{t}$[L T<sup>-1</sup>] is the maximum soil water flux; $Sorpt$  [L T<sup>-1/2</sup>] is the soil sorptivity; $\beta$ [-] is an integral shape parameter, typically fixed at 0.6 ([Haverkamp *et al*., 1994](#_ENDREF_18); [Parlange *et al*., 1982](#_ENDREF_19)); and the slope, $α$ [radian], is the angle between the flow direction of recharge and the vertical axis ($0 ≤ α ≤ π/2$ for inclined flow). 

We use the physically based sorptivity model of [Lassabatere *et al*. (2021)](#_ENDREF_20), which uses a mixed formulation that was validated against analytical expressions for several types of soils using different models of hydraulic functions. The procedure to compute sorptivity is efficient for all types of hydraulic functions and shape parameters with insignificant errors ([Lassabatere *et al*., 2021](#_ENDREF_20)). The sorptivity model is computed as:

$$\begin{equation}
 Sorpt^2\left( \theta _0,\theta _s \right) =\int_0^{\frac{\theta _s-\theta _r}{2}}{\left( \theta _s+\theta -2\theta _0 \right)}D\left( \theta \right) d\theta +\int_0^{\psi \left( \frac{\theta _s-\theta _r}{2} \right)}{\left( \theta _s+\theta \left( \psi \right) -2\theta _0 \right)}K\left( \psi \right) d\psi 
 \end{equation}$$

where $D(θ)$ is the diffusivity function with $D\left( \theta \right) =K\left( \theta \right) \frac{d\psi}{d\theta}$. This equation splits the integral 

$$\begin{equation}
Sorpt^2\left( \theta _0,\theta _s \right) =\int_{\theta _r}^{\theta _s}{\left( \theta _s+\theta -2\theta _0 \right)}D\left( \theta \right) d\theta
\end{equation}$$

into two parts to allow the integration of continuous functions over closed intervals. In the regular expression, the diffusivity function is infinite close to saturation, $θ \to θ_{s}$, which complexes its integration in the vicinity of $θ_{s}$. In the specific equation, the last part of the integration based on $ψ$ is replaced with the integration of the hydraulic conductivity as a function of the water pressure, $K(ψ)$, alleviating the problem of convergence. 

Also, $\varDelta Qmax_{\frac{1}{2}}^{t}$, which is the maximum water infiltrating into the top cell, must be less than the maximum available pore volume of the top cell:

$$\begin{equation}
\begin{cases}	\varDelta Sink_1\left( \psi _{1}^{t-1} \right) =\varDelta Z_i\,\,\varDelta T^t\,\,Sink_1\left( \psi _{1}^{t-1} \right)\\	\varDelta Qmax_{\frac{1}{2}}^{t}=\,\,Min\left[ \varDelta Qmax_{\frac{1}{2}}^{t}; \varDelta Z_1\left[ \theta s_1-\theta _{1}^{t-1} \right] +\varDelta Sink_1\left( \psi _{1}^{t-1} \right) \right]\\\end{cases}
\end{equation}$$

where $θs_{1}$ [L3 L-3] is the saturated volumetric soil water content and $\varDelta Sink_{1}$ [L] is the sink of the top cell.

The amount of water that is not able to infiltrate into the soil can either run off laterally due to the slope or get ponded at the surface, where the ponding $\varDelta H_{pond}^{t}$, [L] is computed as:

$$\begin{equation}
\varDelta H_{pond}^{t}=Max\left\{ \varDelta Pr_{through}^{^t}+\varDelta H_{pond}^{t-1}-\varDelta Qmax_{\frac{1}{2}}^{t};\,\,0 \right\} 
\end{equation}$$

where $\varDelta Pr_{through}^{t}$ [L] is the *throughfall precipitation* (i.e., the amount of water reaching the top cell; computed in [rainfall interception](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYPIX/Interception)).

