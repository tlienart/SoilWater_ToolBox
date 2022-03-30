# HYPIX

Modelling unsaturated flow in highly heterogeneous soils can be accurately performed by solving the [Richards (1931)](#_ENDREF_13) equation (RE), which is commonly adopted by soil vegetation atmosphere transfer models. However, RE is highly nonlinear, and despite numerous efforts over the last decade its solution using numerical methods is demanding, and problems finding techniques to achieve fast and accurate solutions are unresolved (e.g. [Zha *et al*., 2019](#_ENDREF_14)). Here, we propose improvements to RE and implement them in HyPix, using the mixed form of RE, as recommended by [Celia *et al*. (1990)](#_ENDREF_14). The solution of RE is based on [Hassane Maina and Ackerer (2017)](#_ENDREF_15), for which the RE partial differential equation is solved using a *cell-centered finite-volume (implicit finite differences)* scheme for the spatial discretization, with an implicit Euler scheme for the temporal discretization by using the weighted average inter-cell hydraulic conductivity.

##	*Richards equation of HyPix model*

Assuming a rigid solid matrix ([Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp)), the mixed form of RE is written as:

$$\begin{equation}
\frac{\theta _i\left( \psi _{i}^{t} \right) -\theta _i\left( \psi _{i}^{t-1} \right)}{\varDelta T^t}-S_o\frac{\theta _i\left( \psi _{i}^{t} \right)}{\theta _{s_i}}\frac{\left| \psi _{i}^{t} \right|-\left| \psi _{i}^{t-1} \right|}{\varDelta T^t}=\frac{Q_{i-\frac{1}{2}}^{t}-Q_{i+\frac{1}{2}}^{t}}{\varDelta Z_i}-Sink_i\left( \psi _{i}^{t-1} \right)
\end{equation}$$

where $\varDelta T^t$ [T] is the time-step at time $t$; $\varDelta Z_{i}$ [L] is the mesh size of the cell $i$, with the vertical coordinate positive downwards; $θ_{i}$ [L<sup>3</sup> L<sup>-3</sup>] is the volumetric soil water content of the cell $i$; $θ_{s}$ [L<sup>3</sup> L<sup>-3</sup>]  is the saturated volumetric soil water content; $S_{0}$ [L<sup>-1</sup>] is a parameter that accounts for fluid compressibility, which is assumed to be constant with depth; $ψ_{i}$ [L] is the soil water pressure of cell $i$, considering $ψ > 0$ for unsaturated soils; $Q$ [L T<sup>-1</sup>] is the soil water flux based on the extended Darcy’s law, positive downward and negative when water moves upwards; and $Sink_{i}$ [L<sup>3</sup> L<sup>-3</sup>], taken as positive, is the sink term defined as the volume of water per unit time removed from cell $i$ by *evaporation* and *root water uptake*.

![HyPix](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp "Figure 1. Diagram describing the 1D vertical discretization of the Richards equation, where i  is the cell number (cell 1 is the top cell and cell Ni is the bottom cell, therefore cell = i  is below cell = i¬–1); ΔPr  [L] is the precipitation reaching the top of the canopy; ΔPrground [L]  is the precipitation reaching the soil surface (cell = 1); ΔHpond  [L] is the ponding water; and ΔQi+1/2 = Qi+1/2 ΔT [L] is the inter-cell water volume (positive downwards). Water is removed from the soil profile by transpiration, ΔTransp [L], and evaporation, ΔEvapo [L], depending on θ and potential evapotranspiration, ΔPet [L] (partitioned between potential evaporation, ΔPetevap [L], and potential transpiration, ΔPettransp [L]).")
<figcaption align = "center"><b>Figure 1 - Diagram describing the 1D vertical discretization of the Richards equation, where i  is the cell number (cell 1 is the top cell and cell Ni is the bottom cell, therefore cell = i  is below cell = i–1); ΔPr  [L] is the precipitation reaching the top of the canopy; ΔPrground [L]  is the precipitation reaching the soil surface (cell = 1); ΔHpond  [L] is the ponding water; and ΔQi+1/2 = Qi+1/2 ΔT [L] is the inter-cell water volume (positive downwards). Water is removed from the soil profile by transpiration, ΔTransp [L], and evaporation, ΔEvapo [L], depending on θ and potential evapotranspiration, ΔPet [L] (partitioned between potential evaporation, ΔPetevap [L], and potential transpiration, ΔPettransp [L]).</b></figcaption></figure>


### **Computing water infiltration into the soil profile**

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

where $θs_{1}$ [L<sup>3</sup> L<sup>-3</sup>] is the saturated volumetric soil water content and $\varDelta Sink_{1}$ [L] is the sink of the top cell.

The amount of water that is not able to infiltrate into the soil can either run off laterally due to the slope or get ponded at the surface, where the ponding $\varDelta H_{pond}^{t}$, [L] is computed as:

$$\begin{equation}
\varDelta H_{pond}^{t}=Max\left\{ \varDelta Pr_{through}^{^t}+\varDelta H_{pond}^{t-1}-\varDelta Qmax_{\frac{1}{2}}^{t};\,\,0 \right\}
\end{equation}$$

where $\varDelta Pr_{through}^{t}$ [L] is the *throughfall precipitation* (i.e., the amount of water reaching the top cell; computed in [rainfall interception](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYPIX/Interception)).

###	**Computing water fluxes, Q**

   #### *Cells below the top cell: $2 ≤ i ≤ N_{i}$*

   >>The Darcian fluid flux density, $Q_{i-\frac{1}{2}}^{t}$ [L T<sup>-1</sup>], is computed by using the inter-cell hydraulic conductivity between cell $i$ and $i – 1$ ([Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp)), as follows:

   $$\begin{equation}
   \begin{cases}	Q_{i-\frac{1}{2}}^{t}=-K_{i-\frac{1}{2}}\,\,\left[ \frac{\left| \psi _{i-1}^{t} \right|-\left| \psi _{i}^{t} \right|}{\varDelta Z_{i-\frac{1}{2}}}-\cos  \alpha \right]\\	\varDelta Z_{i-\frac{1}{2}}=\frac{\varDelta Z_i+\varDelta Z_{i-1}}{2}\\	K_{i-\frac{1}{2}}=\omega _i\,\,K_i\left( \psi _{i}^{t} \right) +\left[ 1-\omega _i \right] \,\,K_{i-1}\left( \psi _{i-1}^{t} \right)\\	\omega _i=\frac{\varDelta Z_i}{\varDelta Z_i+\varDelta Z_{i-1}}\\\end{cases}
   \end{equation}$$

   >>where $Q_{i-\frac{1}{2}}^{t}$ [L T<sup>-1</sup>] is the flux entering cell $i$ from the top, and $ Q_{i+\frac{1}{2}}^{t}$ [L T<sup>-1</sup>] is the flux exiting cell $i$ from the bottom; $K_{i-\frac{1}{2}}$ [L T<sup>-1</sup>] refers to the weighted average inter-cell hydraulic conductivity (e.g. [Haverkamp and Vauclin, 1979](#_ENDREF_21); [Belfort *et al*., 2013](#_ENDREF_22)), computed with $\omega _i$; and $\varDelta Z_{i-\frac{1}{2}}$ [L] is the distance between cell centres $i$ and $i – 1$, as described in [Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp).


#### *Free drainage bottom condition: $i = N_{i} + 1$*

>>The water flux leaving the bottom cell represents the drainage, and it can be described as a function of the hydraulic conductivity of the bottom cell as:

$$\begin{equation}
Q_{N_{i}+1}^{t}=K_{N_{i}}\left( \psi _{N_{i}}^{t} \right) \,\,\cos \alpha
\end{equation}$$

>>where $K_{N_{i}}\left( \psi_{N_{i}}^{t} \right)$ is the [*unsaturated hydraulic conductivity*](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYDRAULIC_FUNCTIONS/HydraulicFunctions)

>>This free drainage condition is the most widely used when the water table is assumed to be at significant depth. Typically, our model addresses the modelling of water flow in the unsaturated or vadose zone.


##	*Solving the Richards equation using the Newton–Raphson method*

The Picard and Newton–Raphson (NR) iterative methods are the most widely used procedures for solving the RE. NR has been found to be more efficient than the Picard iteration method ([Lehmann and Ackerer, 1998](#_ENDREF_23); [Paniconi and Putti, 1994a](#_ENDREF_24), [1994b](#_ENDREF_25)), so we solved the RE by using the NR algorithm. NR is computed using a first-order Taylor development by solving the Jacobian matrices of the residuals, R [Eq. (12)] (described in Appendix 7.1) in an iterative way that updates $\left[ \psi _i^{t,k+1}-\psi _{i}^{t,k} \right] $ until convergence is achieved. The numerical discretization is a tridiagonal, nonlinear set of equations that needs to be solved for $\left[ \psi _i^{t,k+1}-\psi _{i}^{t,k} \right] $, and for every iteration, $k$:

$$\begin{equation}
\begin{cases}
\psi _i^{t,k=1}=\psi _{i}^{t-1}\\
\psi _i^{t,k+1}=\psi _i^{t,k}-\frac{R\left( \psi _{i}^{t,k} \right)}{\frac{\partial R\left( \psi _i^{t,k} \right)}{\partial \psi _{i-1}^{t,k}}+\frac{\partial R\left( \psi _i^{t,k} \right)}{\partial \psi _{i}^{t,k}}+\frac{\partial R\left( \psi _i^{t,k} \right)}{\partial \psi _{i+1}^{t,k}}}\\
\psi _i^{t,k+1}=\varOmega \,\,\min \left[ \max \left( \psi _i^{t,k+1},0 \right) ,\psi _{\max} \right] +\left( 1-\varOmega \right) \psi _i^{t,k}\\
\end{cases}
\end{equation}$$

where $\psi_{\max} = 10^{7}$ mm is the maximum value of $\psi $; and $Ω [1 ; 0[$ is a parameter used to avoid "overshooting" of the Newton step. Therefore, when $Ω = 1$ there is no reduction of the Newton step. We use $Ω = 0.5$, as recommended by ([Kelley, 2003](#_ENDREF_26)). The feasible range of Ω is set to $[0.2 ; 1.0]$ as indicated in [Table 2](XXXX). In this study the initial $\psi _{i}^{t=0}$ is derived from measured $θ^{t=0}$.

It is expected that for every $k$, $R$ decreases such that $\left| R\left( \psi _{i}^{t,k+1} \right) \right|\,\,\leqslant \left| R\left( \psi _{i}^{t,k} \right) \right|$. The solution of the Jacobian takes place by means of the *tridiagonal function* that is an effective modification of the Gauss algorithm for the solution of the tridiagonal linear set of equations ([Appendix 7.1](XXXXX)). We use the efficient *tridiagonal function* derived from the LAPACK package within Julia to solve the matrix (https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/).

### **Residuals of the Richards equation**

The NR method is used to solve $\psi _{i}^{t}$ such that the residuals of each cell, $R_{i}^{t}$ [L], equal close to 0:

$$\begin{equation}
\begin{cases}	R_{i}^{t}=\varDelta Z_i\left[ \theta _i\left( \psi _{i}^{t} \right) -\theta _i\left( \psi _{i}^{t-1} \right) \right] -\varDelta Z_i\,\,S_o\frac{\theta _i\left( \psi _{i}^{t} \right)}{\theta _{s_i}}\left( \left| \psi _{i}^{t} \right|-\left| \psi _{i}^{t-1} \right| \right) -\varDelta T^t\left[ Q_{i-\frac{1}{2}}^{t}-Q_{i+\frac{1}{2}}^{t} \right] +\varDelta Sink_i\left( \psi _{i}^{t-1} \right)\\	\varDelta Sink_i\left( \psi _{i}^{t-1} \right) =\varDelta Z_i\,\,\varDelta T^t\,\,Sink_i\left( \psi _{i}^{t-1} \right)\\\end{cases}
\end{equation}$$

where $\varDelta Sink_{i}$ [L] is the sink computed for a given $\varDelta T$. To stabilize the numerical scheme, $\varDelta Sink_{i}$ is computed with the $ψ_{i}$ values derived from the previous time-step.

### **Automatic differentiation of the Jacobian with Julia language**

One of the shortcomings of the NR solver is that it requires the mathematical derivatives of $R$ [Eq. (9)]() (Appendix 7.2),  the implementation of which is complicated and time consuming. For example, if we wish to modify the inter-cell hydraulic conductivity [Eq. (7)](), it requires recalculation of the derivatives.
To address this shortcoming, HyPix implements an option whereby the derivatives are analytically derived automatically by using the forward-mode automatic differentiation *[ForwardDiff](https://github.com/JuliaDiff)* in the Julia package ([Revels *et al*., 2016](#_ENDREF_27)). *ForwardDiff* (version 0.10.16) was found to be as accurate as using the mathematical derivatives, and only 10–25% slower compared to using mathematical derivatives.


### **Convergence criterion**

For a given time-step, $\varDelta T$, the iteration of the NR stops either when $k = N_{k}$, where $N_{k}$ is a user-defined maximum number of iterations ([Table 2](XXXX)), or when a convergence criterion is satisfied for the overall water balance of the residuals, $WB_{residual}$ [T<sup>-1</sup>], as described below.

TABLE 2
<figcaption align = "center"><b>Table 2 - Minimum, maximum, and universal values of Richards equation parameters.</b></figcaption></figure>


 $\\$ | $\boldsymbol{\varDelta T_{min}}$ [s] | $\boldsymbol{\varDelta T_{max}}$ [s] | $\boldsymbol{P_{\varDelta θ_{max} Rerun}}$ [-] | $\boldsymbol{N_{k}}$ [-] | $\boldsymbol{WB_{residual}}$ [T <sup>-1</sup>] | $\boldsymbol{\varDelta θ_{max}}$ [L <sup>3</sup>L <sup>-3</sup>] | $\boldsymbol{Ω}$ [-]
--|---|--------------------|-----|----|-----------|-----------|----
**Min** | 1 | $\varDelta T_{min}$ | 1.0 | 10 | $10^{-8}$ | $10^{-2}$ | 0.2
**Max** | $\varDelta T_{max}$ | 3600 | $\varDelta T_{max}/ \varDelta T_{min}$ | 50 | $10^{-20}$ | $10^{-6}$ | 1.0
**Value** | **1** | **3600** | **1.2** | **30** | $\boldsymbol{10^{-11}}$ | $\boldsymbol{10^{-3.3}}$ | **0.5**

We derive a convergence criterion for the NR method based on the Euclidean norm of a vector $R$ ([Driscoll and Braun, 2017](#_ENDREF_28); [Kochenderfer and Wheeler, 2019](#_ENDREF_29); [Kelley, 2003](#_ENDREF_26)), where the iteration stops when the following criterion for $WB_{residual}$ is met:

$$\begin{equation}
\sqrt{\frac{\sum_{i=1}^{N_{\mathrm{i}}}{\left[ \frac{R_{i}^{t}}{\varDelta T^t\,\,\varDelta Z_i} \right] ^2}}{N_{\mathrm{i}}}}\leqslant WB_{\mathrm{residual}}
\end{equation}$$

where Ni is the total number of cells, as described in [Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp), and $WB_{residual}$ is the overall water balance of the residuals. The $R$ is normalized, such that WBresidual is independent of the cell size and the time-step. $WB_{residual}$ considers the overall error of the water balance of the soil profile. The feasible range of $WB_{residual}$ is provided in [Table 2](XXXX).

### **Novel adaptive time-step management**

The time management module optimizes the size of the time-step, $\varDelta T$, such that HyPix uses the largest $\varDelta T$ while meeting the targeted water balance accuracy. There are a large number of heuristic time-stepping methods in the literature (e.g. [Thomas and Gladwell, 1988](#_ENDREF_30); [Kavetski *et al.*, 2001](#_ENDREF_31); [Kavetski and Binning, 2004](#_ENDREF_32); [Miller *et al.*, 2006](#_ENDREF_33)). Among them we selected and improved on the adaptive time-step management of Kirkland et al. (1992) and Ross (2003), because their method is physically based, such that $\varDelta T$ is directly derived from the residuals of the water balance [Eq. (10)]() and requires only one physical fitting physical parameter. $\varDelta T$ is calculated via a maximum increase or decrease of the degree of saturation for each cell, which ensures a higher time resolution when $θ$ variations are large. This improves the convergence rate at the wetting front.

#### **Traditional adaptive time-step management: $\varDelta T$-$\varDelta θ$**

The time-step management of [Kirkland *et al.* (1992)](#_ENDREF_34) and [Ross (2003)](#_ENDREF_35) was developed specifically for the RE based on $θ$. $\varDelta T$-$Δθ$ is derived by rearranging the terms of the residual [Eq. (10)]() assuming that $R ≈ 0$ and $S_{0} ≈ 0$ ($S_{0}$ is by definition small and strictly 0 for non-compressible fluids), leading to:

$$\begin{equation}
\begin{cases}	\varDelta T_i=\frac{\varDelta Z_i\,\,\varDelta \theta _{max_i}+\varDelta Sink_i\left( \psi _{i}^{t-1} \right)}{\left[ Q_{i-\frac{1}{2}}^{t}-Q_{i+\frac{1}{2}}^{t} \right] +\epsilon}\\	\varDelta \theta _{max_i}=\left| \theta _i\left( \psi _{i}^{t} \right) -\theta _i\left( \psi _{i}^{t-1} \right) \right|\\\end{cases}
\end{equation}$$

where $\epsilon$ is a small, positive, infinitesimal quantity used to avoid dividing by 0, and $\varDeltaθ_{max}$ [L<sup>3</sup> L<sup>-3</sup>] is a parameter describing the maximum change ‘allowed’ of $θ$ for a given $\varDelta T$. This is to assure numerical stability as well as to avoid oscillation in the solution. The $\Delta θ_{max}$ feasible range is provided in [Table 2](). The selected $\varDelta T$ is computed by using the Euclidean norm of a vector based on [Driscoll and Braun (2017)](#_ENDREF_36), [Kochenderfer and Wheeler (2019)](#_ENDREF_37), and [Kelley, (2003)](#_ENDREF_38) by using an algorithm similar to [Eq. (11)]():

$$\begin{equation}
\begin{cases}	\varDelta T^t=\sqrt{\frac{\sum_{i=1}^{N_i}{\varDelta T_i^2}}{N_i}}\\	\varDelta T_{min}\leqslant \varDelta T^t\leqslant \varDelta T_{max}\\\end{cases}
\end{equation}$$

where $N_i$ is the total number of cells, as described in [Figure 1](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure1.bmp)), and $\Delta T_{min}$ [T] and $\Delta T_{max}$ [T] are the user-defined minimum and maximum time-steps described in [Table 2]().

The computational time in HyPix decreases when $\varDelta T$ increases when, for example: **(a)** the size of the cell, $\varDelta Z$, increases, **(b)** the hydraulic properties from one cell to another, depicted by $Q$, do not change dramatically (homogeneous soils), and **(c)** the amplitude of the variation of $θ$ decreases.

#### **Novel adaptive time-step management: $\varDelta T$-$\varDelta \psi$**

For the robustness of the solution of the RE [(Eq. (1))]() based on $ψ$ , the change of $\varDelta \psi$ between two consecutive time-steps needs to be small. Nevertheless, as shown in [Figure 2](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure2.bmp), a small change of $θ$ could result in a large change of $\psi$, particularly near saturation and at the dry end of the $θ(ψ)$ curve. To address this shortcoming, we improve the traditional [Ross (2003)](#_ENDREF_39) time-step, $\varDelta T-\varDeltaθ$ [(Eq. (12))](), by allowing $\varDelta θ_{max}$ to vary $\varDelta θ_{max_i}$. For that we introduce a temporary parameter, $\varDelta ψ_{max}$, computed from $θ(ψ)$ and $\varDelta θ_{max}$ [(Eq. (12))]():

$$\begin{equation}
\varDelta \psi _{max_i}\geqslant \left| \psi _{i}^{t}-\psi _{i}^{t-1} \right|
\end{equation}$$

To avoid introducing a second parameter we derive $\varDelta ψ_{max}$ from $\varDelta θ_{max}$. As shown in [Figure 2](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure2.bmp), for a given $θ(ψ)$, the smallest $\varDelta ψ$ occurs around $(θ_s-θ_r)/2$ or at the Kosugi parameter $ψ = ψ_m$. Therefore, $\varDelta ψ_{max}$ is computed for every cell as follows:

$$\begin{equation}
\begin{cases}	\theta _{\frac{1}{2}}=\frac{\theta _{r_i}+\theta _{s_{MacMat_i}}}{2}\\	\varDelta \psi_{\max_i}=\psi \left( \theta _{\frac{1}{2}}-\frac{\varDelta \theta _{max_i}}{2} \right) -\psi \left( \theta _{\frac{1}{2}}+\frac{\varDelta \theta _{max_i}}{2} \right)\\\end{cases}\end{equation}$$

The $\varDelta θ_{max}$ feasible range is provided in [Table 2](), and $\varDelta θ_{max}(ψ)$ is adjusted based on the maximum allowed change of $\varDelta ψ_{max}$, and $ψ$:

$$\begin{equation}
\varDelta \theta _{\max _i}=\frac{\theta _i\left( \max \left\{ \psi _{i}^{t}-\varDelta \psi _{\max _i},0 \right\} \right) -\theta _i\left( \psi _{i}^{t}+\varDelta \psi _{\max _i} \right)}{2}
\end{equation}$$

$\varDelta T_i$ is computed using [Eq. (12)]() but using $\varDelta \theta_{max_i}$ instead of $\varDelta θ_{max}$:

$$\begin{equation}
\varDelta T_i=\frac{\varDelta Z_i\,\,\varDelta \theta_{max _i}+\varDelta Sink_i\left( \psi _{i}^{t-1} \right)}{\left[ Q_{i-\frac{1}{2}}^{t}-Q_{i+\frac{1}{2}}^{t} \right] +\epsilon}
\end{equation}$$

The selected $\varDelta T$ is computed using [Eq. (13)]().

![HyPix](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/Figure2.bmp "Figure 2. Example illustrating that small changes of θ could result in large changes of ψ, particularly near saturation and at the dry end of the unimodal θ(ψ) [Eq. (1)].")
<figcaption align = "center"><b>Figure 2 - Example illustrating that small changes of θ could result in large changes of ψ, particularly near saturation and at the dry end of the unimodal θ(ψ) [Eq. (1) FROM HYDRAULIC FUNCTIONS!!!!].</b></figcaption></figure>

## **Condition to rerun the time-step**
$\varDelta T^t$ is computed by using past pressure, $ψ^{t-1}$, which may not reflect, for example, the passage of a wetting front. Therefore, to assure a good water balance and the stability of the solution of the RE, HyPix recomputes $\varDelta T^t$ if $\varDelta T^t(\psi ^t)\ll \varDelta T^t(\psi ^{t-1})$. HyPix is rerun if the following condition is met:

$$\begin{equation}
\frac{\varDelta T_{}^{t}\left( \psi _{}^{t-1} \right)}{\varDelta T_{}^{t}\left( \psi _{}^{t} \right)}>P_{\varDelta \theta _{max}\_Rerun}
\end{equation}$$

where $P_{\varDelta θ _{max}\_Rerun}$ is a user-defined parameter for which the feasible range is provided in [Table 2]().


### **Accuracy and efficiency of simulations**

#### **Water balance**

The overall water balance, $WB$ [L], of the simulation is derived from the residuals, $R$ [[Eq. (10)]](), and is computed for every time-step as follows:

$$\begin{equation}
WB_{}^{t}=\sum_{i=1}^{N_i}{\varDelta Z_i\left[ \theta _i\left( \psi _{i}^{t} \right) -\theta _i\left( \psi _{i}^{t=0} \right) \right]}-\sum_{i=1}^{N_i}{\sum_{t=1}^{N_t}{\varDelta Z_i\,\,S_o\frac{\theta _i\left( \psi _{i}^{t} \right)}{\theta s_i}\left( \left| \psi _{i}^{t} \right|-\left| \psi _{i}^{t-1} \right| \right)}}\\-\sum_{t=1}^{N_t}{\varDelta T^t\left[ Q_{\frac{1}{2}}^{t}-Q_{N_{i+\frac{1}{2}}}^{t} \right] +\sum_{i=1}^{N_i}{\sum_{t=1}^{N_t}{\varDelta Sink_i\left( \psi _{i}^{t} \right)}}}
\end{equation}$$

where $N_t$ is the final time-step and $N_i$ is the bottom cell.
Because $WB$ increases with the length of the simulation, it is normalized, $WB^*$, to the cumulative infiltration:

$$\begin{equation}
WB^{*t}=\frac{WB_{}^{t}}{\sum_{t=1}^{Nt}{\varDelta Q_{\frac{1}{2}}^{t}}}
\end{equation}$$

An acceptable water balance at the end of the simulation occurs when $WB^{*t}$ is smaller than the uncertainty of measuring precipitation.

#### **Efficiency of solving the Richards equation**

The efficiency of a simulation, $E_{ff}$ [T <sup>-1</sup>], is defined as the average number of iterations, $k$, per day required to achieve a suitable $WB^*$ [[Eq.(20)]](), and it is computed as:

$$\begin{equation}
E_{ff}=\frac{N_{iter}}{\sum_{t=1}^{N_{\mathrm{t}}}{\varDelta T^t}}\,\,
\end{equation}$$

where $N_{iter}$ is the number of iterations. Therefore, the smaller the $E_{ff}$, the faster HyPix would run for a given $WB^*$.
