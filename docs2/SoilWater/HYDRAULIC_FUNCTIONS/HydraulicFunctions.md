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


# Soil hydraulic functions 

## Bimodal Kosugi hydraulic functions

 The representation of the $θ(ψ)$ and $K(θ)$ functions is based on the dual porosity models of [Pollacco *et al*. (2017)](#_ENDREF_40), and we use the bimodal lognormal expressions of the Kosugi model from [Fernández-Gálvez *et al*. (2021)](#_ENDREF_41): 

$$\begin{equation}
\begin{cases}	\theta \left( \psi \right) =\theta _{\mathrm{Mat}}\left( \psi \right) +\theta _{\mathrm{Mac}}\left( \psi \right)\\	\theta _{\mathrm{Mat}}\left( \psi \right) =\frac{1}{2}\left( \theta _{s\mathrm{MacMat}}-\theta _r \right) \,\,erfc\left( \frac{\ln \psi  -\ln  \psi _m}{\sigma \sqrt{2} } \right) +\theta _r\\	\theta _{\mathrm{Mac}}\left( \psi \right) =\frac{1}{2}\left( \theta _s-\theta _{s\mathrm{MacMat}} \right) \,\,erfc\left( \frac{\ln \psi  -\ln  \psi _{m\mathrm{Mac}}}{\sigma _{\mathrm{Mac}}\sqrt{2} } \right)\\\end{cases}
\end{equation}$$

$$\begin{equation}
\begin{cases}	K\left( \psi \right) =K_{\mathrm{Mat}}\left( S_{e\,\,\mathrm{Mat}}\left( \psi \right) \right) +K_{\mathrm{Mac}}\left( S_{e\,\,\mathrm{Mac}}\left( \psi \right) \right)\\	S_e=\frac{\theta -\theta _r}{\theta _s-\theta _r}\\	K_{\mathrm{Mat}}\left( \psi \right) =K_s\frac{\theta _{s\mathrm{MacMat}}-\theta _{\mathrm{r}}}{\theta _s-\theta _r}\,\,\sqrt{S_{e\,\,\mathrm{Mat}}\left( \psi \right)}\,\,\left[ \frac{1}{2}erfc\left( erfc^{-1}\left( 2S_{e\,\,\mathrm{Mat}}\left( \psi \right) \right) +\frac{\sigma}{\sqrt{2}} \right) \right] ^2\\	K_{\mathrm{Mac}}\left( \psi \right) =K_s\frac{\theta _{\mathrm{s}}-\theta _{s\mathrm{MacMat}}}{\theta _s-\theta _r}\,\,\sqrt{S_{e\,\,\mathrm{Mac}}\left( \psi \right)}\,\,\left[ \frac{1}{2}erfc\left( erfc^{-1}\left( 2S_{e\,\,\mathrm{Mac}}\left( \psi \right) \right) +\frac{\sigma _{\mathrm{Mac}}}{\sqrt{2}} \right) \right] ^2\\\end{cases}
\end{equation}$$

where $erfc$ is the complementary error function; $θ$ [L<sup>3</sup> L<sup>-3</sup>] represents the volumetric soil water content and $ψ$ [L] the soil water pressure, considering $ψ > 0$ for unsaturated soils (i.e. matric suction); $θ_s$ [L<sup>3</sup> L<sup>-3</sup>] and $θ_r$ [L<sup>3</sup> L<sup>-3</sup>] are the saturated and residual volumetric soil water content, respectively; ln $ψ_m$ [L] (with the argument of ln in units of length, i.e., $ψ_m$ in [L]) and $σ$ [-] denote the mean and standard deviation of ln $ψ$ [L], respectively in the soil matrix domain; ln $ψ_{mMac}$ [L] and $σ_{Mac}$ [-] denote the mean and standard deviation of ln $ψ$, respectively, in the macropore soil domain; $θ_{sMacMat}$ [L<sup>3</sup> L<sup>-3</sup>] is the volumetric saturated water content that theoretically differentiates inter-aggregate pores (structural macropores) and matrix domains (intra-aggregate micropores), defining the corresponding soil water pressure threshold between macropore and matrix $ψ_{MacMat}$ [L]; $K_s$ [L T<sup>-1</sup>] is the saturated hydraulic conductivity; and $S_{e\ Mat}$ [-] and $S_{e\ Mac}$ [-] denote the effective saturation as a function of $ψ$ in the soil matric and macropore domains, with values between 0 and 1. Finally, $K(ψ)$ refers to the unsaturated hydraulic conductivity, written as a function of $ψ$.

## Constrained soil bimodal Kosugi hydraulic parameters

Estimating the bimodal Kosugi soil hydraulic parameters from observed volumetric soil water content requires the simultaneous estimation of eight parameters ($θ_s$, $θ_r$, $σ$, $ψ_m$, $K_s$, $θ_{sMacMat}$, $σ_{Mac}$, and $ψ_{mMac}$), using limited measurement data. Here, we use the term ‘optimize’ in the estimation of the parameters since the task is to produce a set of values for the parameters that optimizes (typically, minimizes) some objective function. Here, optimization is used in the sense of parameter estimation.

In certain cases where, for example, the measurements have a restricted range [(Pollacco *et al*., 2008b)](#_ENDREF_42), full inversion of the data to produce parametric estimates of the bimodal Kosugi soil hydraulic model is not possible, leading to highly sensitive parameters (eg. [Pollacco *et al*., 2008b](#_ENDREF_43), [2013](#_ENDREF_44); [Pollacco and Mohanty, 2012](#_ENDREF_45)). Reducing the sensitivity of the parameters requires either adding novel independent measurements or incorporating constraints to reduce the effective complexity of the model. Since gathering new measurements is, in this case, not feasible, here we use the set of constraints proposed by [Fernández-Gálvez *et al.* (2021)](#_ENDREF_41), which reduces the number of parameters to be optimized without compromising the fit of the hydraulic functions, while the estimated hydraulic parameters still have physical meaning. 

This set of constraints can be summarized as follows: $θ_r$ is derived from $σ$, $ψ_{mMac}$ and $σ_{Mac}$, which are considered constant from a fixed value of $ψ_{MacMat}$, equal to 100 mm, and $ψ_m$ and $σ$ are dynamically constrained, based on the assumption that $θ(ψ)$ and $K(θ)$ are lognormally distributed. Therefore, the number of hydraulic parameters to be optimized is reduced from eight to five using the principles of soil physics. Although $θ_s$ can be derived from total porosity when soil bulk density and particle density data are available, in this case it is an optimized parameter with a feasible range determined from the maximum observed $θ$ in the corresponding layer. The soil hydraulic parameters are optimized using the bimodal $θ(ψ)$ and $K(θ)$ model described this [section](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYDRAULIC_FUNCTIONS/HydraulicFunctions). This is performed by matching the simulated $θ$ time series from the HyPix hydrological model with the observed $θ$ time series at the corresponding depth. 

The dynamic, physically feasible range of the five optimized bimodal Kosugi hydraulic parameters using the unconstrained conditions and the full set of constraints derived from [Fernández-Gálvez *et al*. (2021)](#_ENDREF_41) is indicated in [Table 1](). In both cases, $ψ_{MacMat}$ is constant, with a value of 100 mm and $P_σ = 3$. Here, $θ_r$ is derived from $σ$ in [Eq. (3)]() and [Table 1](), and is given by:

$$\begin{equation}
\begin{cases}	\sigma^{*}=\frac{\sigma -\sigma_{min}}{\sigma_{max}-\sigma_{min}}\\	\theta_{r}\left( \sigma \right) =\,\,\theta_{rMax}\,\,\frac{1-e^{-\alpha_1\cdot \sigma^{*\,^{\alpha_2}}}}{1-e^{-\alpha_1}}\\\end{cases}
\end{equation}$$



where $θ_{rMax}$ is set at 0.20, the maximum value for $θ_r$ that was found to be satisfactory; $α_1 = 15$ and $α_2 = 4$ are two optimized empirical parameters; $σ^*$ [-] is the normalized $σ$; and $σ_{min}$ [-] and $σ_{max}$ [-] are set at 0.75 and 4.00 from [Table 1](). Therefore, with all these simplifications and additional constraints, we define a model that requires only five parameters:  $θ_s$, $σ$, $ψ_m$, $K_s$, and $θ_{sMacMat}$.

 <figcaption align = "center"><b>Table 1 - Feasible dynamic range of the optimized bimodal Kosugi hydraulic parameters from observed θ. Both the unconstrained and constrained sets of hydraulic parameters have five parameters to be optimized. The difference is that in the dynamically constrained set of hydraulic parameters there is a relationship between σ and ψm (Fernández-Gálvez et al., 2021). The feasible range of Ks is derived from Carsel and Parrish (1988). For both the unconstrained and constrained sets of hydraulic parameters, ψMacMat = 100 mm, and Pσ = 3. θr(σ) is described in [Eq. (3)].</b></figcaption></figure>

 $\\$ | $\boldsymbol{\theta_{s}}$ [m <sup>3</sup>m<sup>-3</sup>] | $\boldsymbol{\theta_{r}}$ [m <sup>3</sup>m<sup>-3</sup>]| $\boldsymbol{σ}$ [-]| $\boldsymbol{\psi _{m}}$ [mm]| $\boldsymbol{K_{s}}$ [cm h<sup>-1</sup>] | $\boldsymbol{\theta_{sMacMat}}$ [m <sup>3</sup>m<sup>-3</sup>]| $\boldsymbol{\psi _{mMac}}$ [mm] | $\boldsymbol{σ_{Mac}}$ [-]
--|---|--------------------|-----|----|-----------|-----------|----|---
$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |
**Min** | $Max(\theta)$ | $\theta_{r}(σ)$ | 0.75 | $\sqrt{\psi_{MacMat}}e^{σ_{min}P_{σ}}$ | 0.02 | $0.75\theta_s$ | $\sqrt{\psi_{MacMat}}$ | $ln \psi_{MacMat}/2P_{σ}$
**Max** | 0.65 | $\theta_{r}(σ)$ | 4.00 | $\sqrt{\psi_{MacMat}}e^{σ_{max}P_{σ}}$ | 30.00 | $\theta_s$ | $\sqrt{\psi_{MacMat}}$ |$ln \psi_{MacMat}/2P_{σ}$
$\\$ |$\\$ |$\\$ |$\\ $ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |
**Min** | $Max(\theta)$ | $\theta_{r}(σ)$ | 0.75 | $\sqrt{\psi_{MacMat}}e^{σP_{σ}}$ | 0.02 | $0.75\theta_s$ | $\sqrt{\psi_{MacMat}}$ | $ln \psi_{MacMat}/2P_{σ}$
**Max** | 0.65 | $\theta_{r}(σ)$ | 4.00 | $\sqrt{\psi_{MacMat}}e^{σP_{σ}}$ | 30.00 | $\theta_s$ | $\sqrt{\psi_{MacMat}}$ | $ln \psi_{MacMat}/2P_{σ}$


## Bimodal van Genuchten hydraulic functions

## Brorks and Corey hydraulic functions

## Clapp and Hornberger hydraulic functions