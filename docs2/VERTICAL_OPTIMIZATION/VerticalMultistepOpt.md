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


# Vertical Multistep Optimization

This section describes the strategy to optimize the hydraulic parameters of the parsimonious bimodal [θ(ψ) and K(θ)](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYDRAULIC_FUNCTIONS/HydraulicFunctions) for each layer of the soil profile by using the [HyPix]() model. Each layer of the profile corresponds to the depths where $θ$ is experimentally measured. 

It is challenging to optimize a large number of parameters simultaneously, and [Pollacco (2005)](#_ENDREF_46) found that optimizing each layer separately, one by one, produces poor results, particularly for a highly heterogeneous soil profile, in the presence of lenses of clay, pebbles, or macropore flow. This is because an optimization method that isolates each layer without considering the overall water flow in the soil profile causes unrepresentative parameters, and thus poor representation of the soil water fluxes. Therefore, we present an inverse modelling algorithm for layered soil. 

When optimizing the profile soil hydraulic parameters, initially the soil profile is considered to be homogeneous. Then a stepwise grouping of local layers (zones defined by the end-user) allows heterogeneous patterns to be addressed. The optimization of the different layers in a specified order and pattern is presented in [Table 1](), where 0 or 1 indicates which soil layer is optimized simultaneously at each specific step. For example, layers containing the number 1 in [Table 1]() show the grouping of different layers (zones) in which the soil hydraulic parameters have the same optimal value (i.e., homogeneous layer).

In the first step (*Opt_1*) it is assumed that the soil is homogeneous, and therefore the whole profile is modelled with five optimized soil hydraulic parameters, and the same values were given for each parameter in the different layers of the soil profile. The derived effective optimal hydraulic parameters will be used for hydrological models requiring only one homogeneous layer. In the second step (*Opt_2*) (optimal hydraulic parameters for models requiring two layers), only the parameters of the upper half of the profile are optimized, maintaining the bottom half (below the root zone) with the value of optimized parameters derived from the previous step; and then the third step (*Opt_3*) (optimal hydraulic parameters for models requiring three layers, etc) operates on the deeper zone, keeping the upper half of the profile with the previously optimized parameters.

Optimizing from top to bottom produce better results than from bottom to top, because water percolates downwards and so a change of the hydraulic parameters of the top layer will affect the lower layer. Note that the top layer is also the layer where the water content has larger variations with time. Each zone of the profile is successively split into two zones from top to bottom, and the optimization is repeated by copying the values of the optimal parameters from the previous optimization step. The vertical multistep optimization assumes that the number of soil layers, $iL$, corresponds to the number of $θ$ observations, located at the centre of each soil layer.

<figcaption align = "center"><b>Table 1a - The vertical multistep optimization method for an odd number of measurement depths is described, where the order in which the optimization is performed for the different groups of layers, iL, is presented for five measurement depths. The layers that correspond to 0 are the cells that keep their values from the previous optimization steps. The N layer splitting ‘mimics’ the number of measurement depths.</b></figcaption></figure>

 Layers | *Opt_1* | *Opt_2*| *Opt_3* | *Opt_4* | *Opt_5* | *Opt_6* | *Opt_7* | *Opt_8* | *Opt_9* 
---|---|---|---|---|---|---|---|---|---
$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$
**iL** | 1 | 2a | 2b | 3a | 3b | 4a | 4b | 5a | 5b
**1**  | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0
**2**  | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 0 | 0 
**3**  | 1 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 0
**4**  | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0
**5**  | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 1

<figcaption align = "center"><b>Table 1b - The vertical multistep optimization method for pair number of measurements depths is described, where the order in which the optimization is performed for the different groups of layers, iL is presented for four measurement depths. The layers that correspond to 0 are the cells that keep their values from the previous optimization steps. The N layer splitting ‘mimics’ the number of measurement depths.</b></figcaption></figure>

 Layers | *Opt_1* | *Opt_2*| *Opt_3* | *Opt_4* | *Opt_5* | *Opt_6* | *Opt_7*  
---|---|---|---|---|---|---|---
$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ |$\\$ 
**iL** | 1 | 2a | 2b | 3a | 3b | 4a | 4b 
**1**  | 1 | 1 | 0 | 1 | 0 | 0 | 0 
**2**  | 1 | 1 | 0 | 0 | 1 | 0 | 0  
**3**  | 1 | 0 | 1 | 0 | 0 | 1 | 0 
**4**  | 1 | 0 | 1 | 0 | 0 | 0 | 1 

## Improvement of the accuracy by increasing the number of layers

Relative errors are presented below. We assume that the HyPix model computed with the greatest number of layers, *iL*, is reference (free of error), *ref*, and models computed with a reduced number of layers are simulated, *sim*. The reduction of error in different steps of computation (a) drainage [Eq. (1)](), (b) evapotranspiration [Eq. (2)]() and (c) root zone soil water content (top 600 mm) [Eq. (3)]() is possible by increasing the number of layers, *iL*, as described in [Table 1]().

Relative error (%) for *drainage*
$$\begin{equation}
\zeta_Q = \frac{\left| \sum_{t=1}^{t=N_t}{\varDelta T^t\,\,Q_{sim_{N_z+\frac{1}{2}}}^{t}-\sum_{t=1}^{t=N_t}{\varDelta T^t Q_{ref_{N_z+\frac{1}{2}}}^{t}}} \right|}{\sum_{t=1}^{t=N}{\varDelta T^t Q_{ref_{N_z+\frac{1}{2}}}^{t}}}
\end{equation}$$


Relative error (%) for *evapotranspiration*
$$\begin{equation}
\begin{cases}	\varDelta Et = \varDelta Evap + \varDelta Rwu \\	\zeta_{et} = \frac{\left| \sum_{t=1}^{t=N_t}{\varDelta Et_{sim}^t-\sum_{t=1}^{t=N_t}{\varDelta Et_{ref}^t}} \right|}{\sum_{t=1}^{t=N_t}{\varDelta Et_{ref}^t}}\\
\end{cases}
\end{equation}$$

Relative error (%) for *root zone soil water content* 
$$\begin{equation}
\zeta _{swc}=\frac{\left| \sum_{t=1}^{N_t}{\sum_{\mathrm{i}=1}^{Z_{Nroot}}{\theta_{sim_{i}}^{t}-}}\sum_{t=1}^{N_t}{\sum_{\mathrm{i}=1}^{Z_{Nroot}}{\theta_{ref_{i}}^{t}}} \right|}{\sum_{t=1}^{N_t}{\sum_{\mathrm{i}=1}^{Z_{Nroot}}{\theta_{ref_{i}}^{t}}}}
\end{equation}$$

## Weighted objective function of the multistep optimization

The global optimizer searches for the optimal hydraulic parameters by minimizing a weighting objective function, $WOF$. We selected the robust global optimiser [BlackBoxOptim v0.5.0](https://github.com/robertfeldt/BlackBoxOptim.jl) written in Julia [(Bezanson *et al*. 2017)](#_ENDREF_47) and selected the *adaptive_de_rand_1_bin_radiuslimited* method. 

We designed the *WOF* to address the issue that observed $θ$, $θ_{obs}$ [L<sup>3</sup> L<sup>-3</sup>], in deeper layers (below the rooting zone), contains greater uncertainties compared to $θ_{obs}$ measurements in the rooting zone. It is also important that simulated $θ$, $θ_{sim}$ [L<sup>3</sup> L<sup>-3</sup>], gives more accurate predictions in the rooting zone compared to below the rooting zone. This is because all the dynamism of soil evaporation and root water uptake takes place in the rooting zone. Prioritizing the root zone by optimizing from top to bottom is performed in the algorithm used by the vertical multistep optimization (see above). This is performed by introducing a weighting, $W$, in the $WOF$ by assuming the measurements placed at different depths are equally spaced. The $WOF$ is computed as follows:

$$\begin{equation}
\begin{cases}	W_{iL}=2\frac{N_{iL}+1-iL}{N_{iL}\left( N_{iL}+1 \right)}\\	\sum_{iL=1}^{N_{iL}}{W_{iL}=1}\\	WOF_{\theta}=\sqrt{\frac{\sum_t^{N_t}{\sum_{iL}^{N_{iL}}{W_{iL}\left| \theta obs_{iL}^{t}-\theta sim_{iL}^{t} \right|^2}}}{N_{iL}\,\,N_t}}\\\end{cases}
\end{equation}$$

where $N_{iL}$ is the total number of layers, $iL$ , where $θ$ is measured; $N_t$ is the total number of time-steps; and $θ_{obs}$ and $θ_{sim}$  are observed and simulated $θ$, respectively.

