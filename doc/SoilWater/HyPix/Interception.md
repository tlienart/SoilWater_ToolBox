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


# INTERCEPTION MODEL

The parsimonious physically based interception model is an improvement of [Pollacco *et al.* (2013a)](#_ENREF_21). The following interception model uses *potential evaporation of a wet canopy* ([$ \varDelta Pet_{int} $ [L]](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYPIX/Potential_evapotranspiration)), $LAI$ [-], and *extinction coefficient for solar radiation* ($K_{g}$ [-]) set to 0.5. The *gross precipitation depth* that falls on top of a canopy, [L], is partitioned following [Rutter *et al.* (1971)](#_ENREF_22) as:

$$\begin{equation} 
	\varDelta Pr\,\,=\,\,\varDelta Pr_{int}+\varDelta Pr_{ground}
\end{equation}$$

where  $ \varDelta Pr_{ground} $ [L] is the fraction of precipitation reaching the soil surface through gaps in the canopy, and $\varDelta Pr_{int}$ [L] is the *intercepted precipitation depth.* They are computed as:

$$\begin{equation}
 \varDelta Pr_{ground}\,\,=G_{apFrac}\,\,\varDelta Pr\\                                                               \varDelta Pr_{int}=\,\,\left[ 1-G_{apFrac} \right] \,\,\varDelta Pr 
 \end{equation}$$

where the gap fraction, $G_{apFrac}$ [-], is calculated using the Beer–Lambert law as:

$$\begin{equation}
G_{apFrac}=1-\,\,e^{-K_g⋅LAI}
 \end{equation}$$

The foliage of the canopy is considered as a water storage filled up to a depth $Sint$ [L], *with a saturated storage capacity*, $Sint_{sat}$ [L]. When the canopy is fully saturated ($Sint=Sint_{sat}$), then any excess of $\varDelta Pr_{int}$ overflows, $\varDelta Pr_{over}$ [L], to the soil surface. The amount of water that reaches the soil surface is the *throughfall precipitation* [L]:

$$\begin{equation}
 \varDelta Pr_{through}=\varDelta Pr_{ground}+\varDelta Pr_{over} 
 \end{equation}$$

The water storage of the canopy is first computed as:

$$\begin{equation}
 Sint^t=Sint^{t-1}+\varDelta Pr_{int} 
 \end{equation}$$

A fraction of the water from $Sint$ will be evaporated at the rate of the *actual evaporation* depth, $\varDelta Evap_{int}$ [L], during and after a rainfall event. The maximum quantity of water that can be evaporated from a wet canopy during a time-step is computed according to [Deardorff (1978)](#_ENREF_3), which assumes that $\varDelta Evap_{int}$ is proportional to the fraction of the canopy that is wet:

$$\begin{equation}
 \varDelta Evap_{int}=\,\,\varDelta Pet_{int}\left[ \frac{Min\left\{ Sint^t;Sint_{sat} \right\}}{Sint_{sat}} \right] ^{Pevap_{int}} 
\end{equation}$$

where $Pevap_{int}$ [-] is a constant parameter for which [Deardorff (1978)](#_ENREF_3) gives a constant value of 2/3.

$ \varDelta Pr_{over} $ is computed as:

$$\begin{equation}
\begin{cases}                                                               \varDelta Pr_{over}=\,\,Min\left[ Sint_{}^{t-1}+\varDelta Pr_{int}-\varDelta Evap_{int}-S_{sat}^{}\,\,;\,\,0 \right]\\                                                       Sint_{}^{t}=Min\left[ Sint_{}^{t-1}+\varDelta Pr_{int}-\varDelta Evap_{int}-\varDelta Pr_{over};0 \right]\\\end{cases} 
\end{equation}$$  

Rainfall interception of gross rainfall loss, $InterceptionLoss$ [0–1], is computed by:

$$\begin{equation}
 InterceptionLoss\,\,=1-\frac{\sum_{\mathrm{t}=1}^{\mathrm{N}_{\mathrm{t}}}{\varDelta Pr_{through}^{\mathrm{t}}}}{\sum_{\mathrm{t}=1}^{\mathrm{N}_{\mathrm{t}}}{\varDelta Pr^{\mathrm{t}}}} 
\end{equation}$$

In the HyPix model the rainfall interception module is run first, followed by the computation of [$\varDelta Pet_{transp}$ and $\varDelta Pet_{evap}$](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYPIX/Interceptionhttps://manaakiwhenua.github.io/SoilWater_ToolBox.jl/HYPIX/Potential_evapotranspiration).
