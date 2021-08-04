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

 The representation of the $θ(ψ)$ and $K(θ)$ functions is based on the dual porosity models of [Pollacco *et al*. (2017)](#_ENDREF_40), and we use the bimodal lognormal expressions of the Kosugi model from [Fernández-Gálvez *et al.*, (2021)](#_ENDREF_41): 

$$\begin{equation}
\begin{cases}	\theta \left( \psi \right) =\theta _{\mathrm{Mat}}\left( \psi \right) +\theta _{\mathrm{Mac}}\left( \psi \right)\\	\theta _{\mathrm{Mat}}\left( \psi \right) =\frac{1}{2}\left( \theta _{s\mathrm{MacMat}}-\theta _r \right) \,\,erfc\left( \frac{\ln \psi  -\ln  \psi _m}{\sigma \sqrt{2} } \right) +\theta _r\\	\theta _{\mathrm{Mac}}\left( \psi \right) =\frac{1}{2}\left( \theta _s-\theta _{s\mathrm{MacMat}} \right) \,\,erfc\left( \frac{\ln \psi  -\ln  \psi _{m\mathrm{Mac}}}{\sigma _{\mathrm{Mac}}\sqrt{2} } \right)\\\end{cases}
\end{equation}$$

$$\begin{equation}
\begin{cases}	K\left( \psi \right) =K_{\mathrm{Mat}}\left( S_{e\,\,\mathrm{Mat}}\left( \psi \right) \right) +K_{\mathrm{Mac}}\left( S_{e\,\,\mathrm{Mac}}\left( \psi \right) \right)\\	S_e=\frac{\theta -\theta _r}{\theta _s-\theta _r}\\	K_{\mathrm{Mat}}\left( \psi \right) =K_s\frac{\theta _{s\mathrm{MacMat}}-\theta _{\mathrm{r}}}{\theta _s-\theta _r}\,\,\sqrt{S_{e\,\,\mathrm{Mat}}\left( \psi \right)}\,\,\left[ \frac{1}{2}erfc\left( erfc^{-1}\left( 2S_{e\,\,\mathrm{Mat}}\left( \psi \right) \right) +\frac{\sigma}{\sqrt{2}} \right) \right] ^2\\	K_{\mathrm{Mac}}\left( \psi \right) =K_s\frac{\theta _{\mathrm{s}}-\theta _{s\mathrm{MacMat}}}{\theta _s-\theta _r}\,\,\sqrt{S_{e\,\,\mathrm{Mac}}\left( \psi \right)}\,\,\left[ \frac{1}{2}erfc\left( erfc^{-1}\left( 2S_{e\,\,\mathrm{Mac}}\left( \psi \right) \right) +\frac{\sigma _{\mathrm{Mac}}}{\sqrt{2}} \right) \right] ^2\\\end{cases}
\end{equation}$$

where $erfc$ is the complementary error function; $θ$ [L<sup>3</sup> L<sup>-3</sup>] represents the volumetric soil water content and $ψ$ [L] the soil water pressure, considering $ψ > 0$ for unsaturated soils (i.e. matric suction); $θ_s$ [L<sup>3</sup> L<sup>-3</sup>] and $θ_r$ [L<sup>3</sup> L<sup>-3</sup>] are the saturated and residual volumetric soil water content, respectively; ln $ψ_m$ [L] (with the argument of ln in units of length, i.e., $ψ_m$ in [L]) and $σ$ [-] denote the mean and standard deviation of ln $ψ$ [L], respectively in the soil matrix domain; ln $ψ_{mMac}$ [L] and $σ_{Mac}$ [-] denote the mean and standard deviation of ln $ψ$, respectively, in the macropore soil domain; $θ_{sMacMat}$ [L<sup>3</sup> L<sup>-3</sup>] is the volumetric saturated water content that theoretically differentiates inter-aggregate pores (structural macropores) and matrix domains (intra-aggregate micropores), defining the corresponding soil water pressure threshold between macropore and matrix $ψ_{MacMat}$ [L]; $K_s$ [L T<sup>-1</sup>] is the saturated hydraulic conductivity; and $S_{e\ Mat}$ [-] and $S_{e\ Mac}$ [-] denote the effective saturation as a function of $ψ$ in the soil matric and macropore domains, with values between 0 and 1. Finally, $K(ψ)$ refers to the unsaturated hydraulic conductivity, written as a function of $ψ$.
