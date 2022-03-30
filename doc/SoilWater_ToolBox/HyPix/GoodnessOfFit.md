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


# Evaluating goodness of fit

The goodness of fit of the measured soil water content ($θ_{obs}$) to the simulated ($θ_{sim}$) at the corresponding depths using the HyPix model was assessed using the root mean squared error, $RMSE$, Nash–Sutcliffe efficiency coefficient, $NSE$, and the refined index of agreement proposed by [Willmotts *et al*. (2012)](#_ENDREF_48), $dr$, as follows:   

$$\begin{equation}
RMSE=\sqrt{\frac{\sum_{i=1}^{N_t}{\left[ \theta_{obs_i}-\theta_{sim_i} \right] ^2}}{N_t}}
\end{equation}$$

$$\begin{equation}
NSE=\max \left\{ 1-\frac{\sum_{i=1}^{N_t}{\left[ \theta_{obs_i}-\theta_{sim_i} \right] ^2}}{\sum_{i=1}^{N_t}{\left[ \theta_{obs_i}-\overline{\theta_{obs}} \right] ^2}};0 \right\}
\end{equation}$$

$$\begin{equation}
dr=\left\{ \begin{array}{c}	1-\frac{\sum_{i=1}^{N_t}{\left| \theta _{sim_i}-\theta _{obs_i} \right|}}{c\sum_{i=1}^{N_t}{\left| \theta _{obs_i}-\overline{\theta _{obs}} \right|}}, \mathrm{when} \sum_{i=1}^{N_t}{\left| \theta _{sim_i}-\theta _{obs_i} \right|}\leqslant c\sum_{i=1}^{N_t}{\left| \theta _{obs_i}-\overline{\theta _{obs}} \right|}\\	\frac{c\sum_{i=1}^{N_t}{\left| \theta _{obs_i}-\overline{\theta _{obs}} \right|}}{\sum_{i=1}^{N_t}{\left| \theta _{sim_i}-\theta _{obs_i} \right|}}-1, \mathrm{when} \sum_{i=1}^{N_t}{\left| \theta _{sim_i}-\theta _{obs_i} \right|}>c\sum_{i=1}^{N_t}{\left| \theta _{obs_i}-\overline{\theta _{obs}} \right|}\\\end{array} \right. 
\end{equation}$$

with $c$ = 2 in [Eq. (3)](). The $dr$ [-] index, with values between −1.0 and 1.0, indicates the sum of the magnitudes of the differences between the model-predicted and observed deviations about the observed mean relative to the sum of the magnitudes of the perfect-model ($θ_{obs}$ = $θ_{sim}$ for all $i$) and observed deviations about the observed mean. This index is used because, in general, more rationally related to model accuracy than are other existing indices [(Willmotts *et al*., 2012)](#_ENDREF_48).