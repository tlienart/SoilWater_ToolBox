<!-- MathJax -->

<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

# International SoilWater-ToolBox 2021
## *Current state and future*


J.A.P. Pollacco <sup> 1 </sup>, J. Fernández-Gálvez  <sup> 2 </sup>, L. Lilburne <sup> 1  </sup>, S. Carrick  <sup> 1  </sup>, S. McNeill  <sup> 1  </sup>, D.A. Peltzer <sup> 1  </sup> B. Belfort <sup> 3  </sup>, P. Ackerer<sup> 3  </sup>, L. Lassabatere <sup> 4  </sup>, R. Angulo-Jaramillo <sup> 4  </sup>, S.C. Zammit <sup> 5  </sup>, C. Rajanayaka <sup> 5  </sup>
---

> <sup> 1 </sup> Manaaki Whenua -- Landcare Research, Lincoln 7608, **New Zealand**

> <sup> 2 </sup> Department of Regional Geographic Analysis and Physical Geography, University of Granada, Granada, **Spain**

> <sup> 3 </sup> Université de Strasbourg, CNRS/EOST, ITES UMR 7063, Institut Terre et Environnement de Strasbourg, Strasbourg, **France**

> <sup> 4 </sup> Univ Lyon, Université Claude Bernard Lyon 1, CNRS, ENTPE, UMR5023 LEHNA, Vaulx en Velin, Lyon 69518, **France**

> <sup> 5 </sup> National Institute of Water and Atmospheric Research, Christchurch, **New Zealand**


---

---

---

## 1. Open source software

The open source **SoilWater-ToolBox** software is written in the performant and readable Julia language ([https://julialang.org/](https://julialang.org/)). It can be downloaded from [https://github.com/manaakiwhenua/SoilWater_ToolBox/](https://github.com/manaakiwhenua/SoilWater_ToolBox/) and is available under the **GP-3.0 License**. This software includes a set of interlinked modules that can also be used independently.

## 2. Mission statement

The aim of the **SoilWater-ToolBox** is to derive soil hydraulic parameters and soil water fluxes using wide range of physically based, cost-effective methods. The estimated hydraulic parameters can be directly implemented into the physically based Hydrological Pixel (**HyPix**) model to compute the soil-water balance. The **HyPix** model can also be used to derive the soil hydraulic parameters from time series of *soil water content* measurements. The **SoilWater-ToolBox** enables the user to perform inter-comparison and sensitivity of the hydraulic parameters computed from different methods on soil-water fluxes of interest.

## 3.  SoilWater-ToolBox based on peer-reviewed publications

The following modules are implemented in the **SoilWater-ToolBox.** They were developed for the following specific scientific purposes, based on peer-reviewed scientific publications:

- **the Intergranular Mixing Particle size distribution (IMP) model,** which derives unimodal hydraulic parameters using particle size distribution (Pollacco et al., 2020)
- **the General Beerkan Estimation of Soil Transfer parameters method,** which derives the unimodal hydraulic parameters from single ring infiltration experiments (Fernández-Gálvez et al., 2019)
- **the HyPix model,** a physically based hydrological model which solves the Richards equation, and computes *soil-water content*, *transpiration*, *evaporation*, and *drainage* when input with climate data and hydraulic parameters (Pollacco et al., 2021)
- **the sorptivity model,** a novel procedure for straightforward computation of sorptivity (Lassabatere et al., 2021), which is implemented into the General Beerkan Estimation of Soil Transfer parameters method (Lassabatere et al., 2021) and in the **HyPix** model (Pollacco et al., 2021)
- **to derive saturated hydraulic conductivity** from unimodal and bimodal *θ*(*ψ*) (Pollacco et al., 2017, 2013)
- **to invert hydraulic parameters** from *soil moisture* time series using the **HyPix** model (Pollacco et al., 2021)
- **to derive unique and physical bimodal Kosugi hydraulic parameters** from inverse modelling (Fernández-Gálvez et al., 2021) using water retention and/or unsaturated hydraulic conductivity data directly measured in the laboratory, or indirectly obtained from inverting *soil moisture* time series using the **HyPix** model (Pollacco et al., 2021).

## 4.  Schematic flow chart

A simplified schematic flow chart of the **SoilWater-ToolBox** is provided below. It shows the complex relationship between the different modules, the required input data and outputs prodcued by the **SoilWater-ToolBox.** .

![HyPix](https://manaakiwhenua.github.io/SoilWater_ToolBox.jl/FIGURE/SoilWater-ToolBox-FlowChart.bmp "SoilWater-ToolBox Flowchart")

## 5.  Current applications

The management of both soil and water resources is of primarily importance. Current pressure on natural resources together with climate trends, increase the need for better understanding and predicting the movement of water in the soil. Characterizing soil properties to model soil water dynamic is also directly link to the water cycle with direct implications into the environment, management of crop production, and related socio-economical aspects.

The **SoilWater-ToolBox** is being used in several projects in New Zealand (Manaaki Whenua -- Landcare Research, NIWA, Plant & Food), as well as in France and Spain, for:

- **laboratory data**, to derive unique sets of physical hydraulic parameters from laboratory data even when key data are missing, such as the *unsaturated hydraulic conductivity*;
- **S-Map-Hydro** **across New-Zealand**, to derive physical hydraulic parameters suitable for a wide range of hydrological models -- the derived hydraulic parameters are **scaled** to the vertical scale of interest and corrected for **stone** content;
- **validating/adjusting S-Map-Hydro** by feeding the hydraulic parameters into a physical hydrological model, **HyPix,** and comparing the *soil-water content* outputs with measured values;
- **inverting hydraulic parameters**, to derive hydraulic parameters by inverting time series *soil-water content* data (e.g. TDR, soil moisture capacitance sensor FDR, neutron probe);
- the **automatic infiltrometer**, to derive physical hydraulic parameters from automatic infiltration tests ;
- **particle size distribution**, to derive physical hydraulic parameters exclusively from *soil particle size distribution*.

## 6.  Potential future applications

Potential future applications include:

- **SoilWater-HyPix-2D:** spatialising the **HyPix** model such that it predicts spatially distributed water balance;
- **wilding pines:** quantifying the impact of the succession of different vegetation types on the hydrological balance;
- **particle size distribution from laser:** deriving bimodal hydraulic parameters that account for the matrix and macropore domains of soils;
- **National Soils Data Repository (NSDR):** automatically feeding soil data from the NSDR into the **SoilWater-ToolBox;**
- **Land Cover Database (LCDB):** deriving vegetation parameters automatically from remote sensing and the **LCDB;**
- **infiltration data from NSDR:** automatically feeding data into the software for soil hydraulic characterization;
- **additional vegetation data from remote sensing:** when available, this will be implemented to improve the transpiration module and its impact on the soil-water balance;
- **SoilWater-HyPix-CenW:** coupling **HyPix** with **CenW** model -- the **HyPix** hydrological model has an advanced unsaturated module which accurately predicts the movement of water in the unsaturated soil; **CenW** has an advanced comprehensive forest growth model based on linked flows of *carbon*, *energy*, *nutrients* and *water* in trees and the soil.

## 7.  Advantages of the Julia language

**The Julia language:**

- runs as fast as C+, and is suitable for hydrological modelling,
- is suitable for use with large data sets and parallel computing,
- supports encoding via Unicode, UTF-8,
- harmonises different packages, making code easy to read,
- allows interoperability with other programming languages, such as C, Fortran, R and Python,
- facilitates package management,
- has a great Julia community for support,
- has mature libraries.

## 9. Further links (tests)

For information on the author [Authors](https://manaakiwhenua.github.io/SoilWater_ToolBox/Authors)

For information on HyPix model [HyPixhttps://github.com/manaakiwhenua/SoilWater_ToolBox](https://manaakiwhenua.github.io/SoilWater_ToolBox/HYPIX/HyPix_Introducing)

## TESTING EQUATIONS

$$
\frac{A}{B}

$$

## 8.  Publications

Fernández-Gálvez, J., Pollacco, J.A.P., Lassabatere, L., Angulo-Jaramillo, R., Carrick, S., 2019. A general Beerkan Estimation of Soil Transfer parameters method predicting hydraulic parameters of any unimodal water retention and hydraulic conductivity curves: Application to the Kosugi soil hydraulic model without using particle size distribution data. Advances in Water Resources 129, 118--130. https://doi.org/10.1016/j.advwatres.2019.05.005Fernández-Gálvez,

J., Pollacco, J.A.P., Lilburne, L., McNeill, S., Carrick, S., Lassabatere, L., Angulo-Jaramillo, R., 2021. Deriving physical and unique bimodal soil Kosugi hydraulic parameters from inverse modelling. Advances in Water Resources 153, 103933. https://doi.org/10/gkbdsxLassabatere,

L., Peyneau, P.-E., Yilmaz, D., Pollacco, J., Fernández-Gálvez, J., Latorre, B., Moret-Fernández, D., Di Prima, S., Rahmati, M., Stewart, R.D., Abou Najm, M., Hammecker, C., Angulo-Jaramillo, R., 2021. Scaling procedure for straightforward computation of sorptivity. Hydrology and Earth System Sciences. https://doi.org/10.5194/hess-2021-150Pollacco,

J.A.P., Fernández-Gálvez, J., Carrick, S., 2020. Improved prediction of water retention curves for fine texture soils using an intergranular mixing particle size distribution model. Journal of Hydrology 584, 124597. https://doi.org/10.1016/j.jhydrol.2020.124597Pollacco,

J.A.P., Fernandez-Galvez, J., Carrick, S., McNeill, S., Peltzer, D.A., Lassabatere, Laurent, Raphael, A.-J., Ackerer, P., Belfort, B., Zammit, C., Channa, R., 2021. HyPix: 1D Richards equation hydrological model in Julia language using a multistep optimization scaling method for flexible soil vertical discretization. Submitted to Environmental Modelling & Software.

Pollacco, J.A.P., Nasta, P., Ugalde, J.M.S., Angulo-Jaramillo, R., Lassabatere, L., Mohanty, B.P., Romano, N., 2013. Reduction of feasible parameter space of the inverted soil hydraulic parameters sets for Kosugi model. Soil Science SS-S-12-00268.

Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725--2737. https://doi.org/10.5194/hess-21-2725-2017
