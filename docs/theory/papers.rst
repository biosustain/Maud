================
Relevant Reading
================

These are some things you might like to read if you are interested in Maud.


Theory of Enzyme kinetics
=========================

- Monod, J., Wyman, J., & Changeux, J. (1965). On the nature of allosteric
  transitions: A plausible model. Journal of Molecular Biology, 12(1),
  88–118. http://dx.doi.org/10.1016/S0022-2836(65)80285-6

  The original version of Maud's allostery framework.

- Popova, S. V., & Sel'kov, E. E. (1975). Generalization of the model by monod,
  wyman and changeux for the case of a reversible monosubstrate reaction. FEBS
  Letters, 53(3), 269–273. http://dx.doi.org/10.1016/0014-5793(75)80034-2

  This paper sets out the generalised MWC framework that Maud uses to represent
  allosteric regulation.

- Saa, P. A., & Nielsen, L. K. (2017). Formulation, construction and analysis
  of kinetic models of metabolism: A review of modelling
  frameworks. Biotechnology Advances, 35(8),
  981–1003. http://dx.doi.org/10.1016/j.biotechadv.2017.09.005

  A review paper covering different approaches to modelling metabolic networks.


Themodynamics
=============

- Mahamkali, V., McCubbin, T., Beber, M., Marcellin, E., & Nielsen,
  L. K. (2020). multiTFA: a Python package for multi-variate
  Thermodynamics-based Flux Analysis. bioRxiv, (),
  2020–12–01–407387. http://dx.doi.org/10.1101/2020.12.01.407387

- Alberty, R. A. (2003). Thermodynamics of biochemical reactions. Hoboken, N.J:
  Wiley-Interscience.

  Chapter 4 explains how to adjust thermodynamic measurements to take into
  accound experimental conditions.

- Noor, E., Haraldsd\'ottir, Hulda S., Milo, R., & Fleming,
  R. M. T. (2013). Consistent Estimation of Gibbs Energy Using Component
  Contributions. PLoS Computational Biology, 9(7),
  1003098. http://dx.doi.org/10.1371/journal.pcbi.1003098

  Original paper for the method behind equilibrator.

- Du, B., Zielinski, D. C., & Palsson, B. O. (2018). Estimating Metabolic
  Equilibrium Constants: Progress and Future Challenges. Trends in Biochemical
  Sciences, 43(12), 960–969. http://dx.doi.org/10.1016/j.tibs.2018.09.009

  Review paper discussing the problem of estimating thermodynamic parameters in
  general.

- Du, B., Zhang, Z., Grubner, S., Yurkovich, J. T., Palsson, B. O., &
  Zielinski, D. C. (2018). Temperature-Dependent Estimation of Gibbs Energies
  Using an Updated Group-Contribution Method. Biophysical journal, 114(11),
  2691–2702. http://dx.doi.org/10.1016/j.bpj.2018.04.030

  Develops the component contribution method in several ways. There is a nice
  dataset in an excel sheet in the supplementary information.


Computational Models of metabolic networks
==========================================

- Shepelin, D., Machado, D., Nielsen, L. K., & Herrg\aard, Markus
  J. (2020). Benchmarking kinetic models of Escherichia coli
  metabolism. bioRxiv, (),
  2020–01–16–908921. http://dx.doi.org/10.1101/2020.01.16.908921

  A review paper that compares various computational models of E. coli

- St. John, P., Strutz, J., Broadbelt, L. J., Tyo, K. E. J., & Bomble,
  Y. J. (2018). Bayesian inference of metabolic kinetics from genome-scale
  multiomics data. bioRxiv, (), . http://dx.doi.org/10.1101/450163

  This paper takes a very similar approach to Maud, but represents reactions
  using lin-log kinetics.

- Kacser, H., & Burns, J. A. (1973). The control of flux.
  Symposia of the Society for Experimental Biology, 27, 65–104.

  Implementation of Metabolic Control Analysis

- Steuer, R., & Junker, B. H. (2008). Computational Models of Metabolism: 
  Stability and Regulation in Metabolic Networks. In S. A. Rice (Ed.), 
  Advances in Chemical Physics (pp. 105–251). John Wiley & Sons, Inc. 
  https://doi.org/10.1002/9780470475935.ch3

  Great description of implementation of MCA as done in Maud.



Statistics
==========

- Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B., Yao,
  Y., Kennedy, L., … (2020). Bayesian Workflow. arXiv:2011.01808 [stat], (), .

  A general discussion of the workflow that Maud is intended to be part of.


Algorithms
==========

- Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt,
  M., Brubaker, M., … (2017). Stan: A Probabilistic Programming
  Language. Journal of Statistical Software, 76(1),
  1–32. http://dx.doi.org/10.18637/jss.v076.i01

  Introduces Stan and the ideas behind it

- Margossian, C. C. (2019). A Review of automatic differentiation and its
  efficient implementation. WIREs Data Mining and Knowledge Discovery, 9(4),
  . http://dx.doi.org/10.1002/WIDM.1305

  Relatively recent and accessible discussion of automatic differentiation.

- Carpenter, B., Hoffman, M. D., Brubaker, M., Lee, D., Li, P., & Betancourt,
  M. (2015). The Stan Math Library: Reverse-Mode Automatic Differentiation in
  C++. arXiv:1509.07164 [cs], (), .

  Explains Stan's automatic differentiation implementation.

- Lyness, J. N., & Moler, C. B. (1967). Numerical Differentiation of Analytic 
  Functions. SIAM Journal on Numerical Analysis, 4(2), 202–210.
  https://doi.org/10.1137/0704019

  Explains the complex number differentiation

