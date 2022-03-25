# Drift Markov Model

[![PLMLab build status](https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR/badges/master/pipeline.svg)](https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR/-/pipelines) [![CRAN versions](https://www.r-pkg.org/badges/version/drimmR)](https://CRAN.R-project.org/package=drimmR) [![CRAN logs](https://cranlogs.r-pkg.org/badges/drimmR)](https://CRAN.R-project.org/package=drimmR) [![RDocumentation](https://api.rdocumentation.org/badges/version/drimmR)](https://www.rdocumentation.org/packages/drimmR)

## Description

Performs the drifting Markov models (DMM) which are 
non-homogeneous Markov models designed for modeling the heterogeneities of
sequences in a more flexible way than homogeneous Markov chains or even 
hidden Markov models. In this context, we developed an R package dedicated to 
the estimation, simulation and the exact computation of associated reliability 
of drifting Markov models. The implemented methods are described in:

* Vergne, N. (2008). Drifting Markov models with polynomial drift and applications to DNA sequences. Statistical applications in genetics and molecular biology, 7(1). doi:10.2202/1544-6115.1326. [Journal version](https://www.degruyter.com/document/doi/10.2202/1544-6115.1326/html)

* Barbu, V. S., & Vergne, N. (2019). Reliability and survival analysis for drifting Markov models: modeling and estimation. Methodology and Computing in Applied Probability, 21(4), 1407-1429. doi:10.1007/s11009-018-9682-8. [Journal version](https://link.springer.com/content/pdf/10.1007/s11009-018-9682-8.pdf)

## Contribute

The official repository is at [PLMLab](https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR/). But to help with issues and contributions, a mirror has been setup at [Github](https://github.com/corentin-dev/drimmR).

## Install

* Install from CRAN:

```R
install.packages('drimmR')
```

* Install latest development version from `git`:

```R
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_git("https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR", dependencies = TRUE, build_vignettes = FALSE)
```

## Getting started

```R
library("drimmR")
data(lambda, package = "drimmR")
dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
simulate(dmm, model_size=100)
```

## Parallelisation

You can use multiple processors for most functions. Add `ncpu=4` to run on 4 cores. If you want R to detect the number of available cores on your system, use `ncpu=-1`.

## License

GPL v3
