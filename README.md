# Drift Markov Model

[![PLMLab build status](https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR/badges/master/pipeline.svg)](https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR/pipelines) [![CRAN versions](https://www.r-pkg.org/badges/version/drimmR)](https://CRAN.R-project.org/package=drimmR) [![CRAN logs](https://cranlogs.r-pkg.org/badges/drimmR)](https://CRAN.R-project.org/package=drimmR)

## Install

Latest CRAN version :

```R
install.packages('drimmR')
```

Or for the latest development version:

```R
devtools::install_git("https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR.git")
```

## Example

```R
library("drimmR")
data(lambda, package = "drimmR")
dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
simulate(dmm, model_size=100)
```

## Parallelisation

You can use multiple processors for most functions. Add `ncpu=4` to run on 4 cores. If you want R to detect the number of available cores on your system, use `ncpu=-1`.

## Get sources

You need to have git. Then:
```
git clone https://plmlab.math.cnrs.fr/lmrs/statistique/drimmR.git
```

## Build

```bash
R CMD build .
R CMD check --as-cran .
```

## License

GPL v3
