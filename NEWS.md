# Changelog

## commits 04/05/2021 :
- version 1.0.1
- correcting package dependencies

## commits 19/03/2021 :
- default ncpu to 2
- adding ncpu to parallel functions
- changing donttest to nothing for tests that runs under 5 seconds
- adding license in repo
- adding a readme
- adding CI CD

## commits 02/03/2021 :
- change the authors format + replace dontrun examples by donttest (WARNING : resulting error in package compilation because examples use more than 2 nodes)

## commits 01/03/2021 :
- change the description of the DESCRIPTION file after first review on 01/03/2021

## commits 28/02/2021 :
- following pre-test notes on CRAN on 26/02/2021 : deletion of foreach and doSNOW package in code and dependance (description file) + deletion of VignetteBuilder 
- WARNING : add VignetteBuilder in description file after adding repertory "vignette" (see devtools_history.R file sent by FileX on 26/02/2021) 

## commits 26/02/2021 :
- updates following meeting on 25/02/2021. Fix k1 parameter in reliability measures (+k1 is null by default)
- add some parameter input tests to several functions

## commits 25/02/2021 :
- same updates as 24/02/2021
- deletion of plot arguments. Add plots in examples on help pages (+ deletion). Add plots inside code.
- deletion of vignette directory (see repertory vignette sent by email on 26/02/2021)

## commits 24/02/2021 :
- update of reference manual (add drimmR-package help page)
- update following Nicolas Vergne's new remarks in mail on 24/02/2021

## commits 23/02/2021 :
- following the remarks of the meeting on 18/02, change of the attribute length of the model : length sequence -1 to model from X_0 to X_n in dmmsum function 
- add a correction on reference manual

## commits 22/02/2021 :
- modifications on reference manual following Vlad Barbu's corrections sent the 20/02/2021 
- add \dontrun on examples code to pass CRAN tests
- add globals.R for global variables to pass CRAN test

## commits 19/02/2021 :
- modifications on reference manual (description page + add of optional argument output_file on several functions + add of references + corrections following tests)

## commits 18/02/2021 :
- modifications on reference manual (description page + add details to getStationaryLaw and getDistribution functions)
- correction of words_probas plot (former version : figure plot was not scaled on the selected frame)

## commits 17/02/2021 :
- add of reference manual
- add of vignette (still latex errors because of bibtex file)
- deletion of word_expect function (not started yet)

## commits 16/02/2021 :
- correction of words_probas roxygen imports (because filter function was not recognized)
- add a return of simulated sequence for simulate function
- update of each help of functions

## commits 15/02/2021 :
- from package version 0.2.0 to 1.0.0 in description file
- add help descriptions and details for dmmsum, errorRATE, availability, maintainability and reliability functions
- correction of custom initial law conditions for customisable vector of initial law (init.estim=vector of length nb.state^order)
- correction of word applications and failure rates (toy object "mod" inside code instead of "x")
- word_expect application still in progress : do not use yet
- add of allgenerics.R file gathering all S3 generic functions. Deletion of former S3 generic function files (ex: simulate.R)


## commits 12/02/2021 :
- deletion of the plot.R file
- add of the words_probas function


## commits 10/02/2021 and 11/02/2021 :
- delete 2 plot functions and insert them into word applications
- MAJ of 3 word applications : word_proba, word_probas and length_probas
