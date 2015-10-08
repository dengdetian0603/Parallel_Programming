0. change to the source directory

1. Install R by:
$	sudo apt-get update
$	sudo apt-get install r-base
$	sudo apt-get install r-base-dev

2. Install required R packages by:
$	Rscript installpackages.R

3. Compile source code by:
$	sudo R CMD INSTALL buckshot

4. Run desired experiment by the following command, replacing path with the directory to data, nthreads with the number of threads to use, maxiter with the maximum iteration allowed, and solver with an integer specifying the optimization method to use: 0-pSCD 1-sCCD 2-pCCD 
$	Rscript run_buckshot.R path nthreads maxiter solver

Due to data confidentiality, fMRI data used for project cannot be released to unauthorized person.