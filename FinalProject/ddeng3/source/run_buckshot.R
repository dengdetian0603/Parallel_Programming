args = commandArgs(TRUE)

datadir = as.character(args[1])
nthreads = as.numeric(args[2])
maxIter = as.numeric(args[3])
solver = as.numeric(args[4])

library(buckshot)
load(datadir)

datY = dat[,1]
datX = dat[,-1]


fitmodel = function(Solver)
{
	buckshot(x=datX, y=datY, type='lasso', lambda=1.6, path.length=Solver, max.iter=maxIter,
         convergence.threshold=1e-5, threads=nthreads, verbose=FALSE)
}

ptm = proc.time()
plassofit = fitmodel(solver)
proc.time()-ptm
