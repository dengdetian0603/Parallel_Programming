# thread 1 2 4 8 16
iter0 = c(6670, 4279, 2275, 798, 379)
thread = c(1,2,4,8,16)
t0 = c(0.61,0.59, 0.52, 0.51, 0.50)
t1 = c(26.91, 28.73, 18.91, 12.56, 9.66)
t2 = c(28.61, 18.38, 12.91, 11.98, 15.35)
t3 = c(33.61, 18.38, 13.35, 11.91, 10.98)

plot(thread, (iter0[1]+2200)/(iter0+2200),type="l",ylab="Speed up: Ts/Tl",xlab="Number of threads")
lines(thread, t1[1]/t1,col="blue")
lines(thread, t2[1]/t2,col="red")
lines(thread, t3[1]/t3,col="purple")
legend("bottomright",legend=c("pSCD simu","sCCD simu","pCCD simu","pCCD fMRI"),col=c("black","blue","red","purple"),lty=1)

t0 = c(1,1.005, 0.785, 0.85, 0.80)
t2 = c(1, 0.85, 0.82, 0.79, 0.66)
t1 = c(1, 0.88, 0.91, 0.78, 0.75)

plot(log(thread,2), t0, type="l",ylab="Scale up: Ts/Tn",xlab="log Number of threads",ylim=c(0.5,1))
lines(log(thread,2), t1, col="blue")
lines(log(thread,2), t2, col="red")
legend("bottomleft",legend=c("pSCD simu","sCCD simu","pCCD simu"),col=c("black","blue","red"),lty=1)


#datX = matrix(rnorm(10000*5000,2,4),nrow=10000)
#Y = datX%*%matrix(rep(c(1,-1,0,0,0),1000),ncol=1)
#datY = Y + rnorm(10000, 1, 3)