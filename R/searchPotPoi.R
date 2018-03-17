searchPotPoi <-
function(k, X, y, dom = c(0,1), plotting=FALSE){
    ## ########################
    ## A Function identifying potential points of impact of Cov(Z_delta, Y)
    ## Input:
    ##  - k: grid indices on domain
    ##  - X: centered X-data; Nxp MMatrix
    ##  - y: centered y-data
    ##  - dom: the domain where X(t) s live on
    ##  - plotting: binary variable if covariance should be plotted
    ## Output:
    ##  - ind.x: index of potential points of impact
    ##X <- X_st; y = Y ; k = 1
    N           <- ncol(X)
    p 		<- nrow(X)              	# number of discretization points 
    a 		<- range(dom)[1]			# left border
    b 		<- range(dom)[2]			# right border
    t 		<- seq(a,b,length.out=p) 	# grid in the domain [a,b]
    
    delta  	<- k*(b-a)/(p-1)    		# translate k to delta
    Jdelta 	<- (k+1):(p-k)      		# only 'inner' grid-point in [a,b]
    ## print(X)
    ## print(y)
    CorXY      <- 1/N*((X%*%y)) #up to factor 1/n
    CorZY      <- apply(CorXY,2,function(X){X[Jdelta]-1/2*( X[Jdelta-k] + X[Jdelta+k])})
    
    COR        <- abs(CorZY) #Criteria
    COR.aux    <- COR
    Jdelta.aux <- Jdelta
    t.delta    <- t[Jdelta]
    t.aux      <- t.delta
    taucand    <- vector() #container for hat(tau)
    ## Estimate tau!
    ##print(COR)
    while(sum(COR.aux>0.000000001)>1){
        tau.aux <- which.max(COR.aux)  
        taucand <- c(taucand,t.aux[tau.aux])
        indi    <- (abs(t.aux-t.aux[tau.aux])>=sqrt(delta)/2)# theoretical threshold
        t.aux   <- t.aux[indi]
        COR.aux <- COR[t.delta%in%t.aux]
    }
    if(plotting){
        plot(t.delta,COR,xlab="",ylab="",yaxt='n',type="l")#main=paste("k =",k, "kappa.est =",round(kappa,2), "   c =", c, "delta= ", round(delta,3)))
        ##abline(v=t[tau],col="black",lwd=1)
        abline(v=taucand,col="red")
    }
    ## Get the indexsets corresponding to the tau-candidates and return them:
    ind.x <- which(t%in%taucand)[rank(taucand)]

    list(ind.x = ind.x, cor = COR)
}
