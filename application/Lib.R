EstimatePTVPCFactorsPenAllBIC <- function(X, Z, r, center=mean, scale=NULL, penAll=TRUE, ptvcste=FALSE)
{
	n <- ncol(X);
	T <- nrow(X);
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))));
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/");
	if (!is.null(center)) Z <- sweep(Z, 2, apply(Z, 2, function(x) do.call(center, list(x))));
	if (!is.null(scale)) Z <- sweep(Z, 2, apply(Z, 2, function(x) do.call(scale, list(x))), FUN="/");
	tmp <- eigen(t(X)%*%X);
	Lambda <- tmp$vec[,1:r];
	F <- X%*%tmp$vec[,1:r];

	F <- scale(F, center=FALSE)

	X <- Z;
	n <- ncol(X);

	res <- list();
	for (i in 1:n)
	{
		print(i);
		y <- X[,i];
		res[[i]] <- ptvfit(y~F, peninit=penAll, postest=FALSE, ptvcste=ptvcste);
	}

	ptvRawCoef <- array(0, c(n, r, T));
	ptvLambda <- array(NA, c(n, r, T));
	ptvRel <- array(0, c(n, r, T));
    intercepts <- c()
	tuningpars <- c();
	for (i in 1:n)
	{		
		tuningpars <- c(tuningpars, res[[i]]$lambda);
		if (ptvcste) {
            intercepts <- cbind(intercepts, as.numeric(diff(res[[i]]$coef[,1])!=0))
            coef <- res[[i]]$coef[,-1]
        } else {
            coef <- res[[i]]$coef
        }
		for (t in 1:T)
		{
			ptvLambda[i, , t] <- coef[t, ]
			if (t>1) ptvRawCoef[i, , t] <- coef[t,] - coef[t-1,]; 
			if (t>1) ptvRel[i, , t] <- (coef[t,] - coef[t-1,])/coef[t-1,]; 
		}
	}		
	return(list(Lambda=ptvLambda, RawCoef=ptvRawCoef, Rel=ptvRel, tuningpars=tuningpars, intercepts=rowSums(intercepts)));
}

EstimatePCFactors <- function(X, r, center=mean, scale=NULL)
{
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))));
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/");
	tmp <- eigen(t(X)%*%X);
	Lambda <- tmp$vec[,1:r];
	F <- X%*%tmp$vec[,1:r];
	return(list(Lambda=Lambda, F=F));
}

