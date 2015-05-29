library(parsimonious)
library(macrods)
library(zoo)
library(xtable)
library(ggplot2)
library(reshape2)
source("Lib.R")

# Number of factors, do not change, is hardcoded in some places!
noF <- 4;
# Should the intercept be time-varying?
ptvcste <- TRUE
 
# Load data
FData <- scale(window(getmacrodata(SW2009, c(2006,4)), start=c(1959,3)))[,SW2009$include==1]
AllData <- scale(window(getmacrodata(SW2009, c(2006,4)), start=c(1959,3)))[,SW2009$include!=0]
ForecastData <- window(getmacrodata(SW2009, c(2006,4), h=4), start=c(1959,3))[,SW2009$include!=0]
PCEst <- EstimatePCFactors(FData, noF)
if (file.exists("PTVEst.RData")) {
   load("PTVEst.RData")
} else {
    PTVEst <- EstimatePTVPCFactorsPenAllBIC(FData, AllData, noF, penAll=FALSE, ptvcste=ptvcste); 
    save(PTVEst,file='PTVEst.RData')
}
F <- ts(PCEst$F, start=c(1959,3), freq=4)
aggidx <- (SW2009$include==2)[SW2009$include!=0]
aggidx[c(7,8,18,19,20,24,32,35,75,79,84,86,95)] <- FALSE
aggidx[c(12,13,29,30,44,51,64,65,66,67,99,100,104,105,106,110,111,112,113,126,127,128,129,130,132,133,134,139,143,144)] <- TRUE

vnames <- cbind(SW2009$shortnames, colnames(SW2009$rawdata), SW2009$tcodes, SW2009$include, SW2009$longnames)
vnames <- vnames[SW2009$include != 0,]
vnames[vnames[,4]=="2",4] <- 0 

if (!file.exists("figures")) dir.create("figures")
if (!file.exists("tables")) dir.create("tables")

# Factor eq:
res1 <- res2 <- res3 <- c()
for (i in 1:144) res1 <- c(res1, sum(PTVEst$RawCoef[i,,]!=0))
for (i in 1:144) res2 <- c(res2, sum(PTVEst$RawCoef[i,,43:102]!=0))
for (i in 1:144) res3 <- c(res3, sum(PTVEst$RawCoef[i,,103:162]!=0))
AllResFac <- cbind(res1, res2, res3)
rownames(AllResFac) <- colnames(AllData)

# Forecasting eq:
res1 <- res2 <- res3 <- allnz <- resstd <- fcerr <- fcols <-  c();
allmse <- c()
allres <-  allolsres <- c()
forecastintercepts <- c()

for (i in 1:ncol(AllData))
{
	print(i);
	mne <- vnames[i,2]
	lna <- vnames[i,5]
	a <- ForecastData[,i]
	b <- AllData[,i]
	z <- cbind(a, scale(cbind(F, b, lag(b, -1), lag(b, -2), lag(b, -3)), center=FALSE))
	Z <- data.frame(window(z, start=c(1960,2), end=c(2005,4)))
	colnames(Z) <- c("Y", paste("F", 1:4, sep=""), paste("Y_", 0:3, sep=""))
	res <- ptvfit(Y~F1+F2+F3+F4+Y_0+Y_1+Y_2+Y_3, Z, peninit=FALSE, postest=FALSE, ptvcste=ptvcste)
	
	lmZ <- lm(Z)
	rmses <- c(sqrt(var(Z[,1])), sqrt(mean(res$residuals^2)), sqrt(mean(lmZ$residuals^2)))
	rmses <- c(rmses, rmses[2]/rmses[3])
	allmse <- rbind(allmse, rmses)
	allres <- cbind(allres, res$residuals)
	allolsres <- cbind(allolsres,lmZ$residuals)
	
	rcoef <- ts(matrix(res$rawcoef,ncol=res$r,byrow=TRUE), start=c(1960,2), freq=4)
    if (ptvcste) {
        forecastintercepts <- cbind(forecastintercepts, as.numeric(rcoef[-1,1]!=0))
        rcoef <- rcoef[,-1]
    }

	# FacCoef Figures
	coefs <- ts(res$coef, start=c(1960,2), freq=4)
    if (ptvcste) coefs <- coefs[,-1]
	pdata <- data.frame(x=(time(coefs)), coefs)
	pdata <- melt(pdata, id.vars=c("x"));
	levels(pdata$variable) <- c(paste("Factor", 1:noF), paste("X Lag", 0:3))
	g <- ggplot(pdata) + geom_line(aes(x=x, y=value)) + theme_bw() + 
			facet_wrap(~variable, ncol=2, scale='free')+
			scale_x_continuous(breaks = seq(round(min(pdata$x)), round(max(pdata$x)), by = 5))
	g <- g + geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + ylab("") + xlab("") + ggtitle(strsplit(lna,',')[1][[1]])
    g <- g + theme(plot.title = element_text(vjust=2, size=12))
	pdf(paste("figures/FacCoefs-", mne, ".pdf", sep=""), width=10, height=6)
	print(g)
	dev.off()

	# All periods
	nz <- colSums(rcoef[-1,]!=0)
	res1 <- rbind(res1, nz)
	# 1970s to 1984
	nz <- colSums(window(rcoef, start=c(1970,1), end=c(1984,4))!=0)
	res2 <- rbind(res2, nz)
	# 1985 to 2000
	nz <- colSums(window(rcoef, start=c(1985,1), end=c(1999,4))!=0)
	res3 <- rbind(res3, nz)
	# All non zero
	allnz <- rbind(allnz,sum(rcoef[-1,]!=0))
}

#rmse plots
resstd <- apply(allres,2,function(x)sqrt(mean(x^2)))
olsstd <- apply(allolsres,2,function(x)sqrt(mean(x^2)))

Frmse <- function(x,resstd)(mean((x)^2))
x<-data.frame(cbind(time(rcoef),apply(allres,1,Frmse,resstd),apply(allolsres,1,Frmse,olsstd)))
colnames(x) <- c('Date','ptv','OLS')
mx <- melt(x,id.vars = 1)
grmse <- ggplot(mx,aes(x=Date,y=value,linetype=variable)) + geom_line() + theme_bw() + ylab("") +  theme(legend.title=element_blank())
print(grmse)
ggsave(filename = 'figures/grmse.pdf', plot = grmse, width = 10, height = 4)

rownames(res1) <- colnames(AllData)
rownames(res2) <- colnames(AllData)
rownames(res3) <- colnames(AllData)

AllRes <- cbind(
	rowSums(res1[,1:4]),
	rowSums(res2[,1:4]),
	rowSums(res3[,1:4]),
	rowSums(res1[,5:8]),
	rowSums(res2[,5:8]),
	rowSums(res3[,5:8]))
AllRes[AllRes == 0] <- NA
AllResFac[AllResFac == 0] <- NA

out <- cbind(vnames[,1], AllResFac, AllRes);
print(xtable(out), sanitize.rownames.function=function(x) {x}, 
	  only.contents=TRUE, hline.after=NULL, include.colnames=FALSE, include.rownames=FALSE, 
	  file="tables/nobreaks.tex", booktabs=TRUE)
print(xtable(out[aggidx,]), sanitize.rownames.function=function(x) {x}, 
	  only.contents=TRUE, hline.after=NULL, include.colnames=FALSE, include.rownames=FALSE, 
	  file="tables/nobreaks_agg.tex", booktabs=TRUE)
print(xtable(out[!aggidx,]), sanitize.rownames.function=function(x) {x}, 
	  only.contents=TRUE, hline.after=NULL, include.colnames=FALSE, include.rownames=FALSE, 
	  file="tables/nobreaks_subagg.tex", booktabs=TRUE)

mseout <- cbind(vnames[,1], formatC(allmse, digits=4, format="f"),allnz)
rownames(mseout) <- NULL
print(xtable(mseout), sanitize.rownames.function=function(x) {x}, 
	  only.contents=TRUE, hline.after=NULL, include.colnames=FALSE, include.rownames=FALSE, 
	  file="tables/mse.tex", booktabs=TRUE)
print(xtable(mseout[aggidx,]), sanitize.rownames.function=function(x) {x}, 
	  only.contents=TRUE, hline.after=NULL, include.colnames=FALSE, include.rownames=FALSE, 
	  file="tables/mse_agg.tex", booktabs=TRUE)
print(xtable(mseout[!aggidx,]), sanitize.rownames.function=function(x) {x}, 
	  only.contents=TRUE, hline.after=NULL, include.colnames=FALSE, include.rownames=FALSE, 
	  file="tables/mse_subagg.tex", booktabs=TRUE)

# Factor figures

a <- a2 <- c()
for (factorno in 1:noF)
{
	b <- b2 <- c()
	for (i in 1:190) b <- c(b, sum((PTVEst$RawCoef[,factorno,i]!=0)))
	for (i in 1:190) b2 <- c(b2, mean(abs(PTVEst$RawCoef[,factorno,i])))
	a <- cbind(a, b)
	a2 <- cbind(a2, b2)
}
b <- b2 <- c()
for (i in 1:190) b <- c(b, sum((PTVEst$RawCoef[,1:factorno,i]!=0)))
for (i in 1:190) b2 <- c(b2, mean(abs(PTVEst$RawCoef[,1:factorno,i])))

d <- ts(a, start=c(1959,3), freq=4)
colnames(d) <- paste("F", 1:noF, sep="")
d2 <- ts(b, start=c(1959,3), freq=4)
pdata <- data.frame(x=(time(d)), AllF=coredata(d2), coredata(d))
pdata <- melt(pdata, id.vars=c("x"));
levels(pdata$variable) <- c("All Factors", paste("Factor", 1:noF))
g <- ggplot(subset(pdata,variable!='All Factors')) + geom_bar(aes(x = x, y = value), stat = "identity") + theme_bw() + 
		facet_wrap(~variable,ncol=2,scales='free_y') +
		scale_x_continuous(breaks = seq(round(min(pdata$x)), round(max(pdata$x)), by = 5))
g <- g + geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + ylab("") + xlab("") + ggtitle('Number of breaks in the loadings of the factor model')
g <- g + theme(plot.title = element_text(vjust=2, size=12))
pdf(paste("figures/FactorFigAll.pdf", sep=""), width=10, height=7)
print(g)
dev.off()

d <- ts(a2, start=c(1959,3), freq=4)
colnames(d) <- paste("F", 1:noF, sep="")
d2 <- ts(b2, start=c(1959,3), freq=4)
pdata <- data.frame(x=(time(d)), AllF=coredata(d2), coredata(d))
pdata <- melt(pdata, id.vars=c("x"));
levels(pdata$variable) <- c("All Factors", paste("Factor", 1:noF))
g <- ggplot(subset(pdata,variable!='All Factors')) + geom_bar(aes(x = x, y = value), stat = "identity") + theme_bw() + 
		facet_wrap(~variable,ncol=2,scales='free_y') +
		scale_x_continuous(breaks = seq(round(min(pdata$x)), round(max(pdata$x)), by = 5)) 
g <- g + geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + ylab("") + xlab("") + ggtitle('Mean absolute change in the loadings of the factor model')
g <- g + theme(plot.title = element_text(vjust=2, size=12))
pdf(paste("figures/FactorFigAllChange.pdf", sep=""), width=10, height=6)
print(g)
dev.off()

# Loadings figures

factorno <- 1:noF
b <- c();
for (i in 1:144) b <- c(b, sum((PTVEst$RawCoef[i,factorno,]!=0)))
cbind(b, colnames(AllData), vnames)->d
for (ii in 1:nrow(d))
{
	varname <- d[ii, 2]
	mne <- d[ii, 4]
	lna <- vnames[ii, 5]
	colnames(AllData) == varname->idx
	loadings <- t(PTVEst$Lambda[idx,,])
	loadings <- ts(loadings, start=c(1959,3), freq=4)
	pdata <- data.frame(x=(time(loadings)), loadings)
	pdata <- melt(pdata, id.vars=c("x"));
	levels(pdata$variable) <- paste("Factor", 1:noF)
	g <- ggplot(pdata) + geom_line(aes(x=x, y=value)) + theme_bw() + 
				facet_wrap(~variable, scales="free",ncol=2)+
				scale_x_continuous(breaks = seq(round(min(pdata$x)), round(max(pdata$x)), by = 5))
	g <- g + geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
	g <- g + ylab("") + xlab("") + ggtitle(strsplit(lna,',')[1][[1]])
    g <- g + theme(plot.title = element_text(vjust=2, size=12))
	pdf(paste("figures/Loadings-", mne, ".pdf", sep=""), width=10, height=4)
	print(g)
	dev.off()
}


# Plotting the estimated factors
colnames(F) <- paste('Factor',1:4,sep=' ')
fdf <- data.frame(x=time(F),F)
fdf[,2] <- fdf[,2]*-1
mfac <- melt(fdf,id.vars=1)
levels(mfac$variable) <- paste("Factor", 1:4)

g <- ggplot(mfac) + geom_line(aes(x=x, y=value)) + theme_bw() + 
			facet_wrap(~variable, scales="free",ncol=2)+
			scale_x_continuous(breaks = seq(round(min(pdata$x)), round(max(pdata$x)), by = 5))
g <- g + geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
g <- g + ylab("") + xlab("") + ggtitle('Estimated Factors')
g <- g + theme(plot.title = element_text(vjust=2, size=12))

pdf("figures/EstimatedFactors.pdf", width=10, height=4)
print(g)
dev.off()

# Time-varying intercepts
if (ptvcste) {
    iforecast <- ts(rowSums(forecastintercepts), start=c(1960,3), freq=4)
    ifactor <- ts(PTVEst$intercepts, start=c(1959, 4), freq=4)
    d <- cbind(ifactor, iforecast) 
    colnames(d) <- c("Factor Model", "Forecasting Model")
    pdata <- data.frame(x=(time(d)), coredata(d))
    pdata <- melt(pdata, id.vars=c("x"));
    levels(pdata$variable) <- c("Factor Model", "Forecasting Model")
    g <- ggplot(pdata) + geom_bar(aes(x = x, y = value), stat = "identity") + theme_bw() + 
    		facet_wrap(~variable,ncol=1,scales='free_y') +
    		scale_x_continuous(breaks = seq(round(min(pdata$x)), round(max(pdata$x)), by = 5))
    g <- g + geom_rect(aes(xmin=1960.25, xmax=1961, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + geom_rect(aes(xmin=1969.75, xmax=1970.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + geom_rect(aes(xmin=1973.75, xmax=1975.0, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + geom_rect(aes(xmin=1980, xmax=1980.5, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + geom_rect(aes(xmin=1981.5, xmax=1982.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + geom_rect(aes(xmin=1990.5, xmax=1991, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + geom_rect(aes(xmin=2001, xmax=2001.75, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.01)
    g <- g + ylab("") + xlab("") + ggtitle('Number of breaks in the intercept of the models')
    g <- g + theme(plot.title = element_text(vjust=2, size=12))
    pdf(paste("figures/PtvIntercepts.pdf", sep=""), width=10, height=7)
    print(g)
    dev.off()
}

