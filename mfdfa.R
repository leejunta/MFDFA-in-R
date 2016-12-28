#MFDFA
#Multifractal Detrended Fluctuation Analysis
#Grinnell College Psychology Department
#Please contact Jun Taek Lee at <leejunta@grinnell.edu> for more information

mfdfa <- function(x,s,q,m,R2) {

#x: data
#s: vector of bin sizes (min:max)
#q: q-order (min:max)
#m: polynomial order
#R2: coefficient of determination parameter
#r2: coefficient of determination of q's
#segments: number of bins
#N: length of data
#xmean: mean of data
#RMS: root mean square
#qRMS: (root mean square)**q
#Fq: qth order fluctuation function, Fq(s)
#hq: Hurst exponent
#tq: tau(q)
#a: Singularity strength/Holder exponent (alpha)
#fa: dimension of series characterized by alpha

options(max.print=1000000)
segments <- vector()
fa_ai <- data.frame()
RMS <- data.frame()
qRMS <- data.frame()
Fq <- data.frame()
log10Fq <- list()
log10s <- list()
hq <- vector()
tq <- vector()
a <- vector()
fa <- vector()
r2 <- vector()

#Standardize data with mean=10, sd=1
N <- length(x)
xmean <- mean(x)
xstd <- sd(x)
x<-(x-xmean)/xstd+10

if (min(s)<(m+2)) stop("s<m+2",call.=FALSE)
if (max(s)>N/4) stop("s>N/4",call.=FALSE)

Y <- as.vector(rep(0,N))
for (i in 1:N) {
	for (k in 1:i) {
		Y[i] <- Y[i]+(x[k]-xmean)
	}
}

for (ns in 1:length(s)) {
	segments[ns] <- floor(N/s[ns])
	for (v in 1: segments[ns]) {
		index <- c(((v-1)*s[ns]+1):(v*s[ns]))
		linearmodel <- lm(Y[index]~poly(index,m,raw=TRUE))
		yv <- linearmodel$fitted.values
		RMS[ns,v] <- sqrt(sum((Y[index]-yv[1:(max(index)+1-min(index))])**2)/s[ns])
	}
	for (nq in 1:length(q)) {
		if (q[nq]==0) {
			qRMS[nq,ns] <- sum(log(RMS[ns,1:segments[ns]]**2))
			Fq[nq,ns] <- exp(0.5*qRMS[nq,ns]/segments[ns])
		} else {
			qRMS[nq,ns] <- sum((RMS[ns,1:segments[ns]]**(q[nq])))
			Fq[nq,ns] <- (qRMS[nq,ns]/segments[ns])**(1/q[nq])
		}
	}
}

nq <- c(1:length(q))
log10s <- lapply(nq, function(nq) log10(s)) 
log10Fq <- lapply(nq, function(nq) as.vector(t(log10(Fq[nq,]))))

for (i in 1:length(q)) {
	inf <- which(is.infinite(log10Fq[[i]]))
	ii <- 0
	repeat {
		ii <- ii + 1
		temp <- inf[ii]
		if (is.na(temp)) break
		log10Fq[[i]][temp] <- NA
		log10s[[i]][temp] <- NA
	}
	log10Fq[[i]] <- log10Fq[[i]][!is.na(log10Fq[[i]])]
	log10s[[i]] <- log10s[[i]][!is.na(log10s[[i]])]
}

tempFq <- log10Fq
temps <- log10s
tempr2 <- vector()

for (nq in 1:length(q)) {
	i <- 0
	repeat {
		i <- i + 1
		rsquared <- summary(lm(tempFq[[nq]]~temps[[nq]]))
		length <- length(tempFq[[nq]][!is.na(tempFq[[nq]])])
		tempr2[i] <- rsquared$r.squared
		if (length < length(s)/4) {
			ii <- which.max(tempr2)
			log10Fq[[nq]] <- head(log10Fq[[nq]],(length(log10Fq[[nq]])-ii))
			log10s[[nq]] <- head(log10s[[nq]],(length(log10s[[nq]])-ii))
			r2[nq] <- max(tempr2)
			break
		}
		if (rsquared$r.squared < R2) {
			tempFq[[nq]][length] <- NA
			temps[[nq]][length] <- NA
		} else {
			r2[nq] <- tempr2[i]
			log10Fq[[nq]] <- tempFq[[nq]][!is.na(tempFq[[nq]])]
			log10s[[nq]] <- temps[[nq]][!is.na(temps[[nq]])]
			break
		}
	}	
}

for (nq in 1:length(q)) {
	linearmodel <- lm(log10Fq[[nq]]~log10s[[nq]])
	hq[nq] <- linearmodel$coefficients[2]
	tq[nq] <- hq[nq]*q[nq]-1
}

for (nq in 1:length(q)-1) {
	a[nq] <- (tq[nq+1]-tq[nq])/(q[nq+1]-q[nq])
	fa[nq] <- q[nq]*a[nq]-tq[nq]
}

fa_ai <- rbind(fa_ai,a)
fa_ai <- rbind(fa_ai,fa)
fa_ai <- t(fa_ai)
colnames(fa_ai) <- c("alpha","f(alpha)")
return(fa_ai) 
}
