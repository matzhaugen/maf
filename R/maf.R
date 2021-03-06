#' Maximum Autocorrelation Factors
#' 
#' This function performs a linear transformation of the data where the new time series
#' are mutually orthogonal with identity covariance matrix.
#' The input is a matrix of column vector time series.
#' @param x A matrix where each column is a time series. The number of 
#' rows must be greater than the number of columns
#' @return An object of class "Maf" containing the following items:
#' \describe{
#' 	\item{\code{x}}{The original data matrix.}
#' 	\item{\code{autocor}}{A vector of autocorrelations for each maf factor.}
#'  \item{\code{rotation}}{A matrix of rotation vectors, each column representing a vector, 
#' 	with the first column containing the first rotation (MAF1), the second representing
#'  the second rotation and so on. Each column's squares sum to 1.}
#'  \item{\code{mafs}}{Contains the matrix with maf factors in the 
#' 	columns, with highest autocorrelation in the first column.}
#' }
#' @examples
#' # Extract mafs from dataset
#' maf.object = maf(treeringTimeseries)
#' # Plot the first 6 mafs with undertainty estimates and estimate number of 
#' # significant mafs contained in the dataset
#' plot(maf.object) 
#-------------
#' @export
maf <- function(x) {

	if (class(x)=="data.frame" || (class(x)=="matrix" && dim(x)[2]>1)){
		p = dim(x)[2]
		n = dim(x)[1]
		x = scale(x)
		svd = svd(cov(x)) #If you want to sum of the timesteps/n to be 1
		#svd = svd(cov(x)*(n-1)) #If you want to sum of the timesteps to be 1

		a = svd$u%*%diag(svd$d^(-0.5))%*%t(svd$u)

		y = x%*%a
		cov.zd = cov(apply(y,2,diff))*(n-2)
		
		svd.zd = svd(cov.zd/(n-1))
		u = svd.zd$u[,p:1]
		aa = a%*%u
		aa = apply(aa,2,function(x) {x/sqrt(sum(x^2))})

		maf = x%*%aa
		neg.time.cor = diag(cov(maf,matrix(rep(1:n,p),n,p)))<0
		mafs = t(apply(maf,1,function(x) {(-neg.time.cor+!neg.time.cor)*x}))
		aa = t(apply(aa, 1,function(x) {(-neg.time.cor+!neg.time.cor)*x}))

		d = svd.zd$d[p:1]
	} else {
		print("x is not a matrix or a data.frame, no MAF transform performed.")
		mafs = as.matrix(x)
		aa = as.matrix(1); a.i = as.matrix(1); d = as.matrix(1)
	}
	out = list(x=x, mafs=maf, rotation=aa, autocor=1 - d/2)
	class(out) = c("Maf")
	out
}

#This function does what "filter" does for 2 sides except that
#the vector length is preserved by renormalizing the filter window with
#the available timesteps. f must be odd length.
#If f is a scalar, the value is assumed to be the number of years per degree of
#freedom. E.g. if f = 10, a 100 length vector will have 10 DoF and 11 knots.
filter2 <- function(x,f) {
	
	x = as.matrix(x)
	l = dim(x)[1]
	p = dim(x)[2]
	if (length(f) == 1) {
		if (f == 1 || f == 0) {
			return(x)
		} else {
			if (f < 1 && f > 0) {
				return(apply(as.matrix(x),2,function(y) {
					loess(y~c(1:l), span=f)$fitted
				}))
			} else {
				return(apply(as.matrix(x),2,function(y) {
					loess(y~c(1:l), span=2*f/l)$fitted
				}))
			}
		}
	}
}

# Function: blockBs(x, block.size, circular = F)
# ---------------------------------
# This function will take in a matrix of size n-by-p and return a resampled version
# with the same predictor order and size but permuted rows. The permutation will be in
# blocks of size "block.size". The circular argument refers to weather if a block randomly
# starts after n-blocksize, in which case the block would have to continue at the
# beginning in order to be complete. Otherwise, such selections would not be possible.
# 7/22/14: Cannot handle d>0 yet.
blockBs <- function(x, blocksize, circular = T, d=0) {
	x = as.matrix(x)
	n = dim(x)[1]
	p = dim(x)[2]
	length = 0
	indeces = NULL
	Ex = colMeans(x)
	if (!d==0) {
		dx = x-x[1,]
		for (i in 1:d) {
			dx = as.matrix(apply(dx, 2, diff))
		}
	}

	#Getting new indeces
	if (!circular) {
		while (length(indeces)<n) {
			startInd = sample(1:(n-blocksize+1-d),1)
			indeces = c(indeces,startInd:(startInd+blocksize-1))
		}
	} else {
		while (length(indeces)<(n-d)) {
			startInd = sample(1:(n-d),1)
			v = startInd:(startInd+blocksize-1)
			if ((startInd+blocksize-1)>(n-d)) {
				indeces = c(indeces,((v-1) %% (n-d))+1) # makes the indeces larger than (n-d)
					# start from the beginning
			} else {
				indeces = c(indeces,v)
			}
		}
	}

	if (!d==0) {
		print(dim(dx))
		print(length(indeces[1:(n-d)]))
		dxstar = as.matrix(dx[indeces[1:(n-d)],])
		ixstar = dxstar
		for (i in 1:d) {
			ixstar = x + dxstar
			xstar = x[1,] + rbind(0, ixstar)
		}

	} else {
		xstar = x[indeces[1:(n-d)],]
	}
	xstar
}

getSnrEmpir <- function(x, f) {
	if (class(x) != "matrix") {
			out = getSnrEmpirSub(x, f)
		} else {
			out = apply(x, 2, getSnrEmpirSub, f)
		}
	out
}
#Get the empirical snr based on filter size f and time series input x.
#outputs a scalar.
getSnrEmpirSub <- function(x, f, plotit=F, ...) {

	signal = filter2(x, f)
	noise = x - signal
	snr = var(signal)/var(noise)
	if (plotit) {
		plot(x, type="l", ...)
		lines(signal, col=3, lwd=2)
	}
	snr
}

# Get the sign of the max of the absolute value of each vector entry
vectorAbsMaxSign <- function(v) {
	output = rep(0, length(v))
	maxAbsV = max(abs(v))
	maxV = max(v)
	minV = min(v)
	if (maxAbsV == maxV) {
		output[which.max(v)] = 1
	} else if (maxAbsV == -minV) {
		output[which.min(v)] = -1
	} 
	output
}

#' Test number of significant maf time series in data
#' 
#' Perform a statistical test for how many signals are hiding in the data. 
#' The output gives p-values for each time series
#' @param maf.object The output of the \code{maf} function.
#' @param smooth.span Fraction between 0-1 that specifies the proportion of timesteps
#' to include in the smoothing window, which is weighted by a tricubic. See ?loess for
#' details. Alternatively, one can specify an integer greater than 1 which refers to the 
#' number of time points to include in the tricubic filter.
#' @param B Number of replications in the confidence interval calculations
#' @param alpha The significance level
#' @param block.size The block size of the resampled residuals, i.e. the number of 
#' contiguous time steps to sample at the time from the set of residuals.
#' @return A list with the following items
#' \describe{
#'  \item{\code{statStar}}{ Empirical SNR for each resampled MAF time series, i.e. a p x B matrix 
#'  where p is the total number of predictors and B is the number of bootstraps/replications}
#'  \item{\code{statObs}}{ The empirical SNR of the original MAF time series.}
#'  \item{\code{pval}}{The p-values of each MAF, where the p-value refers to the number of 
#' resampled MAFs that have a higher empirical SNR than the original MAFs. If there are less
#' than \eqn{\alpha}, e.g. 0.05, resampled MAFs that have a higher empirical SNR than the
#' original MAF, the MAF in question is not significant.}
#' \item{\code{mafStar}}{ The B resampled maf timeseries, in a 3d array.}
#' \item{\code{mafSmooth}}{ The smooth original maf estimates.}
#' }
#' @details 
#' The confidence intervals are obtained using a resampling scheme which extracts the residuals of the original 
#' time series after smoothing
#' with the same filter as the one parametrized by \code{smooth.span}. The residuals are 
#' bootstrapped or block bootstrapped (if there is temporal structure in the residuals)
#' and added back to the smooth, creating a new data set. From the new data set a new set of
#' MAFs are calculated.
#' 
#' The smoothing parameter \code{smooth.span} is also used when calculating the empirical 
#' Signal-to-Noise ratio. Where we subtract the smooth from the time series and calculate
#' the variance of the smooth estimate over the residual variance. This gives an empirical
#' signal-to-noise estimate.
#' 	
#' @examples
#' # Extract mafs from dataset
#' maf.object = maf(treeringTimeseries)
#' # Test for how many mafs are present
#' test.Maf(maf.object) 
#-------------

# There are two methods:
# method 1: Do subtract the smooth and calculate test statistic
# p val = number of test stats that are greater than original
# method 2: Do NOT subtract the smooth and calculate test statistic
#' @export
test.Maf <- function(maf.obj, alpha=0.05, block.size=5, smooth.span=30, B=100) {
	
	mafs = (maf.obj$mafs)
	n = dim(mafs)[1]
	p = dim(mafs)[2]	
	mafStar = array(0, dim=c(B, n, p))
	statStar = matrix(0, p, B)
	x = maf.obj$x		
	mafSmooth = filter2(mafs, smooth.span)[, 1:p]
	xSmooth = filter2(x, smooth.span)
	residuals = x - xSmooth
	statObs = sort(getSnrEmpir(mafs, smooth.span), decreasing=T)
	for (i in 1:B){
		residualsStar = blockBs(residuals, block.size)
		# xStar = residualsStar # method 1
		xStar = xSmooth + residualsStar # method 2
		mafStarUnsigned = as.matrix(maf(xStar)$mafs[, 1:p])
		corMat = cor(mafStarUnsigned, mafSmooth)
		permutation = t(apply(corMat, 1, vectorAbsMaxSign))
		mafStar[i,,] = mafStarUnsigned %*% permutation
		statStar[, i] = sort(getSnrEmpir(maf(residualsStar)$mafs[, 1:p], smooth.span), decreasing=T)
	}

	# calculate p-values
	j = 0
	pval = apply(statStar, 1, function(x) {
		j <<- j + 1
		sum(x > statObs[j]) / B # method 1
		# sum(x > statObs[j]) / B  # method 2
	})
	if (class(statStar) != "matrix") {
		statStar = matrix(statStar, 1, B)
	}
	
	if (all(pval < alpha)) {
		nSignificantMafs = p
	} else {
		nSignificantMafs = min(which(pval > alpha)) - 1 # method 1
		# nSignificantMafs = min(which(pval > alpha)) # method 2
	}
	print(paste("Estimate of # of MAFs in dataset: ", nSignificantMafs))
	output = list(pval=pval, statStar=statStar, statObs=statObs, mafStar=mafStar, mafSmooth=mafSmooth)
}


#' Plot maf object
#' 
#' This function will plot the first few mafs in the dataset overlaid
#' with a smooth version of the same maf and, if desired, confidence intervals. The
#' legend shows the emirical signal-to-noise ratio, i.e. signal variance (across time) divided by noise variance
#' Additionally, an estimate of the number of mafs will be printed.
#' @param maf.object The output of the \code{maf} function.
#' @param smoothSpan Fraction between 0-1 that specifies the proportion of timesteps
#' to include in the smoothing window, which is weighted by a tricubic. See ?loess for
#' details. Alternatively, one can specify an integer greater than 1 which refers to the 
#' number of time points to include in the tricubic filter.
#' @param cexVal The size of the labels in the plot
#' @param nmaf The number of MAFs to plot
#' @param with.wncertainty A logical specifying whether to plot the confidence interval
#' around each maf.
#' @param B Number of replications in the confidence interval calculations
#' @param alpha The significance level
#' @param block.size The block size of the resampled residuals, i.e. the number of 
#' contiguous time steps to sample at the time from the set of residuals.
#' @return A plot of the MAFs and a list containing 
#' \describe{
#'  \item{\code{statStar}}{ Empirical SNR for each resampled MAF time series, i.e. a p x B matrix 
#'  where p is the total number of predictors and B is the number of bootstraps/replications}
#'  \item{\code{statObs}}{ The empirical SNR of the original MAF time series.}
#'  \item{\code{pval}}{The p-values of each MAF, where the p-value refers to the number of 
#' resampled MAFs that have a higher empirical SNR than the original MAFs. If there are less
#' than \eqn{\alpha}, e.g. 0.05, resampled MAFs that have a higher empirical SNR than the
#' original MAF, the MAF in question is not significant.}
#' }
#' @details 
#' The confidence intervals are obtained using a resampling scheme which extracts the residuals of the original 
#' time series after smoothing
#' with the same filter as the one parametrized by \code{smoothSpan}. The residuals are 
#' bootstrapped or block bootstrapped (if there is temporal structure in the residuals)
#' and added back to the smooth, creating a new data set. From the new data set a new set of
#' MAFs are calculated.
#' 
#' The smoothing parameter \code{smooth.span} is also used when calculating the empirical 
#' Signal-to-Noise ratio. Where we subtract the smooth from the time series and calculate
#' the variance of the smooth estimate over the residual variance. This gives an empirical
#' signal-to-noise estimate.
#' 	
#' @examples
#' # Extract mafs from dataset
#' maf.object = maf(treeringTimeseries)
#' # Plot the first 6 mafs with undertainty estimates and estimate number of 
#' # significant mafs contained in the dataset
#' plot(maf.object) 
#-------------
#' @export
plot.Maf <- function(maf.obj, smooth.span=30, cexVal=1.5, 
					 which.maf=c(1:min(dim(maf.obj$x)[2], 3)), 
					 with.uncertainty=TRUE, B=100, alpha=0.05, 
					 block.size=5) {
	
	mafs = maf.obj$mafs
	n = dim(mafs)[1]
	p = dim(mafs)[2]		
	
	if (with.uncertainty) {
		output = test.Maf(maf.obj, alpha, block.size, smooth.span, B)		
		mafSmooth = output$mafSmooth
	} else {
		output = NULL
		mafSmooth = filter2(mafs, smooth.span)[, 1:p]
	}

	snr = round(getSnrEmpir(mafs, smooth.span), digits=3)

	# Plot
	nmaf = length(which.maf)
	if (nmaf < 3) {
			par(mfrow=c(1, nmaf))	
	} else {
		par(mfrow=c(ceiling(nmaf / 3), 3))
	}
	for (i in which.maf) {
		par(mar=c(4,4,4,2))
		plot(mafs[, i], cex.axis=cexVal, cex.main=cexVal, main=paste("MAF", i),
			cex.lab=cexVal, ylim=range(c(mafs[, which.maf])), col=gray(0.5), type="l",
			ylab="", xlab="Time step")
		lines(mafSmooth[, i], lwd=3, col=1)
		if (with.uncertainty) {
			mafStar = output$mafStar
			upper = apply(mafStar[,,i], 2, quantile, 1-alpha/2)
			lower = apply(mafStar[,,i], 2, quantile, alpha/2)	
			lines(upper, lty=2, lwd=1)
			lines(lower, lty=2, lwd=1)			
		}
		legend('topleft', c(paste("SNR:", snr[i])), lwd=2)
	}
}
#' High dimensional MAF extraction
#' 
#' This is a generalization of the MAF algorithm when the sample size (length of time series) 
#' is smaller than the number of time series, i.e. n < p. The idea is to first partition the set of time series 
#' into smaller groups such that the n > p and the problem is well-posed with a non-singular covariance
#' matrix. From each group we extract the first MAF and then 
#' @param mat A matrix of dimension n x p. Tested only for dimensions smaller than n/p~10000.
#' @param alpha Tuning parameter \in [0, 1] that describes the size of the partitioning groups. Viz, 
#' group size <= n * alpha. We find for the cases we tried that the results are not that dependant 
#' on the value of \alpha as long as it's in the range of 0.2 to 0.8.
#' @param type What type of partition to use (default is 'random'). We find that a random partition yields
#' the best results. There are theoretical justifcations for this in an upcoming paper.
#'
#' @return A list containing the following items
#' \describe{
#'  \item{\code{out}}{MAF1 time series}
#'  \item{\code{a}}{The set of weights which when multiplied by the original time series produces MAF1.}
#' }
#' @details 
#' Currently only implemented for the first MAF, for higher MAFs either stay tuned or subtract out MAF1 from the data set,
#' e.g. by linear regression, and then run the algorithm on the residuals.
#' 	
#' @examples
#' # Make a signal
#' n = 150 # number of time steps
#' p = 600 # number of time series
#' rho = 0.4 # noise cross correlation
#' cov_e = matrix(rho, p, p) # noise covariance
#' diag(cov_e) = 1
#' noise = mvrnorm(n=n, mu=rep(0, p), Sigma=cov_e)
#' x=seq(0,4,len=n);
#' signal = (x-1)*(x-2)*(x-3)
#' signal.sd = c(rep(0.2, p/2), rep(0.6, p/2))
#' signals = signal %*% t(signal.sd) / sd(signal)
#' data = signals + noise
#' result = MAFpart(data)
#' plot(result$out, type="l", xlab='Time steps', ylab='Signal (Unitless)', 
#' 	main='High Dimensional Signal Extraction')
#' lines(lm(result$out ~ signal)$fit, lwd=2, col=2)
#' legend('topleft', c('Estimate', 'Truth'), col=c(1,2), lwd=c(1,2))
#-------------
#' @export
MAFpart <- function(mat, alpha=0.2, type="random") {
	p = dim(mat)[2]
	n = dim(mat)[1]
	if (is.vector(mat)) {
		out = mat
		a = 1
	} else if (p<=(n * alpha)) {
		maf = MAF(mat)
		out = as.matrix(maf$out[,1], n, 1)
		a = maf$a[,1]
	} else {
		m = ceil(p/(alpha * n))	
		if (type=="random") {
			v = rep(1:m, ceil(p / m))
			ind = v[sample(v, p)] # shuffle randomly!
		} else if (type=="kmeans") {
			ind = kmeans(cov(mat), m)$clu			
		} else if (type=="kmediod") {
			ind = pam(cov(mat), m, cluster.only=TRUE)
		} else {
			print("not a valid type")
		}		

		myList = list(); length(myList) = m
		out = list(); length(myList) = m
		v = list(); length(myList) = m
		for (i in 1:m) {
			myList[[i]] = mat[, ind == i] 
			maf = MAFpart(myList[[i]], type=type)
			out[[i]] = maf$out
			v[[i]] = maf$a
		}

		new_mat = do.call(cbind, out)
		mafrec = MAFpart(new_mat, type=type)
		out = mafrec$out
		w = mafrec$a
		a = rep(0, p)
		for (i in 1:m) {
			a[ind == i] = v[[i]] * w[i]
		}
		a = a / sqrt(sum(a^2))
	}
	list(out=out, a=a)
}

