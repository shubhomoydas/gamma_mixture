
#--------------------------------------
# NOTE: Make sure that the samples are never ZERO, 
# else, assign a very small value
#--------------------------------------

#######################################################
#
# samples - n x D matrix of numbers where:
#       n - number of instances
#       D - number of algorithms in the ensemble
#   
#######################################################

fnRankAggregationGamma <- function (samples, isDebug=F) {

	n <- nrow(samples)
	D <- ncol(samples)

	converged = F

	# ========== Initialize the Gamma priors ==========
	r <- matrix(100,nrow=2,ncol=D) # number of samples to use for sum(log(samples))
	r[1, ] <- round(n * 0.05)
	s <- matrix(100,nrow=2,ncol=D) # number of samples to use for sum(samples)
	p <- matrix(0,nrow=2,ncol=D) # stores sum(log(samples))
	q <- matrix(0,nrow=2,ncol=D) # stores sum(samples)
	# Process each detector d = 1..D
	for (d in 1:D) {
	  tmp_samples_d <- samples[,d]
	  Xsorted <- tmp_samples_d[order(tmp_samples_d)]
	  # For the variables p,q,r,s the first row corresponds 
	  # to anomalous data parameters and the second row corresponds
	  # to regular data parameters
	  for (i in 1:2) {
	    if (i==1) {
	      # Compute prior for anomalous data instances
	      Xtmp <- Xsorted[(n-r[i,d]+1):n]
	    } else {
	      # Compute prior for regular data instances
	      Xtmp <- Xsorted[((n/2)-(s[i,d]/2)):((n/2)+(s[i,d]/2)-1)]
	    }
	    p[i,d] <- sum(log(Xtmp))
	    s[i,d] <- length(Xtmp) # re-adjust s
	    q[i,d] <- sum(Xtmp)
	  }
	  # Following is for Debug purpose
	  #X_95 <- quantile(samples[,d],probs=c(0.95,1.00))
	  #X_50 <- quantile(samples[,d],probs=c(0.50,0.75))
	}
	# Following are fixed hyper-parameters for Beta prior of Pi
	alpha <- 2
	beta <- 100

	# Following are parameters that need to be inferred
	# The first row is for anomalies and second for regular
	alpha_kd <- matrix(1,nrow=2,ncol=D) # alpha_kd MUST be strictly greater than zero
	# We set initial values of anomaly alphas to have larger values than regular alphas
	alpha_kd[1,] <- 1.0
	alpha_kd[2,] <- 0.01
	beta_kd <- matrix(0,nrow=2,ncol=D)

	# f_1 is the first derivative of log-likelihood
	# f_2 is the second derivative of log-likelihood
	f_1 <- matrix(0,nrow=2,ncol=D)
	f_2 <- matrix(0,nrow=2,ncol=D)

	# Initialize Z to random. Z is latent variable that indicates
	# whether an instance is anomaly or regular
	Z <- matrix(0,nrow=n,ncol=2)
	Z[,1] <- runif(n)
	Z[,2] <- 1-Z[,1]

	log_f <- matrix(0,nrow=n,ncol=2)

	z_sums <- matrix(0, nrow=2, ncol=1)
	z_x_sums <- matrix(0, nrow=2, ncol=D)
	z_logx_sums <- matrix(0, nrow=2, ncol=D)

	# ========== Compute MLEs of Parameters ==========

	MaxEpochs <- 500 # ideally we should only run till convergence
	prev_ll <- -.Machine$double.xmax
	for (epoch in 1:MaxEpochs) {
	  
	  z_sums[,1] <- apply(Z,2,sum)
	  
	  for (d in 1:D) {
	    z_x_sums[1,d] <- sum(Z[,1]*samples[,d]) # anomaly proportions
	    z_x_sums[2,d] <- sum(Z[,2]*samples[,d]) # normal  proportions
	    z_x_sums[,d] <- z_x_sums[,d] + q[,d]
	    z_logx_sums[1,d] <- sum(Z[,1]*log(samples[,d])) # anomaly proportions
	    z_logx_sums[2,d] <- sum(Z[,2]*log(samples[,d])) # normal  proportions
	    z_logx_sums[,d] <- z_logx_sums[,d] + p[,d]
	  }
	  
	  # M-Step
	  # ===========
	  
	  # Compute MLE of Pi
	  Pi <- (z_sums[1,] + alpha-1)/(n + alpha-1 + beta-1)
	  
	  # Compute MLE of alpha_kd
	  for (d in 1:D) {
	    z_sums_ds <- s[,d]+z_sums
	    z_sums_dr <- r[,d]+z_sums
	    
	    # Compute MLE of alpha_kd by Minka's method
	    # This method is again an iterative method
	    for (iter in 1:50) {
	      f_1[,d] <- z_logx_sums[,d] + 
	        z_sums_ds*log(alpha_kd[,d]) - z_sums_ds*log(z_x_sums[,d]/z_sums_ds) -
	        z_sums_dr*digamma(alpha_kd[,d])
	      f_2[,d] <- z_sums_ds*(1/alpha_kd[,d]) - z_sums_dr*trigamma(alpha_kd[,d])
	      tmp_alpha_hat <- (1/alpha_kd[,d]) + (1/(alpha_kd[,d]*alpha_kd[,d]))*(f_1[,d]/f_2[,d])
	      tmp_alpha_hat <- 1/tmp_alpha_hat
	      diff_alpha <- alpha_kd[,d] - tmp_alpha_hat
	      dist_alpha <- sum(diff_alpha*diff_alpha)
	      if (min(tmp_alpha_hat) < 0) {
	        print(tmp_alpha_hat)
	        stop("tmp_alpha_hat is negaive!")
	      }
	      alpha_kd[,d] <- tmp_alpha_hat
	      #print(alpha_kd)
	      if (dist_alpha < 1e-8) {
	        #cat("alpha converged at iter: ",iter,"\n")
	        break;
	      } else {
	        #cat("iter: ",iter,", dist: ",dist_alpha,"\n")
	      }
	    }
	    
	    # Compute MLE of beta_kd
	    beta_kd[,d] <- alpha_kd[,d]*(z_sums_ds/z_x_sums[,d])
	    
	  }
	  
	  # E-Step
	  # ===========
	  
	  log_f[,1:2] <- 0
	  for (d in 1:D) {
	    log_f[,1] <- log_f[,1] + log(apply(as.matrix(dgamma(samples[,d], shape=alpha_kd[1,d], rate=beta_kd[1,d])),1,max,1e-320))
	    log_f[,2] <- log_f[,2] + log(apply(as.matrix(dgamma(samples[,d], shape=alpha_kd[2,d], rate=beta_kd[2,d])),1,max,1e-320))
	  }
	  log_f[,1] <- (log(Pi) + log_f[,1])
	  log_f[,2] <- (log(1-Pi) + log_f[,2])
	  
	  # This part of the log-likelihood should monotonously increase
	  ll <- sum(Z[,1]*log_f[,1] + Z[,2]*log_f[,2]) + (alpha-1)*log(Pi) + (beta-1)*log(1-Pi)
	  
	  for (i in 1:n) {
	    log_f[i,] <- log_f[i,] - mean(log_f[i,])
	    mxi <- max(log_f[i,])
	    if (mxi > 700) {
	      # handle numerical issues
	      log_f[i,] <- log_f[i,] - (mxi - 700)
	    }
	    Z[i,] <- exp(log_f[i,])
	    Z[i,] <- Z[i,] / sum(Z[i,])
	  }
	  
	  diff_ll <- abs(ll - prev_ll)
	  if (isDebug == T) cat("Epoch",epoch,": Loglikelihood: ",ll, "; change: ",diff_ll,"\n")
	  prev_ll <- ll
	  if (diff_ll < 1e-4) {
	    if (isDebug == T) cat("Converged...\n")
	    converged = T
	    break;
	  }
	  
	}
	
	# Print inferred parameters...
	#Pi
	#alpha_kd
	#beta_kd
	params <- list(pi=Pi, alphaKd=alpha_kd, betaKd=beta_kd, z=Z, epochs=epoch, converged=converged, diffLoglik=diff_ll)

}

