source("./scripts/functions/fdexgun.R") # the disfit routines

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
RobustEZ.from.File = function(file.name)
{
#========================================================================== 
# Read in data from file.name
#==========================================================================
	 
  d    = scan(file = file.name)
  N    = d[1]
  rt   = d[2:length(d)]
  ac   = length(rt)/N

  pars = Get.Robust.vaTer(rt, ac, min_U=min(rt), max_U=max(rt),start_p_EG=.95)
 
	return(pars)
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Get.Robust.vaTer = function(dat, Pc, min_U=min(dat), max_U=max(dat), start_p_EG=.95, s=.1)
{
#========================================================================== 
# Get cleaned-up values for EZ pars
#==========================================================================
	meanvarp = Get.Robust.MV(dat, min_U, max_U, start_p_EG)
	MRT      = meanvarp[1]
	VRT      = meanvarp[2]
	p_EG	   = meanvarp[3]
	# MRT & VRT cleaned up, now plug these in to EZ equations
	s2 = s^2
	# The default value for the scaling parameter s equals .1
	if (Pc == 0)
		cat("Oops, Pc == 0!\n")
	if (Pc == 0.5)
		cat("Oops, Pc == .5!\n")
	if (Pc == 1)
		cat("Oops, Pc == 1!\n")
	# If Pc equals 0, .5, or 1, the method will not work, and
	# an edge-correction is required.
	L = qlogis(Pc)
	# The function "qlogis" calculates the logit.
	x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
	v = sign(Pc-0.5)*s*x^(1/4)
	# This gives drift rate.
	a = s2*qlogis(Pc)/v
	# This gives boundary separation.
	y = -v*a/s2
	MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
	Ter = MRT-MDT
	# This gives nondecision time.
	return(c(v, a, Ter, p_EG))
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Get.Robust.MV = function(dat, min_U, max_U, start_p_EG)
{
#========================================================================== 
# Get cleaned-up values for RT mean and RT variance
#==========================================================================
	pars     = Est.ExGaussUnif(dat, min_U, max_U, start_p_EG)
	EG_mean  = pars[1] + pars[3]
	EG_var   = pars[2]^2 + pars[3]^2
	meanvarp = c(EG_mean, EG_var, pars[4])
	return(meanvarp)
}
# to check:
# mu=.4; sigma=.035; tau=.065; p_EG=.8; min_U=0.1; max_U=1.2
# dat  = Gen.ExGaussUnif(1000, mu, sigma, tau, p_EG, min_U, max_U)
# Est.ExGaussUnif(dat,.1,1.2,.9)
# Get.Robust.MV(dat, .1,1.2,.9)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Est.ExGaussUnif = function(dat, min_U, max_U, start_p_EG)
{
#========================================================================== 
# Get Estimates for Ex-Gauss/Unif RT distribution
#==========================================================================
	start.pars = c(MMest.ExGauss(dat), start_p_EG)

	lo   = c(0.1, .010, .010, 0.60)      # lower bounds; note bound on mix proportion 
	up   = c(1,   .300, .300, 0.9999999) # upper bounds
	y    = c(min_U,max_U) # fixed unif pars

# quasi newton
  res = try(optim(start.pars, loglexgm, loglexggm, method="L-BFGS-B",
	            lower=lo, upper=up, hessian=T, x=dat, y=y), TRUE)

  if (class(res) == "try-error") #optimization failed
	{
	  pars = c(NA, NA, NA, NA)   
	}

  if (class(res) != "try-error") #optimization worked
	{
    pars = res$par[1:4]
	}

	return(pars)
}
# to check:
# mu=.4; sigma=.035; tau=.065; p_EG=.8; min_U=.1; max_U=1.2
# dat = Gen.ExGaussUnif(1000, mu, sigma, tau, p_EG, min_U, max_U)
# Est.ExGaussUnif(dat,.1,1.2,.9)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
MMest.ExGauss = function(dat)
{
#========================================================================== 
# Get Method of Moments Estimates for ExGauss distribution
# cf. Heathcote (1996)
#==========================================================================
	M1    = mean(dat)
	M2    = var(dat) #division by n-1
	M3    = (length(dat)/(length(dat)-1)) * mean ( (dat-M1)^3 ) #effectively division by n-1

  if (M3 >= 0)
	 tau  = (0.5*M3)^(1/3)
	if (M3 < 0)
	 tau  = 0
  mu    = M1 - tau
	if ( (M2 - tau^2)>0 )
	{
		sigma = sqrt(M2 - tau^2) 
	}
	if ( ((M2 - tau^2)<=0) || (sigma<=0) || (tau<=0) ) #in case of weird values
	{
		tau   = 0.8 * sd(dat)
		mu    = M1 - tau
		#sigma = sqrt(M2 - tau^2) # this is the same as:
		sigma = 0.6 * sd(dat)	    # (correcting a small mistake in Heathcote 1996)
	}
	pars = c(mu, sigma, tau)
	return(pars)
}
# to check:
# mu=.4; sigma=.035; tau=.065
# dat = Gen.ExGauss(1000, mu, sigma, tau)
# MMest.ExGauss(dat)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Gen.ExGaussUnif = function(N, mu, sigma, tau, p_EG, min_U, max_U)
{
#========================================================================== 
# Generate data from ExGauss distribution Contaminated with a Uniform 
#==========================================================================
	EG = rnorm(round(p_EG*N), mean=mu, sd=sigma) +	rexp(round(p_EG*N), rate=1/tau)	
	Un = runif(N-round(p_EG*N), min=min_U, max=max_U)
	dat = c(EG,Un)
	
	return(dat)
}
# to check:
# mu=.4; sigma=.035; tau=.065; p_EG=.8; min_U=.1; max_U=1.2
# dat = Gen.ExGaussUnif(1000, mu, sigma, tau, p_EG, min_U, max_U)
# plot(dat)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Gen.ExGauss = function(N, mu, sigma, tau)
{
#========================================================================== 
# Generate data from ExGauss distribution
#==========================================================================
	dat = rnorm(N, mean=mu, sd=sigma) +	rexp(N, rate=1/tau)	
	return(dat)
}
# to check:
# mu=.4; sigma=.035; tau=.065
# dat = Gen.ExGauss(1000, mu, sigma, tau)
# mu+tau # theoretical mean
# mean(dat)
# sigma^2 + tau^2 # theoretical variance
# var(dat)
# 2*tau^3 # theoretical 3rd moment
# (length(dat)/(length(dat)-1)) * mean ( (dat-mean(dat))^3 )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
