#===============================
#subfunction
#===============================
#Poisson
#
# Sample size calculation
#
# Suppose that we have the following information.
# total number of genes for testing: m = 17306
# number of prognostic genes: m1 = 175
# number of true rejections: r1 = 140
# FDR level: f = 0.01
# ratio of total number of reads between two groups: w = 0.9
# power: 1-beta = 0.8
# the minimum average read counts among the prognostic genes in the control group: mu0 = 5

power_poisson<-function(n, w=1.0, rho=2.5, mu0=5.0, f, m, m1){

  estimate.power.given.n <- function(r1) {
    alpha_star<-(r1*f)/((m-m1)*(1-f))
    z_alpha<-qnorm(1-alpha_star/2)

    # Target.Func = abs(print(r1) - 1 + pnorm(z_alpha-(rho-1)*sqrt((n*mu0)/(rho/w + 1))))
    Target.Func = abs(r1 - 1 + pnorm(z_alpha-(rho-1)*sqrt((n*mu0)/(rho/w + 1))))
    return(Target.Func)
  }
  # power_w = tryCatch(uniroot(estimate.power.given.n, c(10^-6, 1))$root, error = function(e) 0)
  power_w = tryCatch(optimize(estimate.power.given.n, c(0, 1))$minimum, error = function(e) 0)
  ### Wald test
  # power_w<-1-pnorm( z_alpha-(rho-1)*sqrt( ( n*mu0 )/ ( rho/w + 1) ))

  if(0){
  ### score
  f2<-sqrt( (1 + rho*w )/(w+rho) )
  power_s<-1-pnorm( z_alpha*f2-(rho-1)*sqrt( ( n*mu0 )/ ( rho/w + 1) ))

  ### log transformation 1
  power_lw<-1-pnorm( z_alpha-log(rho)* sqrt( (n*mu0) / ( 1/(rho*w)+1 ) ) )

  ### log transformation  2
  f4<-sqrt( (2+w+1/w)/((1+w*rho)*(1+1/(w*rho))) )
  power_ls<-1-pnorm( z_alpha*f4-log(rho)* sqrt( (n*mu0) / ( 1/(rho*w)+1 ) ) )

  ### variance stabilizing
  power_tp<-1-pnorm( z_alpha*sqrt( (1+w)/(1+w*rho) )- 2*sqrt( n*mu0 +3/8)*( (sqrt(rho)-1)/sqrt(1/w+rho) ) )

  power_LRT<-function(n, d=0.5, rho=2, lambda=0.5, t=1, alpha=0.05, beta=0.2){
    p_value<-function( x0, x1, d=1){
      if (x1==0 & x0==0) {
        x1<-x1+0.5
        x0<-x0+0.5
      }
      if (x1 ==0 & x0!=0 ) x1<-x1+0.5
      if (x1 !=0 & x0==0 ) x0<-x0+0.5
      r1<-( (x1+x0) )/ (x1*(1+1/d))
      r0<-( (x1+x0) )/ (x0*(1+d))
      W<- -2*(x1*log(r1)+x0*log(r0) )
      1*(1-pchisq(W, 1) )
    }
    k1=0;k0=0
    mu0<-n*t*lambda
    mu1<-n*t*lambda*(rho*d)
    power<-0
    repeat{
      k1_incr<-0
      repeat{
        k1 <- k1 + 1
        power_incr <- dpois(k0, mu0 )*dpois(k1, mu1 )*( p_value(k0, k1, d) <= alpha )*1.0
        power <- power + power_incr
        k1_incr <- k1_incr + power_incr
        if ( ( power_incr < 10^(-30) ) &  k1>mu1 ) break
      }
      k1 <- 0
      k0 <- k0 + 1
      if ( ( k1_incr < 10^(-30) ) & k0>mu0 ) break
    }
    return(power)
  }
  power_lrt<-power_LRT(n=n, d=w, rho=rho, lambda=mu0, alpha=f)
  return( cbind(Wald=power_w, Score=power_s, LogTran1=power_lw, LogTran2=power_ls, VS=power_ls, LRT=power_lrt) )
  }
  return(Wald = power_w)
}


Poisson <- function(x, conds, target.n, w = 1, rho = 1.2, FDR = 0.05,
                             min.mu0 = 5, DE.prop = .1){

  x = x[apply(x, 1, mean) > 5, ]
  mean.group = apply(x, 1, function(y) tapply(y, conds, mean))
  fc = abs(log2((mean.group[1, ]+1)/(mean.group[2, ]+1)))
  # fc.prognostic = fc[order(fc, decreasing = TRUE)[1:round(nrow(x) * DE.prop)]]
  mu0 = max(min(round(mean.group[1, order(fc, decreasing = TRUE)[1:round(nrow(x) * DE.prop)]])), min.mu0)

  m<-nrow(x)
  m1<-DE.prop*m

  power=sapply(target.n,function(z) power_poisson(n=z, w=w, rho=rho, mu0=mu0, f=FDR, m=m, m1=m1))

  return(power)
}

#NB
NB.exact <- function(x, conds, target.n, rho = 1.4, DE.prop = 0.1,
                     min.mu0 = 5, valpha = 0.05, FDR = 0.05, w = 1){
  library("RnaSeqSampleSize")
  require("edgeR")
  x = x[apply(x, 1, mean) > 5, ]
  mean.group = apply(x, 1, function(y) tapply(y, conds, mean))
  fc = abs(log2((mean.group[1, ]+1)/(mean.group[2, ]+1)))
  # fc.prognostic=fc[order(fc,decreasing=TRUE)[1:round(nrow(x)*DE.prop)]]
  mu0 = max(min.mu0, min(round(mean.group[1, order(fc, decreasing = TRUE)[1:round(nrow(x)*DE.prop)]])))

  y <- DGEList(counts = x, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  y <- estimateTagwiseDisp(y) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  tagwise.dispersion = y$tagwise.dispersion

  w <- w
  m <- nrow(x)
  m1 <- DE.prop * m
  vphi_0 = max(tagwise.dispersion)

  power=sapply(target.n, function(z) est_power(n = z, alpha = valpha, w = w,
                                               f = FDR, m = m, m1 = m1, rho = rho,
                                               lambda0 = mu0, phi0 = vphi_0))
  return(power)
}

#RNASeqpower
RNASeqPower <- function(m, conds, target.n, n.prop = 1, effect, alpha = 0.05){
  library(RNASeqPower)
  library(edgeR)
  m=m[apply(m,1,mean)>5,]
  y <- DGEList(counts = m, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y, verbose=TRUE) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  y <- estimateTagwiseDisp(y) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  tagwise.dispersion = y$tagwise.dispersion
  ## y <- estimateTagwiseDisp(y, verbose=TRUE) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  ## plotBCV(y)
  BCV = sqrt(quantile(tagwise.dispersion, probs = 0.5))
  Power = sapply(target.n, function(x) rnapower(sum(m)/(ncol(m)*nrow(m)), n = x, n2 = x*n.prop, cv = BCV, effect = effect, alpha = alpha))
  #rnapower(depth, n, n2 = n, cv, cv2 = cv, effect, alpha, power)
  return(Power)
}



#PROPER
PROPER = function(x, conds, n.prop, target.n, DE.prop = 0.1, effect, default = T){
  library(edgeR)
  library(PROPER)

  x = x[apply(x, 1, mean) > 5, ] #filter mean count < 5
  mean.group = apply(x, 1, function(y) tapply(y, conds, mean))

  if(!default){
    fc = log((mean.group[1, ] +1)/(mean.group[2, ]+ 1))
    fc.abs = abs(fc)
    lfc = fc[order(fc.abs, decreasing = TRUE)[1:round(nrow(x) * DE.prop)]]
  }

  ngenes = nrow(x)

  y <- DGEList(counts = x, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateTagwiseDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags

  #log over-dispersion
  lOD = log(y$tagwise.dispersion)
  # lOD = log(1/50)

  #log baseline expression
  lBaselineExpr = log(rowMeans(x))

  #log fold change, since it's pilot data, we don't know, give true directly
  # lfc = log(2^rtruncnorm(ngenes, a=c(-Inf, 0.26), b=c(-0.26, Inf), mean = 0, sd = 0.2))

  # sim.opts = RNAseq.SimOptions.2grp(ngenes = 20000, DE.prop=0.1, lOD = "cheung", lBaselineExpr = "cheung")
  sim.opts = if(default){
    RNAseq.SimOptions.2grp(ngenes = ngenes, p.DE = DE.prop,
                           lOD = lOD, lBaselineExpr = lBaselineExpr)
  } else{
    RNAseq.SimOptions.2grp(ngenes = ngenes, p.DE = DE.prop, lfc = lfc,
                           # lOD = lOD, lBaselineExpr = lBaselineExpr, lfc = lfc)
                           lOD = lOD, lBaselineExpr = lBaselineExpr)
  }
  ptm <- proc.time()
  simres = runSims(Nreps = target.n, Nreps2 = target.n*n.prop, sim.opts = sim.opts, DEmethod = "edgeR", nsims = 20)
  ptm2 = proc.time()
  ptm2 - ptm

  powers = comparePower(simres, alpha.type = "fdr", alpha.nominal = 0.05,
                        stratify.by = "expr", target.by = "lfc",
                        strata = c(0, 5, 10, 20, 40, 80, 160, 320, 640, 1280, Inf),
                        delta = effect, filter.by = "expr", strata.filtered = 1) #filter < 5

  return(summaryPower(powers)[, "Marginal power"])
}

MixtureModel.Fittting.3parML <- function(p.value, init.r = NULL, init.s = s.lower, l.upper = 1, r.upper = .9, s.upper = Inf,
                                         l.lower = 0, r.lower = 0, s.lower = 1){
  library(limma)
  # lambda.Storey = qvalue(p.value)$pi0
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }


  fBUM <- function(z,d) {
    p.value = d
    lambda = z[1]
    r = z[2]
    s = z[3]
    # Target.Func = -sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    Target.Func = sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  if(is.null(init.r)) init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  # init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p

  #test <- function(x){
  #  init.lambda = sum(p.value>x)/(1-x)/length(p.value)
  #  return(init.lambda)
  #}
  init.lambda = c(convest(p.value)[[1]],runif(10, l.lower, l.upper))
  estimate <- function(x){
    print(x)
    init = c(x,max(0, min(init.r,.9)), init.s)
    if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B",
                                    upper = c(l.upper, r.upper, s.upper),
                                    lower=c(l.lower, r.lower, s.lower),
                                    control = list(maxit = 10000, fnscale = -1,
                                    factr = 10^-10))$par)), error = function(e) {print("error")})=="error")){return(init)}
    else{
      pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B",
                  upper = c(l.upper, r.upper, s.upper),
                  lower = c(l.lower, r.lower, s.lower),
                  control = list(maxit = 10000, fnscale = -1, factr = 10^-10))$par
      LL=sum(log(pars[1]+(1-pars[1])*p.value^(pars[2]-1)*(1-p.value)^(pars[3]-1)/beta(pars[2],pars[3])))
      return(c(pars,LL))
    }
  }
  tmp=lapply(init.lambda,estimate)
  tmp.ava=tmp[sapply(tmp,length)==4]
  ind=which.max(sapply(tmp.ava,function(x) x[4]))
  return(tmp.ava[[ind]])
}

#CDD
MixtureModel.Fittting.pi0 <- function(p.value, restrict = TRUE, l.upper = 0.9){
  library(limma)
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z,d) {
    p.value = d
    r = z[1]
    s = z[2]
    # Target.Func = -sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    Target.Func = sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }

  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p

  #test <- function(x){
  #  init.lambda = sum(p.value>x)/(1-x)/length(p.value)
  #  return(init.lambda)
  #}

  # lambda = min(convest(p.value)[[1]],0.9999)

  init = c(max(0,min(init.r,.9)),max(1,init.s))

  if(restrict){
    lambda = max(min(convest(p.value)[[1]],l.upper), 0.7)
    r.upper = .9
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  } else{
    lambda = min(convest(p.value)[[1]],0.9999)
    r.upper = 1
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  }
  if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B",
                                  upper=c(r.upper,s.upper),lower=c(r.lower,s.lower),
                                  control = list(maxit = 10000, fnscale = -1,
                                  factr = 10^-10))$par)),error = function(e) {print("error")})=="error")){return(init)}
  else{
    pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B",
                 upper=c(r.upper,s.upper),lower=c(r.lower,s.lower),
                 control = list(maxit = 10000, fnscale = -1, factr = 10^-10))$par
    LL=sum(log(lambda+(1-lambda)*p.value^(pars[1]-1)*(1-p.value)^(pars[2]-1)/beta(pars[1],pars[2])))
    return(c(lambda,pars,LL))
  }
}

MixtureModel.Fittting.pi0.2 <- function(p.value, restrict = TRUE){ #lambda set to .9
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z,d) {
    p.value = d
    r = z[1]
    s = z[2]
    Target.Func = -sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }

  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p; print(init.r)
  init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p; print(init.s)

  #test <- function(x){
  #  init.lambda = sum(p.value>x)/(1-x)/length(p.value)
  #  return(init.lambda)
  #}

  # lambda = min(convest(p.value)[[1]],0.9999)

  init = c(max(0,min(init.r, 0.9)),max(1,init.s))

  if(restrict){
    lambda = min(convest(p.value)[[1]],0.9)
    r.upper = .9
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  } else{
    lambda = min(convest(p.value)[[1]],0.9)
    r.upper = 1
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  }
  if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B",upper=c(r.upper,s.upper),lower=c(r.lower,s.lower))$par)),error = function(e) {print("error")})=="error")){return(init)}
  else{
    pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B",upper=c(r.upper,s.upper),lower=c(r.lower,s.lower))$par
    LL=sum(log(lambda+(1-lambda)*p.value^(pars[1]-1)*(1-p.value)^(pars[2]-1)/beta(pars[1],pars[2])))
    return(c(lambda,pars,LL))
  }
}

MixtureModel.Fittting.pi0.3 <- function(p.value, restrict = TRUE, l.upper = 0.9){
  library(limma)
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z,d) {
    p.value = d
    r = z[1]
    s = z[2]
    # Target.Func = -sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    Target.Func = -sum(log(lambda+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }

  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p

  #test <- function(x){
  #  init.lambda = sum(p.value>x)/(1-x)/length(p.value)
  #  return(init.lambda)
  #}

  # lambda = min(convest(p.value)[[1]],0.9999)

  init = c(max(0,min(init.r,.9)),max(1,init.s))

  if(restrict){
    lambda = max(min(convest(p.value)[[1]],l.upper), 0.7)
    r.upper = .9
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  } else{
    lambda = min(convest(p.value)[[1]],0.9999)
    r.upper = 1
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  }
  if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B",
                                   upper=c(r.upper,s.upper),lower=c(r.lower,s.lower))$par)),error = function(e) {print("error")})=="error")){return(init)}
  else{
    pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B",
                 upper=c(r.upper,s.upper),lower=c(r.lower,s.lower))$par
    LL=sum(log(lambda+(1-lambda)*p.value^(pars[1]-1)*(1-p.value)^(pars[2]-1)/beta(pars[1],pars[2])))
    return(c(lambda,pars,LL))
  }
}

#two-beta

MixtureModel.two.beta <- function(p.value, init.r = NULL, init.s = s.lower, l.upper = .9, r.upper = .9, s.upper = Inf, r2.upper = Inf, s2.upper = 1,
                                  l.lower = .7, r.lower = 0, s.lower = 1, r2.lower = 1, s2.lower = 0){
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z,d) {
    p.value = d
    lambda = z[1]
    r = z[2]
    s = z[3]
    r2 = z[4]
    s2 = z[5]

    # Target.Func = -sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    Target.Func = sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  if(is.null(init.r)) init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  # if(is.null(init.s)) init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p


  # lambda = min(convest(p.value)[[1]],0.9999)
  lambda = convest(p.value)[[1]]
  # estimate <- function(x){

  try.init.s = lapply(c(1, seq(10, 100, 10)), function(s){
    init = c(lambda, max(0, min(init.r,.9)), s, 1.5, 1)
    # init = c(x, max(0,min(init.r,.9)),max(1,init.s), 1.5, 1)
    # l.upper = .9; r.upper = .9; s.upper = Inf; r2.upper = Inf; s2.upper = 1
    # l.lower = .7; r.lower = 0; s.lower = 1; r2.lower = 1; s2.lower = 0
    # print(init)
    # if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B", upper=c(l.upper, r.upper, s.upper, r2.upper, s2.upper),lower=c(l.lower, r.lower, s.lower, r2.lower, s2.lower))$par)),error = function(e) {print("error")})=="error")){return(init)}
    # else{
    #   pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B", upper=c(l.upper, r.upper, s.upper, r2.upper, s2.upper),lower=c(l.lower, r.lower, s.lower, r2.lower, s2.lower))$par
    #   LL=sum(log(pars[1]*p.value^(pars[4]-1)*(1-p.value)^(pars[5]-1)/beta(pars[4],pars[5])+(1-pars[1])*p.value^(pars[2]-1)*(1-p.value)^(pars[3]-1)/beta(pars[2],pars[3])))
    #   return(c(pars,LL))
    # }
    model = tryCatch(optim(init, fBUM, d = p.value, method = "L-BFGS-B",
                      upper = c(l.upper, r.upper, s.upper, r2.upper, s2.upper),
                      lower = c(l.lower, r.lower, s.lower, r2.lower, s2.lower),
                      control = list(maxit = 10000, fnscale = -1,
                                     factr = 10^-10)),
                      error = function(e) "error")
    return(model)
  })

  i = which.max(sapply(try.init.s, function(x){
    if(all(x == "error")){
      NA
    } else {
      LL = x$value
  }}))

  # if(all(model == "error")){
  #   pars = init
  # } else {
  #   pars = model$par
  #   LL = model$value
  #   return(c(pars, LL))
  # }
  pars = try.init.s[[i]]$par
  LL = try.init.s[[i]]$value
  return(c(pars, LL))
  # }
  # init.lambda = c(convest(p.value)[[1]],runif(10,0.5,1))
  # tmp=lapply(init.lambda,estimate)
  # tmp.ava=tmp[sapply(tmp,length)==4]
  # ind=which.max(sapply(tmp.ava,function(x) x[4]))
  # return(tmp.ava[[ind]])
}

#two-beta + uniform
MixtureModel.3.beta <- function(p.value, init.r = NULL, init.s = s.lower, l.upper = .9, l2.upper = .9, r.upper = .9, s.upper = Inf, r2.upper = Inf, s2.upper = 1,
                                  l.lower = 0.7, l2.lower = 0, r.lower = 0, s.lower = 4, r2.lower = 1, s2.lower = 0){
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z,d) {
    p.value = d
    lambda = z[1]
    r = z[2]
    s = z[3]
    r2 = z[4]
    s2 = z[5]
    lambda2 = z[6]
    a = log(lambda2*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(lambda - lambda2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s))
    if(sum(!is.finite(a))) print(z)
    Target.Func = -sum(log(lambda2*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(lambda-lambda2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  if(is.null(init.r)) init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  # if(is.null(init.s)) init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p

  # lambda = min(convest(p.value)[[1]],0.9999)
  lambda = convest(p.value)[[1]]
  init = c(lambda, max(0,min(init.r,.9)),max(1,init.s), 1.5, 1, 0.1)
  #l.upper = .9; r.upper = .9; s.upper = Inf; r2.upper = Inf; s2.upper = 1
  #l.lower = .7; r.lower = 0; s.lower = 1; r2.lower = 1; s2.lower = 0
  print(init)
  if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B", upper=c(l.upper, r.upper, s.upper, r2.upper, s2.upper, l2.upper),lower=c(l.lower, r.lower, s.lower, r2.lower, s2.lower, l2.lower))$par)),error = function(e) {print("error")})=="error")){return(init)}
  else{
    pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B", upper=c(l.upper, r.upper, s.upper, r2.upper, s2.upper, l2.upper),lower=c(l.lower, r.lower, s.lower, r2.lower, s2.lower, l2.lower))$par
    LL=sum(log(pars[6]*p.value^(pars[4]-1)*(1-p.value)^(pars[5]-1)/beta(pars[4],pars[5])+pars[1]-pars[6]+(1-pars[1])*p.value^(pars[2]-1)*(1-p.value)^(pars[3]-1)/beta(pars[2],pars[3])))
    return(c(pars,LL))
  }
}

#fix lambda
MixtureModel.two.beta.2 <- function(p.value, lambda = NULL, init.r = NULL, init.s = s.lower, l.upper = 0.9, r.upper = .9, s.upper = Inf, r2.upper = Inf, s2.upper = 1,
                                    r.lower = 0, s.lower = 1, r2.lower = 1, s2.lower = 0){
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z, d, lambda) {
    p.value = d
    r = z[1]
    s = z[2]
    r2 = z[3]
    s2 = z[4]

    # Target.Func = -sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    Target.Func = sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  # init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  # init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p
  if(is.null(init.r)) init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  # if(is.null(init.s)) init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p

  # lambda = min(convest(p.value)[[1]],0.9999)
  if(is.null(lambda)) lambda = max(min(convest(p.value)[[1]], l.upper), 0.7)

  try.init.s = lapply(c(1, seq(10, 100, 10)), function(s){

    init = c(max(0,min(init.r,.9)), s, 1.5, .9)
    #l.upper = .9; r.upper = .9; s.upper = Inf; r2.upper = Inf; s2.upper = 1
    #l.lower = .7; r.lower = 0; s.lower = 1; r2.lower = 1; s2.lower = 0


    model = tryCatch(optim(init, fBUM, d = p.value, lambda = lambda, method = "L-BFGS-B",
                           upper = c(r.upper, s.upper, r2.upper, s2.upper),
                           lower = c(r.lower, s.lower, r2.lower, s2.lower),
                           control = list(maxit = 10000, fnscale = -1,
                                          factr = 10^-10)),
                     error = function(e) "error")
    return(model)
  })

  i = which.max(sapply(try.init.s, function(x){
    if(all(x == "error")){
      NA
    } else {
      LL = x$value
  }}))

  pars = try.init.s[[i]]$par
  LL = try.init.s[[i]]$value

  # if(all(model == "error")){
  #   pars = init
  #   return(c(lambda, pars))
  # } else {
  #   pars = model$par
  #   LL = model$value
  #   return(c(lambda, pars, LL))
  # }
  return(c(lambda, pars, LL))
}

MixtureModel.two.beta.3 <- function(p.value, l.upper = .9, r.upper = .9, s.upper = Inf, r2.upper = Inf,
                                  l.lower = .7, r.lower = 0, s.lower = 1, r2.lower = 1){
  if (all(is.na(p.value))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.value
  if (any(is.na(p.value))) {
    p.value <- p.value[!is.na(pvals)]
  }
  if (min(p.value) == 0) {
    min.nonzero <- min(p.value[p.value > 0])
    p.value[p.value == 0] <- min.nonzero/2
  }

  fBUM <- function(z,d) {
    p.value = d
    lambda = z[1]
    r = z[2]
    s = z[3]
    r2 = z[4]
    s2 = 1
    Target.Func = -sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  ## inital value(MME and non-parametric)
  mean.p = mean(p.value)
  var.p = var(p.value)
  init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p

  # lambda = min(convest(p.value)[[1]],0.9999)
  lambda = convest(p.value)[[1]]
  init = c(lambda, max(0,min(init.r,.9)),max(1,init.s), 1.5)

  #l.upper = .9; r.upper = .9; s.upper = Inf; r2.upper = Inf; s2.upper = 1
  #l.lower = .7; r.lower = 0; s.lower = 1; r2.lower = 1; s2.lower = 0

  if(unique(tryCatch(print(c(optim(init, fBUM,d=p.value,method= "L-BFGS-B", upper=c(l.upper, r.upper, s.upper, r2.upper),lower=c(l.lower, r.lower, s.lower, r2.lower))$par)),error = function(e) {print("error")})=="error")){return(init)}
  else{
    pars = optim(init, fBUM,d=p.value,method= "L-BFGS-B", upper=c(l.upper, r.upper, s.upper, r2.upper),lower=c(l.lower, r.lower, s.lower, r2.lower))$par
    s2 = 1
    pars = c(pars, s2)
    LL=sum(log(pars[1]*p.value^(pars[4]-1)*(1-p.value)^(pars[5]-1)/beta(pars[4],pars[5])+(1-pars[1])*p.value^(pars[2]-1)*(1-p.value)^(pars[3]-1)/beta(pars[2],pars[3])))
    return(c(pars,LL))
  }
}

Estimate.EDR.from.pilot <- function(Data, status, group.name = c("Control", "Case"), FDR, M,
                                    target.N, target.R = NULL, target.theta = NULL,
                                    method,  transform.null = F, tol = 0.1,
                                    filter = T, filter.level = 5, resample.DE = F,
                                    fc.cut = 1, p.cut = 10^-4, s.lower = 1,
                                    l.upper = 0.9, know.DE = F, consider.lfc = F){

  N0 = sum(status == group.name[1])
  N1 = sum(status == group.name[2])
  theta = N1/N0
  if(is.null(target.theta)) target.theta = theta

  mean.gene = apply(Data, 1, mean)

  Data.filter = if(filter) Data[which(mean.gene>filter.level),] else Data

  R = mean(apply(Data.filter, 2, sum))
  if(is.null(target.R)) target.R = R

  # ngenes= nrow(Data.filter)

  Pilot.data = Data.filter
  colnames(Pilot.data)=1:ncol(Data.filter)

  library(edgeR)
  # library(qvalue)
  y <- DGEList(counts = Pilot.data, group = status)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  # et <- exactTest(y,dispersion="common")$table$P #Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
  delta = 1/y$common.dispersion
  ##56.37761

  GLM.fit.each <- function(z){
    y0 = z[which(status==group.name[1])]
    y1 = z[which(status==group.name[2])]
    f <- function(z,d) {
      delta = d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1 = d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(delta+y0)-lgamma(delta)-lgamma(y0+1)+y0*log(R/delta*exp(beta.0))-(y0+delta)*log(1+R/delta*exp(beta.0)))-sum(lgamma(delta+y1)-lgamma(delta)-lgamma(y1+1)+y1*log(R/delta*exp(beta.0+beta.1))-(y1+delta)*log(1+R/delta*exp(beta.0+beta.1)))
      return(Target.Func)
    }

    pt1 <- optim(c(0, 0), f, d = list(delta, R, y0, y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[2])
    var.beta.1 = (1/N0)*((1 + theta*r0)/(theta*R*r1) + (1+theta)/(theta*delta))
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }


  OUT.pars = t(apply(round(Pilot.data),1,GLM.fit.each))
  colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")

  model = cbind(OUT.pars[,1:2], delta)

  p.value = 2*pnorm(-abs(OUT.pars[,3]))
  mean.count = rowMeans(Pilot.data)

  if(0){
    pdf("mean count vs p.pdf")
    plot(p.value, mean.gene, col = "#00000060", cex = .5, pch = 16, ylim = c(0, 100))
    abline(h = 10, col = 2)
    dev.off()
  }
  q.value = p.adjust(p.value, method = "BH")

  mean.by.group = apply(Pilot.data,1,function(x) tapply(x, status, mean))
  fold.change = (mean.by.group[1, ] + 1)/(mean.by.group[2, ] + 1)
  lfc = abs(log2(fold.change))
  large.lfc = lfc > log2(fc.cut)

  if(0){
  #p.value.2 = p.value
  if(resample.DE){ #p-value calibration
    resampling.DE = which(lfc < log2(fc.cut) & p.value < p.cut)
    if(length(resampling.DE)){
      p.value.unif = runif(length(resampling.DE), 0, 1)
      p.value[resampling.DE] = p.value.unif
    }
  }
  }
  is_DE = if(know.DE) sapply(strsplit(names(p.value), split = "_"), function(x) x[3] == "T") else NULL

  true.lambda = if(length(is_DE)) mean(!is_DE) else NULL

  Fitted.Model = switch(method, "3par" = MixtureModel.Fittting.3parML(p.value),
                        "CDD" = MixtureModel.Fittting.pi0(p.value),
                        "CDD.95" = MixtureModel.Fittting.pi0(p.value, l.upper = 0.95),
                        "CDDnoR" = MixtureModel.Fittting.pi0(p.value, restrict = F),
                        "TwoBeta" = MixtureModel.two.beta(p.value, s.lower = s.lower, l.upper = l.upper), #s >= 4
                        "TwoBetaNew" = MixtureModel.two.beta.3(p.value), #fix s2 = 1
                        # "TwoBetaCDD" = MixtureModel.two.beta.2(p.value, s2.upper = .5),
                        "TwoBetaCDD" = MixtureModel.two.beta.2(p.value, s.lower = s.lower, l.upper = l.upper), #use CDD estimate
                        "TwoBetaCDD.95" = MixtureModel.two.beta.2(p.value, l.upper = .95),
                        "TwoBetaLambda1" = MixtureModel.two.beta.2(p.value, l.upper = 1),
                        "TwoBetaLambda.95" = MixtureModel.two.beta.2(p.value, l.upper = .95))

  check.if.use.CDD = function(model, uniform.cut = .5){
    r = model[4]
    s = model[5]
    pbeta(1, r, s) - pbeta(uniform.cut, r, s) - uniform.cut
  }

  index.use.CDD = check.if.use.CDD(Fitted.Model)

  if(method == "TwoBeta" & index.use.CDD < tol) {
    print("Use CDD to estimate lambda")
    Fitted.Model = MixtureModel.two.beta.2(p.value, s.lower = s.lower, l.upper = l.upper)
  }
  # Fitted.Model = MixtureModel.Fittting.3parML(p.value)
  # Fitted.Model = MixtureModel.Fittting.pi0(p.value)
  # Fitted.Model = MixtureModel.Fittting.pi0(p.value, restrict = F)
  # Fitted.Model = MixtureModel.Fittting.pi0.2(p.value, restrict = F) #only restrict lambda as 0.9
  # Fitted.Model = MixtureModel.two.beta(p.value); Fitted.Model
  # Fitted.Model = MixtureModel.two.beta(p.value, l.upper = 1)
  p.value.mod <- p.value
  if (any(is.na(p.value.mod))) {
    p.value.mod <- p.value.mod[!is.na(p.value.mod)]
    model = model[!is.na(p.value.mod),]
    large.lfc = large.lfc[!is.na(p.value.mod)]
  }
  if (min(p.value.mod) == 0) {
    min.nonzero <- min(p.value.mod[p.value.mod > 0])
    p.value.mod[p.value.mod == 0] <- min.nonzero/2
  }
  if (max(p.value.mod) == 1) {
    max.nonone <- max(p.value.mod[p.value.mod < 1])
    p.value.mod[p.value.mod == 1] <- 1-(1-max.nonone)/2
  }

  ngenes = length(p.value.mod)

  if(Fitted.Model=="error"){
    return("error")
  }
  else{
    lambda = as.numeric(Fitted.Model[1])
    r=as.numeric(Fitted.Model[2])
    s=as.numeric(Fitted.Model[3])
    if(method == "TwoBeta"){
      r2=as.numeric(Fitted.Model[4])
      s2=as.numeric(Fitted.Model[5])
      posterior = (lambda*dbeta(p.value.mod,r2,s2))/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda*dbeta(p.value.mod,r2,s2)) #prob be non-DE
    } else posterior = lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda) #prob be non-DE

    ## posterior FDR
    # sample.size=c(n:12,20,30,40,50,100)
    # posterior.two.beta = (lambda*dbeta(p.value.mod,r2,s2))/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda*dbeta(p.value.mod,r2,s2)) #prob be non-DE
    # posterior.CDD=lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda) #prob be non-DE
    # posterior.CDD.unrestrict=lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda) #prob be non-DE
    # posterior.CDD.unrestrict.new=lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda) #prob be non-DE
    #
    # posterior = posterior.two.beta
    # posterior = posterior.CDD
    # posterior = posterior.CDD.unrestrict
    # posterior = posterior.CDD.unrestrict.new

    #check posterior prob
    if(0){
      #p-value vs. posterior
      # pdf("p-value vs posterior CDD vs. 2 beta mixture new.pdf")
      pdf("p-value.pdf")

      p.hist.compare.posterior.prob = function(p.value.mod, posterior){
        hist(p.value.mod, xlab = "p-value", main = "")
        par(new = T)
        plot(p.value.mod, 1-posterior.CDD, cex = .5, pch = 15, col = 2, axes = F, ylab = "", xlab = ""); axis(4)
        points(p.value.mod, 1-posterior.two.beta, cex = .5, pch = 15, col = 4)
        # points(p.value.mod, 1-posterior.CDD.unrestrict, cex = .5, pch = 15, col = 3)
        points(p.value.mod, 1-posterior.CDD.unrestrict.new, cex = .5, pch = 15, col = 3)
        # legend("top", c("Two beta mixture", "CDD restricted", "CDD unrestricted", "CDD unrestricted2"), col = c(4, 2, 3, 5), pch = 16, title = "DE posterior prob.")
        legend("top", c("Two beta mixture", "CDD restricted", "CDD unrestricted"), col = c(4, 2, 3), pch = 16, title = "DE posterior prob.")
      }
      dev.off()
    }
    if(lambda>=0.99){
      return(Fitted.Model)
    }
    else{
      parameter = list(n.old = N0, n.new = target.N,
                       R.old = R, R.new = target.R,
                       theta.old = theta, theta.new = target.theta)
      Resampling <- function(target.N, target.R){
        # DE_status_bootstrap = sample(c(FALSE,TRUE),ngenes,prob = c(lambda,1-lambda),replace=TRUE)
        DE_status_posterior = sapply(1:length(posterior),function(x) sample(c(FALSE,TRUE),1,prob = c(posterior[x],1-posterior[x]),replace=TRUE))
        #####################################################
        ## Posterior
        #########################################
        # transform.p.value.old <- function(p.value.each, DE_status_each, n.old, n.new, transform.null = F){
        #   if(DE_status_each){
        #     statistc.old = qnorm(p.value.each/2,lower.tail=FALSE)
        #     statistc.new = statistc.old*sqrt(n.new)/sqrt(n.old)
        #     p.value.star.star = (1-pnorm(abs(statistc.new)))*2
        #     return(p.value.star.star)
        #   } else if(transform.null) return(rep(runif(1),length(n.new)))
        #   else return(rep(p.value.each, length(n.new)))
        #   # } else{return(rep(p.value.each,length(n.new)))}
        # }

        transform.p.value <- function(p.value.each, DE_status_each, parameter, model, transform.null = F){
          n.old = parameter[[1]]
          n.new = parameter[[2]]
          R.old = parameter[[3]]
          R.new = parameter[[4]]
          theta.old = parameter[[5]]
          theta.new = parameter[[6]]

          adjust = function(n.old, n.new, R.old, R.new, theta.old, theta.new){
            adjust.pilot = (1/n.old)*((1+theta.old*exp(beta.estimate[2]))/(theta.old*R.old*exp(beta.estimate[1]+beta.estimate[2]))+(1+theta.old)/(theta.old*delta))

            adjust.matrix = outer(R.new, n.new, function(x, y){
              adjust.target = (1/y)*((1+theta.new*exp(beta.estimate[2]))/(theta.new*x*exp(beta.estimate[1]+beta.estimate[2]))+(1+theta.new)/(theta.new*delta))
            })
            return(sqrt(adjust.pilot/adjust.matrix))
          }
          if(DE_status_each){
            beta.estimate = model[1:2]
            delta = model[3]
            statistic.old = qnorm(p.value.each/2,lower.tail=FALSE)
            # statistc.new = statistic.old*sqrt(n.new)/sqrt(n.old)
            # adjust.old = (1/n.new)*((1+theta.new*exp(beta.estimate[2]))/(theta.new*R.new*exp(beta.estimate[1]+beta.estimate[2]))+(1+theta.new)/(theta.new*delta))
            # adjust.new = (1/n.old)*((1+theta.old*exp(beta.estimate[2]))/(theta.old*R.old*exp(beta.estimate[1]+beta.estimate[2]))+(1+theta.old)/(theta.old*delta))
            statistic.adjust = adjust(n.old, n.new, R.old, R.new, theta.old, theta.new)
            # statistc.new = statistic.old*sqrt(adjust.new/adjust.old)
            statistic.new = statistic.old*statistic.adjust

            p.matrix = (1-pnorm(abs(statistic.new)))*2
            colnames(p.matrix) = n.new
            rownames(p.matrix) = R.new
            return(p.matrix)
          # } else if(transform.null) return(rep(runif(1),length(n.new)))
          } else if(transform.null){
            p.matrix = matrix(rep(runif(1),length(n.new)*length(R.new)), nrow = length(R.new))
            colnames(p.matrix) = n.new
            rownames(p.matrix) = R.new
            return(p.matrix)
          }
          # else return(rep(p.value.each, length(n.new)))
          else {
            p.matrix = matrix(rep(p.value.each, length(n.new)*length(R.new)), nrow = length(R.new))
            colnames(p.matrix) = n.new
            rownames(p.matrix) = R.new
            return(p.matrix)
          }
          # } else{return(rep(p.value.each,length(n.new)))}
        }

        # p.value.star.posterior.old = sapply(1:ngenes,function(x) transform.p.value.old(p.value.mod[x], DE_status_posterior[x], N0,
        #                                                                            sample.size, transform.null = transform.null))

        p.value.updated = lapply(1:ngenes,function(x) transform.p.value(p.value.mod[x], DE_status_posterior[x], parameter = parameter,
                                                                               model = model[x,], transform.null = transform.null))

        result = lapply(1:length(target.R), function(i){ #for each depth
          p.value.star.posterior = t(sapply(p.value.updated, function(x) x[i,]))
          #for each gene, take out the i-th depth, the length is as n.new
          #target.N x G
          # p.value.star.posterior = matrix(p.transform.by.n, ncol = length(target.N), byrow = T)
          # p.value.star.posterior = p.transform.by.n

          Estimate_Posterior <- function(p.value.star.star){
            if (min(p.value.star.star) == 0) {
              min.nonzero <- min(p.value.star.star[p.value.star.star > 0])
              p.value.star.star[p.value.star.star == 0] <- min.nonzero/2
            }
            if (max(p.value.star.star) == 1) {
              max.non1 <- max(p.value.star.star[p.value.star.star <1])
              p.value.star.star[p.value.star.star == 1] <- (max.non1+1)/2
            }
            ## empirical FDR control
            p.DE=p.value.star.star[DE_status_posterior]
            p.nonDE=p.value.star.star[!DE_status_posterior]
            p.unique=sort(unique(p.value.star.star))

            FDR.unique <- vector(length=length(p.unique))
            for(i in 1:length(p.unique)){
              FDR.unique[i]=sum(p.nonDE<=p.unique[i])/sum(p.value.star.star<=p.unique[i])
            }
            index = which(FDR.unique<=FDR)
            p.value.cut = if(length(index)){
              p.unique[max(index)]
            } else NA #everything is beyond FDR 0.05
            if(min(p.value.star.star)>p.value.cut | is.na(p.value.cut)){Declare_status=rep("nonDE",ngenes)}
            else{
              Declare_status = rep("nonDE", length(DE_status_posterior))
              # if(consider.lfc) Declare_status[which(p.value.star.star <= p.value.cut & large.lfc)] = "DE"
              # else Declare_status[which(p.value.star.star <= p.value.cut)] = "DE"
              Declare_status[which(p.value.star.star <= p.value.cut)] = "DE"
              # Declare_status[-which(p.value.star.star<=p.value.cut)] = "nonDE"
            }
            A = sum((Declare_status=="nonDE")*(!DE_status_posterior))
            # B = sum((Declare_status=="nonDE")*(DE_status_posterior))

            if(consider.lfc){
              B = sum((Declare_status=="nonDE")*(DE_status_posterior)*large.lfc)
              D = sum((Declare_status=="DE")*(DE_status_posterior)*large.lfc)
              D/(B+D)
            } else {
              B = sum((Declare_status=="nonDE")*(DE_status_posterior))
              D = sum((Declare_status=="DE")*(DE_status_posterior))
              D/(B+D)
            }
            C = sum((Declare_status=="DE")*(!DE_status_posterior))
            if((C+D)==0){
              TP_hat_post=0
              ## no declared genes
            }
            else{
              TP_hat_post = D/(C+D)
            }
            if((A+B)==0){
              TN_hat_post=0
              ## no declared genes
            }
            else{
              TN_hat_post = A/(A+B)
            }
            # if(mean(p.value)>0.5+1.96*sqrt(1/(12*length(p.value)))){EDR_post=D/length(p.value)}
            # else{EDR_post = D/(B+D)}
            EDR_post = D/(B+D)
            Declare_post  = sum(Declare_status=="DE")
            # return(c(TP_hat_post,TN_hat_post,EDR_post,Declare_post, B = B, D = D, p.value.cut = p.value.cut))
            return(c(TP = TP_hat_post, TN = TN_hat_post, EDR = EDR_post, DE = Declare_post, B = B, D = D,
                     p.value.cut = p.value.cut, FDR = FDR.unique[max(index)],
                     null.num = sum(p.nonDE<=p.unique[max(index)]), all.num = sum(p.value.star.star<=p.unique[max(index)])))
          }

          Estimate.Posterior.Result = matrix(apply(p.value.star.posterior,2,Estimate_Posterior), ncol = length(target.N)); round(Estimate.Posterior.Result, 3)

          #debug, check p-value distribution after transformation
          if(0){
            temp = cbind(p.value.mod, p.value.star.posterior)

            # pdf("check p-value transformation two beta mixture.pdf", width = 14, height = 14)
            # pdf("check p-value transformation CDD restricted.pdf", width = 14, height = 14)
            # pdf("check p-value transformation CDD restricted.pdf", width = 14, height = 14)
            pdf("check p-value transformation CDD unrestricted new initial r.pdf", width = 14, height = 14)
            par(mfrow = c(3, 3))
            for(i in 1:ncol(temp)){
              if(i == 1) hist(temp[,i], main = "Pilot N = 2", xlab = "p-value", ylim = c(0, 1300), breaks = 20)
              else{
                hist(temp[,i], main = paste("Predicted N = ", target.N[i-1], sep = ""), xlab = "p-value", ylim = c(0, 1300), breaks = 20)
                abline(v = Estimate.Posterior.Result[7,i-1], col = 2, lty = 1, lwd = 2)
              }
            }
            dev.off()
          }

          row.names(Estimate.Posterior.Result) = c("TP", "TN", "EDR", "DeclareDE", "B", "D", "FDRcut", "FDR", "null.num", "all.num")
          colnames(Estimate.Posterior.Result) = target.N
          return(Estimate.Posterior.Result)
        })
        return(result)
      }
      Result = lapply(1:M, function(x){
        # print(x)
        Resampling(target.N = parameter[[2]], target.R = parameter[[4]])
      })
      return(list(Result = Result, Model = Fitted.Model, p = p.value,
                  fold.change = fold.change,
                  Number.DE = sum(q.value < FDR), mean.count = mean.count,
                  index.use.CDD = index.use.CDD, #num.p.calibrate = length(resampling.DE),
                  true.lambda = true.lambda, is_DE = is_DE))
    }
  }
}

Estimate.true.EDR = function(Data, status, group.name = c("Control", "Case"),
                             target.theta = NULL, FDR, filter = T,
                             filter.level = 5, resample.DE = F,
                             fc.cut = 1, p.cut = 10^-4, consider.lfc = F){
  N0 = sum(status == group.name[1])
  N1 = sum(status == group.name[2])
  theta = N1/N0

  mean.gene = apply(Data, 1, mean)

  ngenes = nrow(Data)

  Data.filter = if(filter) Data[which(mean.gene > filter.level),] else Data

  mean.gene.filter = apply(Data.filter, 1, mean)

  ngenes.after.filter= nrow(Data.filter)

  is_DE = sapply(strsplit(row.names(Data.filter), split = "_"), function(x) x[3] == "T")

  num.DE = sum(is_DE)

  True.lambda = mean(!is_DE)

  R = mean(apply(Data.filter, 2, sum))

  Pilot.data = Data.filter
  colnames(Pilot.data)=1:ncol(Data.filter)

  mean.by.group = apply(Pilot.data,1,function(x) tapply(x, status, mean))
  fold.change = (mean.by.group[1, ] + 1)/(mean.by.group[2, ] + 1)

  library(edgeR)
  y <- DGEList(counts = Pilot.data, group = status)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  # et <- exactTest(y,dispersion="common")$table$P #Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
  delta = 1/y$common.dispersion

  GLM.fit.each <- function(z){
    y0 = z[which(status==group.name[1])]
    y1 = z[which(status==group.name[2])]
    f <- function(z,d) {
      delta = d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1 = d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(delta+y0)-lgamma(delta)-lgamma(y0+1)+y0*log(R/delta*exp(beta.0))-(y0+delta)*log(1+R/delta*exp(beta.0)))-sum(lgamma(delta+y1)-lgamma(delta)-lgamma(y1+1)+y1*log(R/delta*exp(beta.0+beta.1))-(y1+delta)*log(1+R/delta*exp(beta.0+beta.1)))
      return(Target.Func)
    }

    pt1 <- optim(c(0, 0), f, d = list(delta, R, y0, y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[2])
    var.beta.1 = (1/N0)*((1 + theta*r0)/(theta*R*r1) + (1+theta)/(theta*delta))
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }

  OUT.pars = t(apply(round(Pilot.data), 1, GLM.fit.each))
  colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")

  p.value = 2*pnorm(-abs(OUT.pars[,3]))

  lfc = abs(log2(fold.change))
  large.lfc = lfc > log2(fc.cut)

  if(0){
  #p.value.2 = p.value
  if(resample.DE){ #p-value calibration
    resampling.DE = which(lfc < log2(fc.cut) & is_DE)
    # resampling.DE = which(lfc < log2(fc.cut) & p.value < p.cut)
    if(length(resampling.DE)){
      p.value.unif = runif(length(resampling.DE), 0, 1)
      p.value[resampling.DE] = p.value.unif
      is_DE[resampling.DE] = F
    }
  }
  }
  True.lambda.post = mean(!is_DE)


  ## empirical FDR control
  p.DE=p.value[is_DE]
  p.nonDE=p.value[!is_DE]
  p.unique=sort(unique(p.value))

  FDR.unique <- vector(length=length(p.unique))
  for(i in 1:length(p.unique)){
    FDR.unique[i]=sum(p.nonDE<=p.unique[i])/sum(p.value<=p.unique[i])
  }

  if(min(FDR.unique)>FDR){return("Cannot control FDR")}
  else{
    index = which(FDR.unique<=FDR)
    p.value.cut = if(length(index)){
      p.unique[max(index)]
    } else NA #everything is beyond FDR 0.05
    if(min(p.value)>p.value.cut | is.na(p.value.cut)){Declare_status=rep("nonDE",ngenes)}
    else{
      Declare_status = rep("nonDE", length(is_DE))
      # if(consider.lfc) Declare_status[which(p.value <= p.value.cut & large.lfc)] = "DE"
      # else Declare_status[which(p.value <= p.value.cut)] = "DE"
      Declare_status[which(p.value <= p.value.cut)] = "DE"
      # Declare_status[-which(p.value<=p.value.cut)] = "nonDE"
    }
    A = sum((Declare_status=="nonDE")*(!is_DE))
    # B = sum((Declare_status=="nonDE")*(is_DE))
    if(consider.lfc){
      B = sum((Declare_status=="nonDE")*(is_DE)*(large.lfc))
      D = sum((Declare_status=="DE")*(is_DE)*(large.lfc))
    } else{
      B = sum((Declare_status=="nonDE")*(is_DE))
      D = sum((Declare_status=="DE")*(is_DE))
    }
    C = sum((Declare_status=="DE")*(!is_DE))
    # D = sum((Declare_status=="DE")*(is_DE))

    if((C+D)==0){
      TP_hat_true=0
      ## no declared genes
    }
    else{
      TP_hat_true = D/(C+D)
    }
    if((A+B)==0){
      TN_hat_true=0
      ## no declared genes
    }
    else{
      TN_hat_true = A/(A+B)
    }
    EDR_true = D/(B+D)
    Declare_true  = sum(Declare_status=="DE")

    return(list(p.value = p.value, lfc = lfc, is_DE = is_DE,
          Result = c(TP_hat_true = TP_hat_true, TN_hat_true = TN_hat_true,
                     EDR_true = EDR_true, Declare_true = Declare_true, B = B, D = D,
              p.value.cut = p.value.cut, FDR = FDR.unique[max(index)],
              null.num = sum(p.nonDE <= p.unique[max(index)]),
              all.num = sum(p.value <= p.unique[max(index)]),
              ngenes.after.filter = ngenes.after.filter, True.lambda = True.lambda,
              True.lambda.post = True.lambda.post)))
  }

}


DE <- function(Data, status, group.name = c("Control", "Case"),
               method = "GLM", FDR, filter.by.5 = T,
               resample.DE = T, fc.cut = 1.4, p.cut = 10^-4){
  N0 = sum(status == group.name[1])
  N1 = sum(status == group.name[2])
  theta = N1/N0

  mean.gene=apply(Data,1,mean)

  Data.filter = if(filter.by.5) Data[mean.gene > 5,] else Data

  R=mean(apply(Data.filter,2,sum))

  Pilot.data = Data.filter
  colnames(Pilot.data)=1:ncol(Data.filter)

  library(edgeR)
  library(qvalue)
  y <- DGEList(counts = Pilot.data, group = status)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  if(method == "edgeR") p.value <- exactTest(y,dispersion="common")$table$P #Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
  delta = 1/y$common.dispersion
  ##56.37761

  GLM.fit.each <- function(z){
    y0 = z[which(status==group.name[1])]
    y1 = z[which(status==group.name[2])]
    f <- function(z,d) {
      delta = d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1 = d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(delta+y0)-lgamma(delta)-lgamma(y0+1)+y0*log(R/delta*exp(beta.0))-(y0+delta)*log(1+R/delta*exp(beta.0)))-sum(lgamma(delta+y1)-lgamma(delta)-lgamma(y1+1)+y1*log(R/delta*exp(beta.0+beta.1))-(y1+delta)*log(1+R/delta*exp(beta.0+beta.1)))
      return(Target.Func)
    }

    pt1 <- optim(c(0, 0), f, d = list(delta, R, y0, y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[2])
    var.beta.1 = (1/N0)*((1 + theta*r0)/(theta*R*r1) + (1+theta)/(theta*delta))
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }

  if(method == "GLM"){
    OUT.pars = t(apply(round(Pilot.data), 1, GLM.fit.each))
    colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")
    p.value = 2*pnorm(-abs(OUT.pars[,3]))
  }

  q.value = p.adjust(p.value, method = "BH")

  mean.by.group = apply(Pilot.data,1,function(x) tapply(x, status, mean))
  fold.change = (mean.by.group[1, ] + 1)/(mean.by.group[2, ] + 1)
  lfc = abs(log2(fold.change))

  #p.value.2 = p.value
  if(resample.DE){ #p-value calibration
    resampling.DE = which(lfc < log2(fc.cut) & p.value < p.cut)
    if(length(resampling.DE)){
      p.value.unif = runif(length(resampling.DE), 0, 1)
      p.value[resampling.DE] = p.value.unif
    }
  }

  return(list(p = p.value, fold.change = fold.change, dispersion = delta,
              Number.DE = sum(q.value < FDR), method = method))
}

Each.EDR <- function(Result, target.N = NULL, target.N.true = target.N, method = "SeqDesign", True = NULL, variation.True = F, sd.True = NULL, pilot.n = 2, method.mean = "mean", True.upper = NULL, True.lower = NULL, output.MSE = F, Power.max = NULL){
  #Result = HIP.n.2
  #Result = HIP.n.8
  if(0){
    l = sapply(Result, function(x) x[[2]][1])
    r = sapply(Result, function(x) x[[2]][2])
    s = sapply(Result, function(x) x[[2]][3])

    lambda = l[4]; r = r[4]; s = s[4]
    p.value.mod = p.value = Result[[4]][[3]]
    Fitted.Model = MixtureModel.two.beta(p.value)
    posterior=lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda) #prob be non-DE

    lambda = as.numeric(Fitted.Model[1])
    r=as.numeric(Fitted.Model[2])
    s=as.numeric(Fitted.Model[3])
    r2=as.numeric(Fitted.Model[4])
    s2=as.numeric(Fitted.Model[5])
    posterior.2B = (lambda*dbeta(p.value.mod,r2,s2))/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda*dbeta(p.value.mod,r2,s2)) #prob be non-DE


    hist(p.value.mod, xlab = "p-value", main = "")
    par(new = T)
    plot(p.value.mod, 1-posterior, cex = .5, pch = 15, col = 2, axes = F, ylab = "", xlab = ""); axis(4)
    points(p.value.mod, 1-posterior.2B, cex = .5, pch = 15, col = 4, axes = F, ylab = "", xlab = "")

    pdf(paste("Check p-value hist 2-12.pdf", sep = ""), width = 14)
    for(n in c(2, 4, 6, 8, 10)){
      Result = get(paste("HIP.n.", n, sep = ""))
      p.value = lapply(Result,function(y) y[[3]])
      pi = sapply(Result, function(x) x[[2]][1])
      model.CDD = sapply(p.value, function(x) MixtureModel.Fittting.pi0(x))
      pi.CDD = model.CDD[1,]
      posterior = lapply(Result, function(x) {
        Fitted.Model = x[[2]]
        p.value.mod = x[[3]]

        lambda = as.numeric(Fitted.Model[1])
        r=as.numeric(Fitted.Model[2])
        s=as.numeric(Fitted.Model[3])
        r2=as.numeric(Fitted.Model[4])
        s2=as.numeric(Fitted.Model[5])
        posterior.2B = (lambda*dbeta(p.value.mod,r2,s2))/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda*dbeta(p.value.mod,r2,s2)) #prob be non-DE
        return(posterior.2B)
      })

      posterior.CDD = lapply(1:10, function(i) {
        Fitted.Model = model.CDD[,i]
        p.value.mod = p.value[[i]]

        lambda = as.numeric(Fitted.Model[1])
        r=as.numeric(Fitted.Model[2])
        s=as.numeric(Fitted.Model[3])
        posterior = lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda)
        return(posterior)
      })

      # pdf(paste("Check p-value hist N", n, ".pdf", sep = ""), width = 14)
      par(mfrow = c(3, 4))
      for(i in 1:10){
        # hist(p.value[[i]], main = paste(round(pi[i], 3), round(pi.CDD[i], 3)))
        hist(p.value[[i]], main = paste(round(pi[i], 3), ", ", round(pi.CDD[i], 3), sep = ""))
        par(new = T)
        plot(p.value[[i]], 1-posterior[[i]], cex = .5, pch = 15, col = 2, axes = F, ylab = "", xlab = ""); axis(4)
        points(p.value[[i]], 1-posterior.CDD[[i]], cex = .5, pch = 15, col = 4, axes = F, ylab = "", xlab = "")
      }
      plot(1, type = "n", axes = F, xlab = "", ylab = ""); text(1, 1, paste("N = ", n, sep = ""), cex = 3)
      plot(1, type = "n", axes = F, xlab = "", ylab = "")
    }
    Result = get(paste("HIP.n.", 12, sep = ""))
    pi = Result[[2]][1]

    Fitted.Model = Result[[2]]
    p.value.mod = p.value = Result[[3]]
    lambda = as.numeric(Fitted.Model[1])
    r=as.numeric(Fitted.Model[2])
    s=as.numeric(Fitted.Model[3])
    r2=as.numeric(Fitted.Model[4])
    s2=as.numeric(Fitted.Model[5])
    posterior = (lambda*dbeta(p.value.mod,r2,s2))/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda*dbeta(p.value.mod,r2,s2)) #prob be non-DE

    model.CDD = MixtureModel.Fittting.pi0(p.value)
    pi.CDD = lambda = as.numeric(model.CDD[1])
    r=as.numeric(model.CDD[2])
    s=as.numeric(model.CDD[3])

    posterior.CDD = lambda/(dbeta(p.value.mod,r,s)*(1-lambda)+lambda)

    hist(p.value, main = paste(round(pi, 3), round(pi.CDD, 3)))
    par(new = T)
    plot(p.value, 1-posterior, cex = .5, pch = 15, col = 2, axes = F, ylab = "", xlab = ""); axis(4)
    points(p.value, 1-posterior.CDD, cex = .5, pch = 15, col = 4, axes = F, ylab = "", xlab = "")

    plot(1, type = "n", axes = F, xlab = "", ylab = ""); text(1, 1, paste("N = ", 12, sep = ""), cex = 3)

    dev.off()
  }

  if(method == "SeqDesign"){
    p.value=lapply(Result,function(y) y[[3]])
    Result=lapply(Result,function(y) y[[1]])
    # Result=Result[(sapply(Result,length)>1)&(sapply(p.value,function(x) mean(x)<0.5))]

    # Result=Result[(sapply(Result,length)>1)]
    # EDR=sapply(Result,function(z) (Reduce('+', z)/length(z))[3,])
    EDR=sapply(Result,function(z) apply(sapply(z, function(x) x[[1]][3,]), 1, median))
  } else{
    sample.size = target.N
    EDR = Result
    if(class(EDR) == "list") EDR = do.call("cbind", EDR)
    rownames(EDR) = sample.size
  }
  # EDR=cbind(EDR,matrix(0,nrow(EDR),length(Result)-ncol(EDR)))

  library(ggplot2)
  library(drc)
  library(DEoptim)

  mean.EDR=apply(EDR,1,function(x) mean(na.omit(x)))
  sd.EDR=apply(EDR,1,function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))

  if(is.null(True)){
    # mean.True=True
    # sd.True=0
    EDR.CI=cbind(mean.EDR,sapply(mean.EDR-1.96*sd.EDR,function(z) max(z,0)),sapply(mean.EDR+1.96*sd.EDR,function(z) min(z,1)),c(as.numeric(names(mean.EDR))+0.5))
    colnames(EDR.CI)=c("EDR","lower","upper","N")
    # True.CI=cbind(mean.True,sapply(mean.True-1.96*sd.True,function(z) max(z,0)),sapply(mean.True+1.96*sd.True,function(z) min(z,1)),c(12, 20, 30, 40, 50, 100))
    # colnames(True.CI)=c("EDR","lower","upper","N")
    # row.names(True.CI)=c(12, 20, 30, 40, 50, 100)
    # data=data.frame(cbind(rbind(EDR.CI,True.CI),Type=c(rep("Predicted EDR", nrow(EDR.CI)), rep("True EDR", nrow(True.CI)))))
    data=data.frame(cbind(EDR.CI,Type=c(rep("Predicted EDR", nrow(EDR.CI)))))
    # data=data.frame(EDR.CI)
    pd <- position_dodge(width=0.2)
    data[,1]=as.numeric(as.matrix(data[,1]))
    data[,2]=as.numeric(as.matrix(data[,2]))
    data[,3]=as.numeric(as.matrix(data[,3]))
    data[,4]=as.numeric(as.matrix(data[,4]))
    out=ggplot(data, aes(N,EDR, color=Type)) + ggtitle(paste("N = ", pilot.n, sep = "")) +
      geom_point(aes(shape=Type),size=3, position=pd) +
      # scale_color_manual(name="Method",values=c("red","blue")) +
      scale_color_manual(name = "Method", values=c("red")) +
      scale_shape_manual(name = "Method", values=c(17)) +
      # scale_fill_manual(values=c("red"))+
      theme_bw() +theme(legend.position="none") +
      geom_errorbar(aes(ymin=lower,ymax=upper),width=1)+ ylim(0, 1)+geom_line(aes(group=Type))+ geom_ribbon(data=data, aes(ymin = lower,ymax=upper), alpha=0.4)
  } else {
    mean.True = if(variation.True) apply(True, 1, function(x) if(method.mean == "mean") mean(na.omit(x)) else if(method.mean == "median") median(na.omit(x))) else True
    # names(mean.True) = c(12,20,30,40,50,100)
    names(mean.True) = target.N.true
    if(is.null(sd.True)){
      sd.True = if(variation.True) apply(True, 1,function(x) sd(na.omit(x))/sqrt(length(na.omit(x)))) else 0
    }

    if(output.MSE){
      Each.N <- function(EDR.each){
        Fit.coef=drm(EDR.each~sample.size, fct = LL.5())$coefficients
        Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
        if(Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(10000)-log(Fit.coef[4]))))^Fit.coef[5]<Power.max){N.hat.true=10000}
        else{
          N.hat.true=round(uniroot(Find.N, c(0, 10000))$root)
        }
        return(N.hat.true)
      }

      N.hat=apply(EDR,2,function(x) Each.N(x))

      Fit.coef=drm(mean.True~sample.size, fct = LL.5())$coefficients
      Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
      N.hat.true=round(uniroot(Find.N, c(0, 100))$root)

      # RMSE = sqrt(mean(((mean.EDR-mean.True)^2)[sample.size>10]))
      RMSE = sqrt(mean(((mean.EDR-mean.True)^2)))

      RMSE.N = sqrt(mean((N.hat-N.hat.true)^2))
    }
    EDR.CI=cbind(mean.EDR,sapply(mean.EDR-1.96*sd.EDR,function(z) max(z,0)),sapply(mean.EDR+1.96*sd.EDR,function(z) min(z,1)),c(as.numeric(names(mean.EDR))+0.5))
    colnames(EDR.CI)=c("EDR","lower","upper","N")
    if(is.null(True.upper)){
      True.CI=cbind(mean.True,sapply(mean.True-1.96*sd.True,function(z) min(max(z,0), 1)),sapply(mean.True+1.96*sd.True,function(z) max(min(z,1), 0)),c(as.numeric(names(mean.True))+0.5))
    } else{
      True.CI=cbind(mean.True, True.upper, True.lower, c(as.numeric(names(mean.True))+0.5))
    }
    colnames(True.CI)=c("EDR","lower","upper","N")
    row.names(True.CI)=c(as.numeric(names(mean.True))+0.5)
    data=data.frame(cbind(rbind(EDR.CI,True.CI),Type=c(rep("Predicted EDR", nrow(EDR.CI)), rep("True EDR", nrow(True.CI)))))
    # data=data.frame(EDR.CI)
    pd <- position_dodge(width=0.2)
    data[,1]=as.numeric(as.matrix(data[,1]))
    data[,2]=as.numeric(as.matrix(data[,2]))
    data[,3]=as.numeric(as.matrix(data[,3]))
    data[,4]=as.numeric(as.matrix(data[,4]))
    out=ggplot(data, aes(N,EDR, color=Type)) + xlim(0, 101) +
      geom_point(aes(shape=Type),size=3, position=pd) +
      scale_color_manual(name="Method",values=c("black", "gray40")) +
      scale_shape_manual(name = "Method", values=c(15, 17)) +
      theme_bw() +
      theme(axis.title = element_blank())+
      theme(axis.text=element_blank())+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
      theme(axis.ticks = element_blank())+
      scale_fill_manual(values=c("red"))+
      theme(legend.position="none") +
      geom_errorbar(aes(ymin=lower,ymax=upper, linetype = Type),width=1)+ ylim(0, 1) + geom_line(aes(group=Type, linetype = Type))+ geom_ribbon(data=data, aes(ymin = lower,ymax=upper, linetype = Type), alpha=0.4)
  }
  # return(out)
  if(output.MSE) return(list(out, RMSE, RMSE.N)) else return(out)
}

Each.EDR.other <- function(EDR, True, Power.max){
  sample.size=c(5,10,20,30,40,50,100)
  Each.N <- function(EDR.each){
    Fit.coef=drm(EDR.each~sample.size, fct = LL.5())$coefficients
    Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
    #par(mfrow=c(1,2))
    #plot(sample.size,EDR.each)
    #points(sample.size,Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(sample.size)-log(Fit.coef[4]))))^Fit.coef[5],type="l",ylim=c(0,1),ylab="EDR",xlab="N")

    library(DEoptim)
    #fitting=curve.fitting(sample.size,EDR.each)
    #plot(sample.size,EDR.each)
    #points(sample.size,1-fitting[1]*sample.size^(-fitting[2]),type="l")

    if(Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(10000)-log(Fit.coef[4]))))^Fit.coef[5]<Power.max){N.hat.true=10000}
    else{
      N.hat.true=round(uniroot(Find.N, c(0, 10000))$root)
    }
    return(N.hat.true)
  }

  N.hat=apply(EDR,2,function(x) Each.N(x))

  mean.EDR=apply(EDR,1,function(x) mean(na.omit(x)))
  sd.EDR=apply(EDR,1,function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))

  mean.True=sapply(True,function(x) mean(na.omit(x)))
  sd.True=0

  RMSE=sqrt(mean(((mean.EDR-mean.True)^2)[sample.size>10]))
  RMSE.old=sqrt(mean(((mean.EDR-mean.True)^2)))
  RMSE.new=mean(apply((EDR-mean.True),1,function(x) sqrt(mean(x^2)))[c(5,10,20,30,40,50,100)>10])

  RMSE.N=sqrt(mean((N.hat-N.hat.true)^2))
  EDR.CI=cbind(mean.EDR,sapply(mean.EDR-1.96*sd.EDR,function(z) max(z,0)),sapply(mean.EDR+1.96*sd.EDR,function(z) min(z,1)),c(5,10,20,30,40,50,100))
  colnames(EDR.CI)=c("EDR","lower","upper","N")
  True.CI=cbind(mean.True,sapply(mean.True-1.96*sd.True,function(z) max(z,0)),sapply(mean.True+1.96*sd.True,function(z) min(z,1)),c(5,10,20,30,40,50,100))
  colnames(True.CI)=c("EDR","lower","upper","N")
  row.names(True.CI)=c(5,10,20,30,40,50,100)
  data=data.frame(cbind(rbind(EDR.CI,True.CI),Type=rep(c("Predicted EDR","True EDR"),each=7)))
  pd <- position_dodge(width=0.2)
  data[,1]=as.numeric(as.matrix(data[,1]))
  data[,2]=as.numeric(as.matrix(data[,2]))
  data[,3]=as.numeric(as.matrix(data[,3]))
  data[,4]=as.numeric(as.matrix(data[,4]))
  out=ggplot(data, aes(N,EDR, color=Type)) + xlim(0, 101) +
    geom_point(aes(shape=Type),size=3, position=pd) +
    scale_color_manual(name="Method",values=c("black", "gray40")) +
    scale_shape_manual(name = "Method", values=c(15, 17)) +
    theme_bw() +
    theme(axis.title = element_blank())+
    theme(axis.text=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.ticks = element_blank())+
    scale_fill_manual(values=c("red"))+
    theme(legend.position="none") +
    geom_errorbar(aes(ymin=lower,ymax=upper, linetype = Type),width=1)+ ylim(0, 1) + geom_line(aes(group=Type, linetype = Type))+ geom_ribbon(data=data, aes(ymin = lower,ymax=upper, linetype = Type), alpha=0.4)
  return(list(out, RMSE, RMSE.N, RMSE.new, RMSE.old))
}

Each.EDR.SeqDEsign <- function(Result, True,method, Power.max){
  for (i in 1:length(Result)){
    for (j in 1:length(Result[[i]])){
      length.list=sapply(Result[[i]][[j]],function(x) length(x)==7)
      if(sum(length.list)>0){Result[[i]][[j]]=Result[[i]][[j]][!length.list]}
    }
  }
  Method <- lapply(Result,function(y) sapply(y[1:5],function(x) Reduce('+', x)/length(x)))
  Result.by.Method=lapply(1:5,function(x) sapply(Method,function(y) y[,x]))
  names(Result.by.Method)=c("3Par(PT)","Storey(PT)","EM(PT)","BUM(PT)","CDD(PT)")
  for(i in 1:length(Result.by.Method)){
    row.names(Result.by.Method[[i]]) = rep(c(5,10,20,30,40,50,100),each=5)
  }


  FDR.predicted=lapply(Result.by.Method,function(x) x[1+c(0:6)*5,])
  EDR.predicted=lapply(Result.by.Method,function(x) x[3+c(0:6)*5,])

  EDR=as.vector(EDR.predicted[[method]])
  sample.size=rep(c(5,10,20,30,40,50,100),50)

  sample.size.N=c(5,10,20,30,40,50,100)
  Each.N <- function(EDR.each){
    Fit.coef=drm(EDR.each~sample.size.N, fct = LL.5())$coefficients
    Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
    #plot(1:2000,Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(1:2000)-log(Fit.coef[4]))))^Fit.coef[5])
    if(Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(10000)-log(Fit.coef[4]))))^Fit.coef[5]<Power.max){N.hat.true=10000}
    else{
      N.hat.true=round(uniroot(Find.N, c(0, 10000))$root)
    }
    return(N.hat.true)
  }
  N.hat=apply(EDR.predicted[[method]],2,function(x) Each.N(x))
  RMSE.N=sqrt(mean((N.hat-N.hat.true)^2))

  Power.max=Power.max
  Fit.coef=drm(EDR.predicted[[method]][,1]~sample.size.N, fct = LL.5())$coefficients
  Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
  plot(5:100,Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(5:100)-log(Fit.coef[4]))))^Fit.coef[5],type="l",ylim=c(0,1),ylab="EDR",xlab="N")
  for (i in 2:50){
    Fit.coef=drm(EDR.predicted[[method]][,i]~sample.size.N, fct = LL.5())$coefficients
    Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
    points(5:100,Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(5:100)-log(Fit.coef[4]))))^Fit.coef[5],type="l")
  }
  abline(h=Power.max,col="red",lty=2)
  for(i in 1:50){
    abline(v=N.hat[i],col="blue")
  }
  library(ggplot2)
  a <- data.frame(N=sample.size,EDR=EDR)
  b=data.frame(N=unlist(lapply(1:length(True),function(x) rep(c(5,10,20,30,40,50,100)[x],length(True[[x]])))),EDR=unlist(True))
  data=cbind(rbind(a,b),c(rep("Predicted",nrow(a)),rep("True",nrow(b))))
  colnames(data)[3]="group"
  #out=ggplot(data, aes(x=N,y=EDR, group = group, colour = group)) +
  #  geom_point(aes(colour = group),position = "jitter") +
  #  geom_smooth(aes(fill = group)) + ylim(0, 1)+ ggtitle(paste(paste(c("FoldChange:",foldchange,";Dispersion:",Dispersion),collapse="")))


  mean.EDR=apply(EDR.predicted[[method]],1,function(x) mean(na.omit(x)))
  sd.EDR=apply(EDR.predicted[[method]],1,function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))

  mean.True=sapply(True,function(x) mean(na.omit(x)))
  sd.True=0

  RMSE=sqrt(mean(((mean.EDR-mean.True)^2)[c(5,10,20,30,40,50,100)>10]))
  RMSE.new=mean(apply((EDR.predicted[[method]]-mean.True),1,function(x) sqrt(mean(x^2)))[c(5,10,20,30,40,50,100)>10])

  EDR.CI=cbind(mean.EDR,sapply(mean.EDR-1.96*sd.EDR,function(z) max(z,0)),sapply(mean.EDR+1.96*sd.EDR,function(z) min(z,1)),c(5,10,20,30,40,50,100)+0)
  colnames(EDR.CI)=c("EDR","lower","upper","N")
  True.CI=cbind(mean.True,sapply(mean.True-1.96*sd.True,function(z) max(z,0)),sapply(mean.True+1.96*sd.True,function(z) min(z,1)),c(5,10,20,30,40,50,100)-0)
  colnames(True.CI)=c("EDR","lower","upper","N")
  row.names(True.CI)=c(5,10,20,30,40,50,100)
  data=data.frame(cbind(rbind(EDR.CI,True.CI),Type=rep(c("Predicted EDR","True EDR"),each=7)))
  pd <- position_dodge(width=0.2)
  data[,1]=as.numeric(as.matrix(data[,1]))
  data[,2]=as.numeric(as.matrix(data[,2]))
  data[,3]=as.numeric(as.matrix(data[,3]))
  data[,4]=as.numeric(as.matrix(data[,4]))
  out=ggplot(data, aes(N,EDR, color=Type)) + xlim(0, 101) +
    geom_point(aes(shape=Type),size=3, position=pd) +
    scale_color_manual(name="Method",values=c("black", "gray40")) +
    scale_shape_manual(name = "Method", values=c(15, 17)) +
    theme_bw() +
    theme(axis.title = element_blank())+
    theme(axis.text=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.ticks = element_blank())+
    scale_fill_manual(values=c("red"))+
    theme(legend.position="none") +
    geom_errorbar(aes(ymin=lower,ymax=upper, linetype = Type),width=1)+ ylim(0, 1) + geom_line(aes(group=Type, linetype = Type))+ geom_ribbon(data=data, aes(ymin = lower,ymax=upper, linetype = Type), alpha=0.4)
  return(list(out,RMSE,RMSE.N,N.hat,RMSE.new))
}

Each.EDR.SeqDEsign.2 <- function(Result, True, method, Power.max){
  p.value=lapply(Result,function(y) y[[3]])
  Result=lapply(Result,function(y) y[[1]])
  EDR.matrix=sapply(Result,function(z) apply(sapply(z, function(x) x[3,]), 1, median))

  EDR=as.vector(EDR.matrix)
  sample.size=rep(c(5,10,20,30,40,50,100),50)

  if(0){
  sample.size.N=c(5,10,20,30,40,50,100)
  Each.N <- function(EDR.each){
    if(all(EDR.each != 1)){
      Fit.coef=drm(EDR.each~sample.size.N, fct = LL.5())$coefficients
      Find.N <- function(x) Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(x)-log(Fit.coef[4]))))^Fit.coef[5]-Power.max
      #plot(1:2000,Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(1:2000)-log(Fit.coef[4]))))^Fit.coef[5])
      if(Fit.coef[2]+(Fit.coef[3]-Fit.coef[2])/(1+exp(Fit.coef[1]*(log(10000)-log(Fit.coef[4]))))^Fit.coef[5]<Power.max){N.hat.true=10000}
      else{
        N.hat.true=round(uniroot(Find.N, c(0, 10000))$root)
      }
    } else N.hat.true = min(sample.size)
    return(N.hat.true)
  }
  N.hat=apply(EDR.matrix,2,function(x) Each.N(x))
  RMSE.N=sqrt(mean((N.hat-N.hat.true)^2))
  }
  mean.EDR=apply(EDR.matrix,1,function(x) mean(na.omit(x)))
  sd.EDR=apply(EDR.matrix,1,function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))

  mean.True=sapply(True,function(x) mean(na.omit(x)))
  sd.True=0

  RMSE=sqrt(mean(((mean.EDR-mean.True)^2)[c(5,10,20,30,40,50,100)>10]))
  RMSE.new=mean(apply((EDR.matrix-mean.True),1,function(x) sqrt(mean(x^2)))[c(5,10,20,30,40,50,100)>10])

  EDR.CI=cbind(mean.EDR,sapply(mean.EDR-1.96*sd.EDR,function(z) max(z,0)),sapply(mean.EDR+1.96*sd.EDR,function(z) min(z,1)),c(5,10,20,30,40,50,100)+0)
  colnames(EDR.CI)=c("EDR","lower","upper","N")
  True.CI=cbind(mean.True,sapply(mean.True-1.96*sd.True,function(z) max(z,0)),sapply(mean.True+1.96*sd.True,function(z) min(z,1)),c(5,10,20,30,40,50,100)-0)
  colnames(True.CI)=c("EDR","lower","upper","N")
  row.names(True.CI)=c(5,10,20,30,40,50,100)
  data=data.frame(cbind(rbind(EDR.CI,True.CI),Type=rep(c("Predicted EDR","True EDR"),each=7)))
  pd <- position_dodge(width=0.2)
  data[,1]=as.numeric(as.matrix(data[,1]))
  data[,2]=as.numeric(as.matrix(data[,2]))
  data[,3]=as.numeric(as.matrix(data[,3]))
  data[,4]=as.numeric(as.matrix(data[,4]))
  out=ggplot(data, aes(N,EDR, color=Type)) + xlim(0, 101) +
    geom_point(aes(shape=Type),size=3, position=pd) +
    scale_color_manual(name="Method",values=c("black", "gray40")) +
    scale_shape_manual(name = "Method", values=c(15, 17)) +
    theme_bw() +
    theme(axis.title = element_blank())+
    theme(axis.text=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.ticks = element_blank())+
    scale_fill_manual(values=c("red"))+
    theme(legend.position="none") +
    geom_errorbar(aes(ymin=lower,ymax=upper, linetype = Type),width=1)+ ylim(0, 1) + geom_line(aes(group=Type, linetype = Type))+ geom_ribbon(data=data, aes(ymin = lower,ymax=upper, linetype = Type), alpha=0.4)
  # return(list(out,RMSE,RMSE.N,N.hat,RMSE.new))
  return(list(out,RMSE,RMSE.new))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

curve.fitting = function(n, EDR, error.type = "original", EDR.bound = 10000, b.bound = 10000, c.lower = 0, b.lower = 0){
  library(DEoptim)
  lkh<-function(x){
    #       return(sum((-x[2]*r^(-x[3])+x[1]-EDR)^2))
    # if(any(-x[1]*n^(-x[2])+x[3] < 0)) print(x)
    if(error.type == "original") return(sum((-x[1]*n^(-x[2])+x[3]-EDR)^2))
    else if(error.type == "log"){
      if(any(-x[1]*n^(-x[2])+x[3] < 0)) return(1000000)
      else return(sum((log(-x[1]*n^(-x[2])+x[3])-log(EDR))^2))
    }
  }
  #     outOptim<-DEoptim(lkh,lower=c(0, 0, 0),upper=c(5, 5, 5), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
  #   outOptim<-DEoptim(lkh,lower=c(0, 0),upper=c(5, 5), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
  outOptim<-DEoptim(lkh, lower = c(b.lower, c.lower, 1500),upper=c(b.bound, 10, EDR.bound), DEoptim.control(trace = F, itermax = 10000))
  ###stochastic fitting
  return(outOptim$optim$bestmem)
}

curve.fitting.2 = function(n, EDR, error.type = "original", EDR.bound = 10000, b.bound = 10^7, b.lower = 0, c.lower = 0, c.upper = 10){
  library(DEoptim)
  lkh<-function(x){
    #       return(sum((-x[2]*r^(-x[3])+x[1]-EDR)^2))
    # return(sum((-x[1]*n^(-x[2])-x[4]*n^(-x[5]^2)+x[3]-EDR)^2))
    # return(sum((-x[1]*n^(-x[2])+x[3]-EDR)^2))
    if(error.type == "original") return(sum((-x[1]*x[2]^(-n)+x[3]-EDR)^2))
    else if(error.type == "log") return(sum((log(-x[1]*n^(-x[2])+x[3])-log(EDR))^2))
  }
  #     outOptim<-DEoptim(lkh,lower=c(0, 0, 0),upper=c(5, 5, 5), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
  #   outOptim<-DEoptim(lkh,lower=c(0, 0),upper=c(5, 5), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
  # outOptim<-DEoptim(lkh, lower = c(0, c.lower, 0, 0, c.lower),upper=c(b.bound, 10, EDR.bound, b.bound, 10), DEoptim.control(trace = F, itermax = 10000))###stochastic fitting
  outOptim<-DEoptim(lkh, lower = c(b.lower, c.lower, 0), upper=c(b.bound, c.upper, EDR.bound), DEoptim.control(trace = F, itermax = 20000))###stochastic fitting
  print(outOptim$optim$bestval)
  print(outOptim$optim$bestmem)
  # init = c(1000, 0, 9000)
  # # optim(init, fn = lkh, n = n, EDR = EDR, method= "L-BFGS-B", upper=c(b.bound, 10, EDR.bound),
  # #       lower=c(0, c.lower, 0))$par
  # a = optim(init, fn = lkh, n = n, EDR = EDR, method= "Nelder-Mead")
  # a$par
  # a$value
  print(outOptim$optim$bestval)
  return(outOptim$optim$bestmem)
}

curve.fitting.new = function(n, EDR, EDR.bound = 10000, b.bound = 10^7, c.lower = 0){
  lkh<-function(x){
    #       return(sum((-x[2]*r^(-x[3])+x[1]-EDR)^2))
    return(sum((-x[1]*n^(-x[2])+x[3]-EDR)^2))
  }

  # pdf("Log fitting model.pdf", width = 14, height = 14)
  pdf("Log fitting model standardize residual.pdf", width = 14, height = 14)
  par(mfrow = c(4, 5), mar = c(3, 3, 2, 1)+ 0.1, mgp = c(2, 0.6, 0), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2, tck = -0.005, font.lab = 2)
  z = seq(2000, 10000, 100)
  result = lapply(z, function(a){

    y = log(a - EDR)
    x = -log(n)
    plot(x, y, xlab = "-log(N)", ylab = "log(a - #DE)", main = a)
    model = lm(y~x)
    lines(x, fitted(model), col = 2)
    cf = coef(model)
    fitted = a - exp(cf[1])*n^(-cf[2])
    res = sqrt(sum((EDR - fitted)^2))
    # res = sqrt(sum((log(EDR) - log(fitted))^2))
    # sum.square.res = sum(rstandard(model)^2)
    # return(list(coef = coefficients(model), sum.square.res))
    return(list(coef = cf, residual = res))
  })
  par(mfrow = c(1, 1))
  plot(z, sapply(result, function(x) x[[2]]), type = "l", xlab = "Converging DE number", ylab = "Goodness of fitting")
  abline(h = 0, col = 2)
  dev.off()
  outOptim<-DEoptim(lkh, lower = c(0, c.lower, 0),upper=c(b.bound, 10, EDR.bound), DEoptim.control(trace = F, itermax = 10000))###stochastic fitting
  return(outOptim$optim$bestmem)
}

fitted.curve = function(n, fit.model) fit.model[3]-fit.model[1]*n^(-fit.model[2])
# fitted.curve.2 = function(n, fit.model) fit.model[3]-fit.model[1]*n^(-fit.model[2])-fit.model[4]*n^(-fit.model[5]^2)
fitted.curve.2 = function(n, fit.model) fit.model[3]-fit.model[1]*fit.model[2]^(-n)

DE.plot = function(N.DE, FDR = 0.05, LFC = 0.4, dataset = "ER", y.min = 500,
                   y.max = 2500, q.full = NULL, lfc.full = NULL,
                   DE.num.FDR.adjust = F){
  # num.sig = sapply(1:14, function(i){ #for different target sample size
  # num.sig = sapply(8:14, function(i){ #for different target sample size
  # num.sig = sapply(6:14, function(i){ #for different target sample size
  #   x = if(i <= 6){
  #     get(paste("Stage.Predict.", i, ".2B", sep = ""))
  #   } else {
  #     get(paste("Stage.DE.", i, sep = ""))
  #   }
  #   num.DE = sapply(x, function(y) y$Number.DE)
  #
  # })

  num.sig = sapply(1:length(N.DE), function(i){
    x = get(paste(dataset, ".DE.", i, sep = ""))
    # num.DE = sapply(x, function(y) y$Number.DE)
    if(0){
      i = 10; j = 14
      i = 1; j = 10
      x = get(paste(dataset, ".DE.", i, sep = ""))

      z = get(paste(dataset, ".DE.", j, sep = ""))

      # pdf(paste("lfc vs q-value ", i, " vs ", j, ".pdf", sep = ""), width = 21, height = 14)
      pdf(paste("lfc vs q-value ", i, " vs ", j, " edgeR.pdf", sep = ""), width = 21, height = 14)
      par(mfrow = c(4, 5), mar = c(3, 3, 2, 1)+ 0.1, mgp = c(2, 0.6, 0), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2, tck = -0.005, font.lab = 2)
      for(i in 1:20){
        y = x[[i]]
        z2 = z[[i]]

        lfc = log2(y$fold.change); lfc2 = log2(z2$fold.change)
        q = p.adjust(y$p, method = "BH"); q2 = p.adjust(z2$p, method = "BH")

        plot(lfc, -log10(q),
             ylim = c(0, 50),
             xlim = c(-2, 2),
             xlab = "lfc", cex = .1, pch = 16, col = "#00000080")
        points(lfc2, -log10(q2), col = "#A020F080", cex = .1, pch = 16)
        abline(v = c(-LFC, LFC), h = -log10(0.05), col = 2)
        # print(sum(abs(lfc) > 0.26 & q < 0.05))
      }
      for(i in 1:20){
        y = x[[i]]
        z2 = z[[i]]

        lfc = log2(y$fold.change); lfc2 = log2(z2$fold.change)
        q = p.adjust(y$p, method = "BH"); q2 = p.adjust(z2$p, method = "BH")

        a = 5
        lfc[lfc > a] = a; lfc[lfc < -a] = -a
        lfc2[lfc2 > a] = a; lfc2[lfc2 < -a] = -a
        hist(lfc, breaks = seq(-a, a, length.out = 100), col = "#FF000080")
        hist(lfc2, breaks = seq(-a, a, length.out = 100), col = "#0000FF80", add = T)
        abline(v = c(-LFC, LFC), col = 2)

      }

      for(i in 1:20){
        y = x[[i]]
        z2 = z[[i]]

        lfc = log2(y$fold.change); lfc2 = log2(z2$fold.change)
        q = -log10(p.adjust(y$p, method = "BH")); q2 = -log10(p.adjust(z2$p, method = "BH"))

        a = 5
        q[q > a] = a
        q2[q2 > a] = a
        hist(q, breaks = seq(0, a, length.out = 40), col = "#FF000080", xlab = "-log10(q)", main = "")
        hist(q2, breaks = seq(0, a, length.out = 40), col = "#0000FF80", add = T)
      }

      dev.off()
    }


    num.DE = sapply(x, function(y){
      lfc = log2(y$fold.change)
      q = p.adjust(y$p, method = "BH")
      a = sum(q < FDR & abs(lfc) > LFC)
      if(DE.num.FDR.adjust) a = a * (1-FDR)
      return(a)
    })
  }); num.sig

  # num.sig = cbind(num.sig, num.sig[, 12], num.sig[, 12], num.sig[, 12], num.sig[, 12], num.sig[, 12], num.sig[, 12])
  # num.sig = cbind(num.sig, num.sig[, 16:30], num.sig[, 16:30])
  # N.DE = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100:150)
  if(!is.null(q.full))DE.full = sum(q.full < FDR & abs(lfc.full) > LFC)

  median.DE = apply(num.sig, 2, median)
  par(mfrow = c(1, 1))
  ylim = if(!is.null(q.full)) c(min(num.sig), max(DE.full, num.sig)) else c(y.min, y.max)
  plot(rep(N.DE, each = nrow(num.sig)), num.sig, main = paste("FDR = ", FDR, ", LFC = ", LFC, sep = ""), xlab = expression("N"[0]), ylab = "Number of DE", ylim = ylim)
  lines(N.DE, median.DE, col = 2)

  if(!is.null(q.full)){
    abline(h = DE.full, col = 3)
    return(list(Num.DE = num.sig, Median.DE = median.DE, DE.full = sum(q.full < FDR & abs(lfc.full) > LFC)))
  } else return(list(Num.DE = num.sig, Median.DE = median.DE))
}

DE.fit.plot = function(N.DE, median.DE, num.sig, drop.index = NULL, DE.bound = 10000,
                       DE.full = NULL, c.lower = 0, b.lower = 0, b.bound = 10^7, FDR = FDR, LFC = LFC,
                       method = 1, error.type = "original", add = F, y.min = 500, y.max = 2500){
  # drop.index = c(1:3, 5, 7, 9)
  if(!is.null(drop.index)){
    N.DE = N.DE[-drop.index]
    median.DE = median.DE[-drop.index]
    num.sig = num.sig[, -drop.index]
  }

  if(method == 1){
    fitting = curve.fitting(N.DE, median.DE, error.type = error.type,
                            EDR.bound = DE.bound, b.bound = b.bound, b.lower = b.lower,
                            c.lower = c.lower); print(fitting)
    fit.value = fitted.curve(N.DE, fitting); fit.value
  } else {
    fitting = curve.fitting.2(N.DE, median.DE, error.type = error.type,
                              EDR.bound = DE.bound, b.bound = b.bound, b.lower = b.lower,
                              c.lower = c.lower); print(fitting)
    fit.value = fitted.curve.2(N.DE, fitting); fit.value
  }

  if(!add){
    par(mfrow = c(1, 1))
    ylim = if(!is.null(DE.full)) c(min(num.sig), max(DE.full, num.sig)) else c(y.min, y.max)
    plot(rep(N.DE, each = nrow(num.sig)), num.sig,
         main = paste("FDR = ", FDR, ", LFC = ", LFC, ", converge at ",
                      round(fitting[3]), ", b = ", round(fitting[1]),
                      ", c = ", round(fitting[2], 3), sep = ""),
         xlab = expression("N"[0]), ylab = "Number of DE", ylim = ylim)
    abline(h = fitting[3], col = 3)
    lines(N.DE, median.DE, col = 2)
    lines(N.DE, fit.value, col = 4)
    if(!is.null(DE.full)){
      abline(h = DE.full, col = 3)
    }
  } else lines(N.DE, fit.value, col = "purple")
  return(fitting)
}


check.convergence.plot = function(model, method = 1, x.max = 10000, x.high = 10000, add = F){
  x = c(2, 4, 6, 8, seq(10, x.max, 10))
  if(!add){
    y = if(method == 1) fitted.curve(x, model) else fitted.curve.2(x, model)
    plot(x, y, #xlim = c(10, 200),
         xlab = expression("N"[0]), ylab = "DE number", main = round(model[3]),
         pch = 16, cex = .5,
         xlim = c(0, x.high),
         ylim = c(min(y), max(y)))
    abline(h = model[3], col = 2)
  } else points(x, fitted.curve.2(x, model), col = 2, pch = 16, cex = .5)
}

p.hist.compare.posterior.prob = function(p.value, model, method = "2B", breaks = 100, add = F, col = 2, add.col = "green", separate.DE.null = F, verbose = F, ...){
  lambda = as.numeric(model[1])
  r=as.numeric(model[2])
  s=as.numeric(model[3])
  if(method == "2B"){
    r2=as.numeric(model[4])
    s2=as.numeric(model[5])
    posterior = (lambda*dbeta(p.value,r2,s2))/(dbeta(p.value,r,s)*(1-lambda)+lambda*dbeta(p.value,r2,s2)) #prob be non-DE
  } else if(method == "CDD") posterior = lambda/(dbeta(p.value,r,s)*(1-lambda)+lambda)

  if(!add){
    if(separate.DE.null){
      DE_status_posterior = sapply(1:length(posterior),function(x) sample(c(FALSE,TRUE),1,prob = c(posterior[x],1-posterior[x]),replace=TRUE))
      p.DE=p.value[DE_status_posterior]
      p.nonDE=p.value[!DE_status_posterior]
      hist(p.DE, col = "#ADD8E680", main = "", xlab = "p-value", breaks = breaks, ...)
      hist(p.nonDE, col = "#90EE9080", add = T, breaks = breaks, ...)

      par(new = T)
      plot(p.value, 1-posterior, cex = .5, pch = 15, col = 2, axes = F, ylab = "", xlab = ""); axis(4)

    } else hist(p.value, xlab = "p-value", main = "", breaks = breaks, ...)
    par(new = T)
    plot(p.value, 1-posterior, cex = .5, pch = 15, col = col, axes = F, ylab = "", xlab = ""); axis(4)
  } else{
    points(p.value, 1-posterior, cex = .5, pch = 15, col = add.col)
  }

  if(verbose) return(1-posterior)
}

variable.pdf = function(p.value, model) {
  lambda = as.numeric(model[1])
  r=as.numeric(model[2])
  s=as.numeric(model[3])
  r2=as.numeric(model[4])
  s2=as.numeric(model[5])
  f = (dbeta(p.value,r,s)*(1-lambda)+lambda*dbeta(p.value,r2,s2)) #prob be non-DE
  return(f)
}

pdf.prop = function(p.value, model) {
  lambda = as.numeric(model[1])
  r=as.numeric(model[2])
  s=as.numeric(model[3])
  r2=as.numeric(model[4])
  s2=as.numeric(model[5])
  f1 = dbeta(p.value, r, s)*(1-lambda) #prob be non-DE
  f2 = lambda*dbeta(p.value, r2, s2)
  return(list(f1, f2))
}

draw.pdf = function(model, col = 4){
  lambda = as.numeric(model[1])
  r=as.numeric(model[2])
  s=as.numeric(model[3])
  r2=as.numeric(model[4])
  s2=as.numeric(model[5])
  f = function(p.value) (dbeta(p.value,r,s)*(1-lambda)+lambda*dbeta(p.value,r2,s2)) #prob be non-DE

  v = seq(0.0001, .99999, 0.01)
  par(new = T)
  plot(v, sapply(v, function(x) f(x)), col = col, axes = F, ylim = c(0, max(f(v))));
  axis(4)
}

Parameter.Estimate <- function(Normalized.HIP.data, R.each){
  ## filtering by mean counts
  n = ncol(Normalized.HIP.data) / 2
  Mean.count = apply(Normalized.HIP.data, 1, mean)
  Mean.count = R.each * Mean.count / mean(Mean.count)

  ## No filtering
  Normalized.HIP.data.filter = Normalized.HIP.data
  dim(Normalized.HIP.data.filter)

  conds = rep(c("HIV", "F344"), each = 12)

  mean.by.group = t(apply(Normalized.HIP.data.filter, 1, function(x) tapply(x, factor(conds), mean)))
  lfc = log2((mean.by.group[, 2] + 1) / (mean.by.group[, 1] + 1))

  return(list(MeanCount = Mean.count))
}

NB.simu.setting.gamma <- function(ngenes, par.setting){
  library(truncnorm)
  m1 = par.setting$Prop.DE
  mean.lfc = par.setting$lfc.mean.sd[1]
  sd.lfc = par.setting$lfc.mean.sd[2]
  boundary.lower.lfc = par.setting$lfc[1]
  boundary.upper.lfc = par.setting$lfc[2]

  # a = c(-Inf, boundary.upper.lfc)
  # b = c(boundary.lower.lfc, Inf)

  set.seed(12345)
  is_DE <- runif(ngenes) < m1

  ## two-sides
  # lfc = rtruncnorm(ngenes, a = a, b = b, mean = mean.lfc, sd = sd.lfc)
  lfc.positive = rtruncnorm(ngenes/2, a = boundary.upper.lfc, b = Inf, mean = mean.lfc, sd = sd.lfc)
  lfc.negative = rtruncnorm(ngenes/2, a = -Inf, b = boundary.lower.lfc, mean = -mean.lfc, sd = sd.lfc)
  lfc = c(lfc.positive, lfc.negative)

  lfc[!is_DE]=0

  q0 <- sample(par.setting$MeanCount, ngenes, replace = TRUE)

  q0A <- q0
  q0A[is_DE]=q0[is_DE] * 2^(lfc[is_DE]/2)
  q0B <- q0
  q0B[is_DE]=q0[is_DE] * 2^(-lfc[is_DE]/2)

  Data.list = list(
    q0A = q0A,
    q0B = q0B,
    Read = sum(q0),
    lfc = log2(q0A/q0B),
    is_DE = is_DE,
    dipsersion = par.setting$dispersion
  )
  return(Data.list)
}

Generate.Data.Based.On.Setting <- function(n, Model){
  dispersion = Model$dipsersion
  ngenes = length(Model$q0A)
  conds <- c(rep("Case", n), rep("Control", n))
  q0A = Model$q0A
  q0B = Model$q0B
  is_DE = Model$is_DE
  Pilot.data <- t(sapply(seq_len(ngenes), function(i) sapply(1:(2*n), function(j) rnbinom(1, mu = ifelse(conds[j] == "Case", q0A[i], q0B[i]), size = dispersion))))
  colnames(Pilot.data) <- paste("Sample", 1:(2 * n), sep = " ")
  rownames(Pilot.data) <- paste("Gene", seq_len(ngenes), ifelse(is_DE,"T", "F"), sep = "_")
  return(Pilot.data)
}

Generate.data.Setting.adjust <- function(sample.size, num.repeat, ngenes, par.setting){
  Model = NB.simu.setting.gamma(ngenes, par.setting)

  if(num.repeat==1){
    Data.HIP = Generate.Data.Based.On.Setting(sample.size, Model)
  }
  else{
    # sfExport("num.repeat")
    # sfExport("Generate.Data.Based.On.Setting")
    sfExport("sample.size")
    sfExport("Model")
    Data.HIP = sfLapply(1:num.repeat, function(x)
      Generate.Data.Based.On.Setting(sample.size, Model))
  }
  return(Data.HIP)
}
