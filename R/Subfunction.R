######################################################################################################################
# RNASeqDesign: Sample size and power calculation for RNA-Seq
# Version : 0.1.0
# Authors : Chien-Wei Lin, Ge Liao, Mei-Ling Ting Lee, Yong Seok Park and George C. Tseng
# latest update : 09/27/2018
######################################################################################################################
#Two-beta
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

    Target.Func = sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  mean.p = mean(p.value)
  var.p = var(p.value)
  if(is.null(init.r)) init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p


  lambda = convest(p.value)[[1]]

  try.init.s = lapply(c(1, seq(10, 100, 10)), function(s){
    init = c(lambda, max(0, min(init.r,.9)), s, 1.5, 1)
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

  pars = try.init.s[[i]]$par
  LL = try.init.s[[i]]$value
  return(c(pars, LL))
}

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

    Target.Func = sum(log(lambda*p.value^(r2-1)*(1-p.value)^(s2-1)/beta(r2,s2)+(1-lambda)*p.value^(r-1)*(1-p.value)^(s-1)/beta(r,s)))
    return(Target.Func)
  }


  mean.p = mean(p.value)
  var.p = var(p.value)
  if(is.null(init.r)) init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  if(is.null(lambda)) lambda = max(min(convest(p.value)[[1]], l.upper), 0.7)

  try.init.s = lapply(c(1, seq(10, 100, 10)), function(s){

    init = c(max(0,min(init.r,.9)), s, 1.5, .9)
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

  return(c(lambda, pars, LL))
}

##' Predict EDR based on a pilot data.
##'
##'
##' @title Predict EDR
##' @param Data Input pilot data.
##' @param status A vector of group labels.
##' @param group.name A vector of length two. First element is the group name of first group, and the second element is of second group.
##' @param FDR FDR level.
##' @param M Number of iterations for parametric bootstrap.
##' @param target.N Targeted sample size.
##' @param target.R Targeted sequencing depth.
##' @param target.theta Targeted proportion of sample size in second group compared to first group.
##' @param tol Threshold for deciding if use CDD to estimate the proportion of DE genes.
##' @param tagwiseDisp Use tag-wise dispersion setting or not. Default is F, which uses common dispersion for all genes.
##' @param filter Filter genes based on mean counts?
##' @param filter.level Filter the genes with mean counts less than this level.
##' @return An object of class list is returned:
##' Result is a list of length M (depends on how many times of parametric bootstrap you run), and each component is also a list of length of variable target.R.
##' Each component contains a matrix with target N in row and summary statistics in columns. Summary statistics provided are TP, TN, EDR.
##' @author Chien-Wei Lin
##' @export
Estimate.EDR.from.pilot <- function(Data, status, group.name = c("Control", "Case"), FDR, M,
                                    target.N, target.R = NULL, target.theta = NULL,
                                    tol = 0.1, tagwiseDisp=F,
                                    filter = T, filter.level = 5){
  method = "TwoBeta"
  N0 = sum(status == group.name[1])
  N1 = sum(status == group.name[2])
  theta = N1/N0
  if(is.null(target.theta)) target.theta = theta

  mean.gene = apply(Data, 1, mean)

  Data.filter = if(filter) Data[which(mean.gene>filter.level),] else Data

  R = mean(apply(Data.filter, 2, sum))
  if(is.null(target.R)) target.R = R

  Pilot.data = Data.filter
  colnames(Pilot.data)=1:ncol(Data.filter)

  library(edgeR)
  y <- DGEList(counts = Pilot.data, group = status)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  if(tagwiseDisp==F){
    y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
    delta = rep(1/y$common.dispersion,nrow(Pilot.data))
  } else if (tagwiseDisp==T){
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
    delta = 1/y$tagwise.dispersion
  }
  


  GLM.fit.each <- function(z,delta){
    y0 = z[which(status==group.name[1])]
    y1 = z[which(status==group.name[2])]
    f <- function(z,d) {
      delta_i = d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1 = d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(delta_i+y0)-lgamma(delta_i)-lgamma(y0+1)+y0*log(R/delta_i*exp(beta.0))-(y0+delta_i)*log(1+R/delta_i*exp(beta.0)))-sum(lgamma(delta_i+y1)-lgamma(delta_i)-lgamma(y1+1)+y1*log(R/delta_i*exp(beta.0+beta.1))-(y1+delta_i)*log(1+R/delta_i*exp(beta.0+beta.1)))
      return(Target.Func)
    }

    pt1 <- optim(c(0, 0), f, d = list(delta, R, y0, y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[2])
    var.beta.1 = (1/N0)*((1 + theta*r0)/(theta*R*r1) + (1+theta)/(theta*delta))
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }


  #OUT.pars = t(apply(round(Pilot.data),1,GLM.fit.each))
  OUT.pars = t(sapply(1:nrow(Pilot.data),function(i) GLM.fit.each(round(Pilot.data)[i,],delta=delta[i])))
  colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")
  rownames(OUT.pars)=rownames(Pilot.data)

  model = cbind(OUT.pars[,1:2], delta)

  p.value = 2*pnorm(-abs(OUT.pars[,3]))
  mean.count = rowMeans(Pilot.data)

  q.value = p.adjust(p.value, method = "BH")

  mean.by.group = apply(Pilot.data,1,function(x) tapply(x, status, mean))
  fold.change = (mean.by.group[1, ] + 1)/(mean.by.group[2, ] + 1)

  Fitted.Model = MixtureModel.two.beta(p.value, s.lower = 1, l.upper = 0.9)

  check.if.use.CDD = function(model, uniform.cut = .5){
    r = model[4]
    s = model[5]
    pbeta(1, r, s) - pbeta(uniform.cut, r, s) - uniform.cut
  }

  index.use.CDD = check.if.use.CDD(Fitted.Model)

  if(method == "TwoBeta" & index.use.CDD < tol) {
    print("Use CDD to estimate lambda")
    Fitted.Model = MixtureModel.two.beta.2(p.value, s.lower = 1, l.upper = 0.9)
  }
  p.value.mod <- p.value
  if (any(is.na(p.value.mod))) {
    p.value.mod <- p.value.mod[!is.na(p.value.mod)]
    model = model[!is.na(p.value.mod),]
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

    #check posterior prob
    if(lambda>=0.99){
      return(Fitted.Model)
    }
    else{
      parameter = list(n.old = N0, n.new = target.N,
                       R.old = R, R.new = target.R,
                       theta.old = theta, theta.new = target.theta)
      Resampling <- function(target.N, target.R){
        DE_status_posterior = sapply(1:length(posterior),function(x) sample(c(FALSE,TRUE),1,prob = c(posterior[x],1-posterior[x]),replace=TRUE))

        transform.p.value <- function(p.value.each, delta, DE_status_each, parameter, model, transform.null = F){
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
            statistic.adjust = adjust(n.old, n.new, R.old, R.new, theta.old, theta.new)
            statistic.new = statistic.old*statistic.adjust

            p.matrix = (1-pnorm(abs(statistic.new)))*2
            colnames(p.matrix) = n.new
            rownames(p.matrix) = R.new
            return(p.matrix)
          } else if(transform.null){
            p.matrix = matrix(rep(runif(1),length(n.new)*length(R.new)), nrow = length(R.new))
            colnames(p.matrix) = n.new
            rownames(p.matrix) = R.new
            return(p.matrix)
          }
          else {
            p.matrix = matrix(rep(p.value.each, length(n.new)*length(R.new)), nrow = length(R.new))
            colnames(p.matrix) = n.new
            rownames(p.matrix) = R.new
            return(p.matrix)
          }
        }

        p.value.updated = lapply(1:ngenes,function(x) transform.p.value(p.value.mod[x], delta[x], DE_status_posterior[x], parameter = parameter,
                                                                        model = model[x,], transform.null = F))

        result = lapply(1:length(target.R), function(i){ #for each depth
          p.value.star.posterior = t(sapply(p.value.updated, function(x) x[i,]))
          #for each gene, take out the i-th depth, the length is as n.new
          #target.N x G

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
              Declare_status[which(p.value.star.star <= p.value.cut)] = "DE"
            }
            A = sum((Declare_status=="nonDE")*(!DE_status_posterior))

            B = sum((Declare_status=="nonDE")*(DE_status_posterior))
            D = sum((Declare_status=="DE")*(DE_status_posterior))

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
            EDR_post = D/(B+D)
            Declare_post  = sum(Declare_status=="DE")
            return(c(TP = TP_hat_post, TN = TN_hat_post, EDR = EDR_post))
          }

          Estimate.Posterior.Result = matrix(apply(p.value.star.posterior,2,Estimate_Posterior), ncol = length(target.N)); round(Estimate.Posterior.Result, 3)

          row.names(Estimate.Posterior.Result) = c("TP", "TN", "EDR")
          colnames(Estimate.Posterior.Result) = target.N
          return(Estimate.Posterior.Result)
        })
        return(result)
      }
      Result = lapply(1:M, function(x){
        Resampling(target.N = parameter[[2]], target.R = parameter[[4]])
      })
      return(list(Result = Result,delta=delta))
    }
  }
}
##' Calculating true EDR based on a given data with given true DE information.
##'
##'
##' @title Calculating true EDR
##' @param Data Input pilot data.
##' @param status A vector of group labels.
##' @param group.name A vector of length two. First element is the group name of first group, and the second element is of second group.
##' @param is_DE A vector of True or False which indicate the DE status of each gene.
##' @param FDR FDR level.
##' @param tagwiseDisp Use tag-wise dispersion setting or not. Default is F, which uses common dispersion for all genes.
##' @param filter Filter genes based on mean counts?
##' @param filter.level Filter the genes with mean counts less than this level.
##' @return An object of class list is returned:
##' The object contains a vector of three summary statistics, which are TP_true, TN_true, and EDR_true
##' Each component contains a matrix with target N in row and summary statistics in columns. Summary statistics provided are TP, TN, EDR.
##' @author Chien-Wei Lin
##' @export
Estimate.true.EDR = function(Data, status, group.name = c("Control", "Case"), is_DE,
                             FDR, filter = T,tagwiseDisp=F, filter.level = 5){
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
  if(tagwiseDisp==F){
    y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
    delta = rep(1/y$common.dispersion,nrow(Pilot.data))
  } else if (tagwiseDisp==T){
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
    delta = 1/y$tagwise.dispersion
  }
  
  
  
  GLM.fit.each <- function(z,delta){
    y0 = z[which(status==group.name[1])]
    y1 = z[which(status==group.name[2])]
    f <- function(z,d) {
      delta_i = d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1 = d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(delta_i+y0)-lgamma(delta_i)-lgamma(y0+1)+y0*log(R/delta_i*exp(beta.0))-(y0+delta_i)*log(1+R/delta_i*exp(beta.0)))-sum(lgamma(delta_i+y1)-lgamma(delta_i)-lgamma(y1+1)+y1*log(R/delta_i*exp(beta.0+beta.1))-(y1+delta_i)*log(1+R/delta_i*exp(beta.0+beta.1)))
      return(Target.Func)
    }

    pt1 <- optim(c(0, 0), f, d = list(delta, R, y0, y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[2])
    var.beta.1 = (1/N0)*((1 + theta*r0)/(theta*R*r1) + (1+theta)/(theta*delta))
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }
  
  
  #OUT.pars = t(apply(round(Pilot.data),1,GLM.fit.each))
  OUT.pars = t(sapply(1:nrow(Pilot.data),function(i) GLM.fit.each(round(Pilot.data)[i,],delta[i])))
  colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")
  rownames(OUT.pars)=rownames(Pilot.data)
  
  p.value = 2*pnorm(-abs(OUT.pars[,3]))

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
      Declare_status[which(p.value <= p.value.cut)] = "DE"
    }

    A = sum((Declare_status=="nonDE")*(!is_DE))
    B = sum((Declare_status=="nonDE")*(is_DE))
    C = sum((Declare_status=="DE")*(!is_DE))
    D = sum((Declare_status=="DE")*(is_DE))

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

    return(list(Result = c(TP_true = TP_hat_true, TN_true = TN_hat_true,
                           EDR_true = EDR_true),delta=delta))
  }

}

##' .. content for \description{} (no empty lines) ..
##'
##'
##' @title
##' @param n Sample size.
##' @param r Sequencing depth.
##' @param EDR Predicted EDR.
##' @return A vector of four parameters are returned.
##' @author Chien-Wei Lin
curve.fitting = function(n, r, EDR){
  library(DEoptim)
  lkh<-function(x){
    return(sum((-x[3]*r^(-x[4])-x[1]*n^(-x[2])+1-EDR)^2))
  }
  outOptim <- DEoptim(lkh,lower=c(0, 0, 0, 0),upper=c(100, 100, 10^10, 100), DEoptim.control(trace = F, itermax = 1000))###stochastic fitting
  return(outOptim$optim$bestmem)
}


##' Summarize the result from function Estimate.EDR.from.pilot and output the summarized EDR table and estimated parameters of two-dimensional surface fitting.
##'
##'
##' @title Summary EDR result and surface fitting
##' @param Result Object from the output of function Estimate.EDR.from.pilot.
##' @param target.N Target sample size used to predict EDR.
##' @param target.R Target total sequencing reads used to predict EDR.
##' @return An object of class list is returned:
##' The object contains a matrix of predicted EDR, where the rows are different target N and columns are different target R, and a vector of four parameters from two-dimensional surface fitting.
##' @author Chien-Wei Lin
##' @export
Summarize.RNASeqDesign = function(Result, target.N, target.R){
  Result2 = Result$Result

  M = length(Result2)

  ps = sapply(1:M, function(j){ #each posterior sampling
    sapply(Result2[[j]], function(x) x[3,]) #for each depth
  })
  ps.median = matrix(apply(ps, 1, median), nrow = length(target.N))

  rownames(ps.median) = target.N; colnames(ps.median) = target.R
  predict.N = rep(target.N, length(target.R))
  predict.R = rep(target.R, each = length(target.N))

  surface.para = curve.fitting(n = predict.N, r = predict.R, EDR = ps.median)

  return(list(EDR = ps.median, Surface.para = surface.para))
}

##' Predict EDR at any sample size and sequencing depth after two-dimensional surface fitting.
##'
##'
##' @title Predict EDR from surface fitting
##' @param target.N Target sample size.
##' @param target.R Target sequencing total reads.
##' @param Surface.para Parameters of two dimensional EDR surface, estimated from function Summarize.RNASeqDesign.
##' @return Predicted EDR.
##' @author Chien-Wei Lin
##' @export
Predict.EDR = function(target.N, target.R, Surface.para){
  EDR.curve = function(n, r, x) 1-x[1]*n^(-x[2])-x[3]*r^(-x[4])

  return(EDR.curve(target.N, target.R, Surface.para))
}


Estimate.lambda.from.pilot <- function(Data, status, group.name = c("Control", "Case"), FDR, M,
                                    target.N, target.R = NULL, target.theta = NULL,
                                    tol = 0.1, tagwiseDisp=F,
                                    filter = T, filter.level = 5){
  method = "TwoBeta"
  N0 = sum(status == group.name[1])
  N1 = sum(status == group.name[2])
  theta = N1/N0
  if(is.null(target.theta)) target.theta = theta
  
  mean.gene = apply(Data, 1, mean)
  
  Data.filter = if(filter) Data[which(mean.gene>filter.level),] else Data
  
  R = mean(apply(Data.filter, 2, sum))
  if(is.null(target.R)) target.R = R
  
  Pilot.data = Data.filter
  colnames(Pilot.data)=1:ncol(Data.filter)
  
  library(edgeR)
  y <- DGEList(counts = Pilot.data, group = status)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  if(tagwiseDisp==F){
    y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
    delta = rep(1/y$common.dispersion,nrow(Pilot.data))
  } else if (tagwiseDisp==T){
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
    delta = 1/y$tagwise.dispersion
  }
  
  
  
  GLM.fit.each <- function(z,delta){
    y0 = z[which(status==group.name[1])]
    y1 = z[which(status==group.name[2])]
    f <- function(z,d) {
      delta_i = d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1 = d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(delta_i+y0)-lgamma(delta_i)-lgamma(y0+1)+y0*log(R/delta_i*exp(beta.0))-(y0+delta_i)*log(1+R/delta_i*exp(beta.0)))-sum(lgamma(delta_i+y1)-lgamma(delta_i)-lgamma(y1+1)+y1*log(R/delta_i*exp(beta.0+beta.1))-(y1+delta_i)*log(1+R/delta_i*exp(beta.0+beta.1)))
      return(Target.Func)
    }
    
    pt1 <- optim(c(0, 0), f, d = list(delta, R, y0, y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[2])
    var.beta.1 = (1/N0)*((1 + theta*r0)/(theta*R*r1) + (1+theta)/(theta*delta))
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }
  
  
  #OUT.pars = t(apply(round(Pilot.data),1,GLM.fit.each))
  OUT.pars = t(sapply(1:nrow(Pilot.data),function(i) GLM.fit.each(round(Pilot.data)[i,],delta=delta[i])))
  colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")
  rownames(OUT.pars)=rownames(Pilot.data)
  
  model = cbind(OUT.pars[,1:2], delta)
  
  p.value = 2*pnorm(-abs(OUT.pars[,3]))
  mean.count = rowMeans(Pilot.data)
  
  q.value = p.adjust(p.value, method = "BH")
  
  mean.by.group = apply(Pilot.data,1,function(x) tapply(x, status, mean))
  fold.change = (mean.by.group[1, ] + 1)/(mean.by.group[2, ] + 1)
  
  #Use two method to estimate lamdba
  #mle
  Fitted.Model.MLE = MixtureModel.two.beta(p.value, s.lower = 1, l.upper = 0.99)
  #print("Use CDD to estimate lambda")
  Fitted.Model.CDD = MixtureModel.two.beta.2(p.value, s.lower = 1, l.upper = 0.99)
  
  check.if.use.CDD = function(model, uniform.cut = .5){
    r = model[4]
    s = model[5]
    pbeta(1, r, s) - pbeta(uniform.cut, r, s) - uniform.cut
  }
  
  index.use.CDD = check.if.use.CDD(Fitted.Model.MLE)
  
  index.use.CDD.emp<-sum(p.value>0.5,na.rm=T)/(sum(p.value>0,na.rm=T)*0.85)
  
  return(list(index=index.use.CDD,index.emp=index.use.CDD.emp,MLE=Fitted.Model.MLE,cdd=Fitted.Model.CDD,pvalue=p.value))
  
}

