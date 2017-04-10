if(!"RNASeqPower" %in% list.files(.Library)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("RNASeqPower")
}

if(!"truncnorm" %in% list.files(.Library)){
  install.packages("truncnorm")
}

if(!"ggplot2" %in% list.files(.Library)){
  install.packages("ggplot2")
}

if(!"drc" %in% list.files(.Library)){
  install.packages("drc")
}

Parameter.Estimate <- function(Normalized.HIP.data, R.each){
  n = ncol(Normalized.HIP.data) / 2
  Mean.count = apply(Normalized.HIP.data, 1, mean)
  Mean.count = R.each * Mean.count / mean(Mean.count)

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

  set.seed(12345)
  is_DE <- runif(ngenes) < m1

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
    sfExport("sample.size")
    sfExport("Model")
    Data.HIP = sfLapply(1:num.repeat, function(x)
      Generate.Data.Based.On.Setting(sample.size, Model))
  }
  return(Data.HIP)
}

Subsampling <- function(Data.all, n, conds.all){
  ind = if(n == ncol(Data.all)) 1:ncol(Data.all) else c(sample(which(conds.all%in%"Case"),n,replace=FALSE),sample(which(conds.all%in%"Control"),n,replace=FALSE))
  Data.sub=Data.all[,ind]
  mean.gene=apply(Data.sub,1,mean)

  Data.filter=Data.sub[-which(mean.gene<5),]
  return(Data.filter)
}

Sampling.function <- function(Sample.size, Data, Status, filter.cutoff, sample.ratio = 1){
  Indicator <- list()
  Status.1=unique(Status)[1]
  Status.2=unique(Status)[2]

  for (i in length(Sample.size):1){
    if(i==length(Sample.size)){
      Indicator[[i]]=c(sample(which(Status%in%Status.1),Sample.size[i]*sample.ratio),sample(which(Status%in%Status.2),Sample.size[i]))
    }
    else{
      Indicator[[i]]=c(sample(Indicator[[i+1]][1:(sample.ratio*length(Indicator[[i+1]])/(sample.ratio + 1))], Sample.size[i]*sample.ratio),sample(Indicator[[i+1]][-(1:(sample.ratio*length(Indicator[[i+1]])/(sample.ratio + 1)))],Sample.size[i]))
    }
  }
  tmp.mean=apply(Data[,Indicator[[1]]],1,mean)
  ind.gene=which(tmp.mean<filter.cutoff)
  return(list(Indicator,ind.gene))
}

my.QQ <- function(Data,Data.name){
  R=apply(Data,2,sum)
  mean.gene=apply(Data,1,mean)

  Data.SRA.HIP.filter=Data[mean.gene>5,]
  ## 11310
  mean.by.group=apply(Data.SRA.HIP.filter,1,function(x) tapply(x,colnames(Data.SRA.HIP.filter),mean))
  lfc.SRA=log(mean.by.group[1,]/mean.by.group[2,])


  Data.filter=Data.SRA.HIP.filter
  n = ncol(Data.filter)/2
  conds <- c(rep("Case",n),rep("Control",n))

  ngenes= nrow(Data.filter)
  Pilot.data = Data.filter
  colnames(Pilot.data)=1:24
  library(edgeR)
  y <- DGEList(counts = Pilot.data, group = conds)
  y <- calcNormFactors(y,method="TMM") #Calculate normalization factors to scale the raw library sizes
  normalize.factor=y[[2]][,3]
  nc <- cpm(y, normalized.lib.sizes=FALSE)

  y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  y <- estimateTagwiseDisp(y) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  et <- exactTest(y,dispersion="common")$table$P #Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
  design <- model.matrix(~conds)
  y <- DGEList(counts = Pilot.data, group = conds)
  y <- calcNormFactors(y,method="TMM") #Calculate normalization factors to scale the raw library sizes
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  f <- glmFit(y, design,dispersion=y$common.dispersion)
  lrt <- glmLRT(f, coef=2)$table$PValue


  GLM.fit.each <- function(z){
    y0 = z[which(conds=="Control")]
    y1 = z[which(conds=="Case")]
    f <- function(z,d) {
      theta=d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1=d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(theta+y0)-lgamma(theta)-lgamma(y0+1)+y0*log(R/theta*exp(beta.0))-(y0+theta)*log(1+R/theta*exp(beta.0)))-sum(lgamma(theta+y1)-lgamma(theta)-lgamma(y1+1)+y1*log(R/theta*exp(beta.0+beta.1))-(y1+theta)*log(1+R/theta*exp(beta.0+beta.1)))
      return(Target.Func)
    }

    pt1 <- optim(c(0, 0), f,d=list(theta,R,y0,y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[1])
    var.beta.1 = (1/n)*((1/R)*(1/r0+1/r1)+2*1/theta)
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }
  theta = 1/y$common.dispersion
  Pilot.data.normalized=sapply(1:ncol(Pilot.data),function(x) Pilot.data[,x]/(apply(Pilot.data,2,sum)/mean(apply(Pilot.data,2,sum)))[x])
  R=mean(apply(Pilot.data.normalized,2,sum))
  OUT.pars = t(apply(Pilot.data.normalized,1,GLM.fit.each))
  p.value = 2*pnorm(-abs(OUT.pars[,3]))

  qqplot(et,lrt,main="Exact test vs. Likelihood ratio test")
  abline(0,1,col="red",lwd=3,lty=2)
  qqplot(et,p.value,main="Exact test vs. Wald test")
  abline(0,1,col="red",lwd=3,lty=2)

  return(list(row.names(Pilot.data)[p.adjust(p.value,"BH")<0.05],p.value,et,lrt))
}

Compare.TEST = function(Data){
  ngenes= nrow(Data)
  n = ncol(Data)/2
  conds <- c(rep("Case",n),rep("Control",n))

  mean.gene = apply(Data,1,mean)
  m = Data[mean.gene>5,]

  status=sapply(strsplit(row.names(m),"_"),function(x) x[3])
  status[status=="F"]=FALSE
  status[status=="T"]=TRUE

  R=mean(apply(m, 2,sum))

  y <- DGEList(counts = m, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  y <- estimateTagwiseDisp(y) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  et <- exactTest(y)$table$P #Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
  design <- model.matrix(~conds)
  y <- DGEList(counts = m, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  f <- glmFit(y, design)
  lrt <- glmLRT(f, coef=2)$table$PValue
  theta = 1/y$common.dispersion
  GLM.fit.each <- function(z){
    y0 = z[which(conds=="Control")]
    y1 = z[which(conds=="Case")]
    f <- function(z,d) {
      theta=d[[1]]
      R = d[[2]]
      y0 = d[[3]]
      y1=d[[4]]
      beta.0 = z[1]
      beta.1 = z[2]
      Target.Func = -sum(lgamma(theta+y0)-lgamma(theta)-lgamma(y0+1)+y0*log(R/theta*exp(beta.0))-(y0+theta)*log(1+R/theta*exp(beta.0)))-sum(lgamma(theta+y1)-lgamma(theta)-lgamma(y1+1)+y1*log(R/theta*exp(beta.0+beta.1))-(y1+theta)*log(1+R/theta*exp(beta.0+beta.1)))
      return(Target.Func)
    }
    pt1 <- optim(c(0, 0), f,d=list(theta,R,y0,y1))$par
    r1 = exp(sum(pt1))
    r0 = exp(pt1[1])
    var.beta.1 = (1/n)*((1/R)*(1/r0+1/r1)+2*1/theta)
    statistics = pt1[2]/sqrt(var.beta.1)
    return(c(pt1,statistics))
  }
  OUT.pars = t(apply(m,1,GLM.fit.each))
  p.value = 2*pnorm(-abs(OUT.pars[,3]))
  return(list(cbind(et,lrt,p.value), is_DE = status))
}


power_poisson<-function(n, w=1.0, rho=2.5, mu0=5.0, f, m, m1){

  estimate.power.given.n <- function(r1) {
    alpha_star<-(r1*f)/((m-m1)*(1-f))
    z_alpha<-qnorm(1-alpha_star/2)

    Target.Func = abs(r1 - 1 + pnorm(z_alpha-(rho-1)*sqrt((n*mu0)/(rho/w + 1))))
    return(Target.Func)
  }
  power_w = tryCatch(optimize(estimate.power.given.n, c(0, 1))$minimum, error = function(e) 0)

  return(Wald = power_w)
}

#Poisson
Poisson <- function(x, conds, target.N, w = 1, rho = 1.2, FDR = 0.05,
                    min.mu0 = 5, DE.prop = .1){

  x = x[apply(x, 1, mean) > 5, ]
  mean.group = apply(x, 1, function(y) tapply(y, conds, mean))
  fc = abs(log2((mean.group[1, ]+1)/(mean.group[2, ]+1)))
  mu0 = max(min(round(mean.group[1, order(fc, decreasing = TRUE)[1:round(nrow(x) * DE.prop)]])), min.mu0)

  m<-nrow(x)
  m1<-DE.prop*m

  power=sapply(target.N,function(z) power_poisson(n=z, w=w, rho=rho, mu0=mu0, f=FDR, m=m, m1=m1))

  return(power)
}

#NB
NB.exact <- function(x, conds, target.N, rho = 1.4, DE.prop = 0.1,
                     min.mu0 = 5, alpha = 0.05, FDR = 0.05, w = 1){
  if(!"RnaSeqSampleSize" %in% list.files(.Library)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("RnaSeqSampleSize")
  }
  library("RnaSeqSampleSize")
  require("edgeR")
  x = x[apply(x, 1, mean) > 5, ]
  mean.group = apply(x, 1, function(y) tapply(y, conds, mean))
  fc = abs(log2((mean.group[1, ]+1)/(mean.group[2, ]+1)))
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

  power=sapply(target.N, function(z) est_power(n = z, alpha = alpha, w = w,
                                               f = FDR, m = m, m1 = m1, rho = rho,
                                               lambda0 = mu0, phi0 = vphi_0))
  return(power)
}

#RNASeqpower
RNASeqPower <- function(x, conds, target.N, n.prop = 1, effect, alpha = 0.05){
  if(!"RNASeqPower" %in% list.files(.Library)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("RNASeqPower")
  }
  library(RNASeqPower)
  library(edgeR)
  x=x[apply(m,1,mean)>5,]
  y <- DGEList(counts = x, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y, verbose=TRUE) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  y <- estimateTagwiseDisp(y) #Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
  tagwise.dispersion = y$tagwise.dispersion
  BCV = sqrt(quantile(tagwise.dispersion, probs = 0.5))
  Power = sapply(target.N, function(x) rnapower(sum(x)/(ncol(x)*nrow(x)), n = x, n2 = x*n.prop, cv = BCV, effect = effect, alpha = alpha))
  return(Power)
}


#PROPER
PROPER = function(x, conds, target.N, n.prop, DE.prop = 0.1, effect){
  if(!"PROPER" %in% list.files(.Library)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("PROPER")
  }
  library(edgeR)
  library(PROPER)

  x = x[apply(x, 1, mean) > 5, ] #filter mean count < 5
  mean.group = apply(x, 1, function(y) tapply(y, conds, mean))

  ngenes = nrow(x)

  y <- DGEList(counts = x, group = conds)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateTagwiseDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags

  lOD = log(y$tagwise.dispersion)

  lBaselineExpr = log(rowMeans(x))

  sim.opts = RNAseq.SimOptions.2grp(ngenes = ngenes, p.DE = DE.prop,
                                    lOD = lOD, lBaselineExpr = lBaselineExpr)
  ptm <- proc.time()
  simres = runSims(Nreps = target.N, Nreps2 = target.N*n.prop, sim.opts = sim.opts, DEmethod = "edgeR", nsims = 20)
  ptm2 = proc.time()
  ptm2 - ptm

  powers = comparePower(simres, alpha.type = "fdr", alpha.nominal = 0.05,
                        stratify.by = "expr", target.by = "lfc",
                        strata = c(0, 5, 10, 20, 40, 80, 160, 320, 640, 1280, Inf),
                        delta = effect, filter.by = "expr", strata.filtered = 1) #filter < 5

  return(summaryPower(powers)[, "Marginal power"])
}

DE <- function(Data, status, group.name = c("Control", "Case"),
               method = "GLM", FDR, filter.by.5 = T){
  N0 = sum(status == group.name[1])
  N1 = sum(status == group.name[2])
  theta = N1/N0

  mean.gene=apply(Data,1,mean)

  Data.filter = if(filter.by.5) Data[mean.gene > 5,] else Data

  R=mean(apply(Data.filter,2,sum))

  Pilot.data = Data.filter
  colnames(Pilot.data)=1:ncol(Data.filter)

  library(edgeR)
  y <- DGEList(counts = Pilot.data, group = status)
  y <- calcNormFactors(y) #Calculate normalization factors to scale the raw library sizes
  y <- estimateCommonDisp(y) #Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
  if(method == "edgeR") p.value <- exactTest(y,dispersion="common")$table$P #Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
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

  if(method == "GLM"){
    OUT.pars = t(apply(round(Pilot.data), 1, GLM.fit.each))
    colnames(OUT.pars) = c("Beta0", "Beta1", "Statistics")
    p.value = 2*pnorm(-abs(OUT.pars[,3]))
  }

  q.value = p.adjust(p.value, method = "BH")

  mean.by.group = apply(Pilot.data,1,function(x) tapply(x, status, mean))
  fold.change = (mean.by.group[1, ] + 1)/(mean.by.group[2, ] + 1)
  lfc = abs(log2(fold.change))

  return(list(p = p.value, fold.change = fold.change, dispersion = delta,
              Number.DE = sum(q.value < FDR), method = method))
}


Each.EDR <- function(Result, target.N = NULL, target.N.true = target.N, method = "SeqDesign", True = NULL, variation.True = F, sd.True = NULL, pilot.n = 2, method.mean = "mean", True.upper = NULL, True.lower = NULL, output.MSE = F, Power.max = NULL){

  if(method == "SeqDesign"){
    Result=lapply(Result,function(y) y[[1]])
    EDR=sapply(Result,function(z) apply(sapply(z, function(x) x[[1]][3,]), 1, median))
  } else{
    sample.size = target.N
    EDR = Result
    if(class(EDR) == "list") EDR = do.call("cbind", EDR)
    rownames(EDR) = sample.size
  }

  library(ggplot2)
  library(drc)
  library(DEoptim)

  mean.EDR=apply(EDR,1,function(x) mean(na.omit(x)))
  sd.EDR=apply(EDR,1,function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))

  if(is.null(True)){
    EDR.CI=cbind(mean.EDR,sapply(mean.EDR-1.96*sd.EDR,function(z) max(z,0)),sapply(mean.EDR+1.96*sd.EDR,function(z) min(z,1)),c(as.numeric(names(mean.EDR))+0.5))
    colnames(EDR.CI)=c("EDR","lower","upper","N")
    data=data.frame(cbind(EDR.CI,Type=c(rep("Predicted EDR", nrow(EDR.CI)))))
    pd <- position_dodge(width=0.2)
    data[,1]=as.numeric(as.matrix(data[,1]))
    data[,2]=as.numeric(as.matrix(data[,2]))
    data[,3]=as.numeric(as.matrix(data[,3]))
    data[,4]=as.numeric(as.matrix(data[,4]))
    out=ggplot(data, aes(N,EDR, color=Type)) + ggtitle(paste("N = ", pilot.n, sep = "")) +
      geom_point(aes(shape=Type),size=3, position=pd) +
      scale_color_manual(name = "Method", values=c("red")) +
      scale_shape_manual(name = "Method", values=c(17)) +
      theme_bw() +theme(legend.position="none") +
      geom_errorbar(aes(ymin=lower,ymax=upper),width=1)+ ylim(0, 1)+geom_line(aes(group=Type))+ geom_ribbon(data=data, aes(ymin = lower,ymax=upper), alpha=0.4)
  } else {
    mean.True = if(variation.True) apply(True, 1, function(x) if(method.mean == "mean") mean(na.omit(x)) else if(method.mean == "median") median(na.omit(x))) else True
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
  if(output.MSE) return(list(out, RMSE, RMSE.N)) else return(out)
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

DE.plot = function(N.DE, FDR = 0.05, LFC = 0.4, dataset = "ER", y.min = 500,
                   y.max = 2500, q.full = NULL, lfc.full = NULL,
                   DE.num.FDR.adjust = F){

  num.sig = sapply(1:length(N.DE), function(i){
    x = get(paste(dataset, ".DE.", i, sep = ""))

    num.DE = sapply(x, function(y){
      lfc = log2(y$fold.change)
      q = p.adjust(y$p, method = "BH")
      a = sum(q < FDR & abs(lfc) > LFC)
      if(DE.num.FDR.adjust) a = a * (1-FDR)
      return(a)
    })
  }); num.sig

  if(!is.null(q.full)) DE.full = sum(q.full < FDR & abs(lfc.full) > LFC)

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

Cost.N.R = function(N, R) N*(Cost.N + Cost.R*R)

EDR.curve = function(n, r, x) 1-x[1]*n^(-x[2])-x[3]*r^(-x[4])

cost.R.function = function(N) (Cost/(N)-Cost.N)/Cost.R
power.R.function = function(N) ((Power-1+true.para[1]*N^(-true.para[2]))/(-true.para[3]))^(-1/true.para[4])

find.admissible = function(Cost.discrete, EDR.discrete){
  a = cbind(Cost.discrete, EDR.discrete); index.original = 1:nrow(a)
  b = a[order(a[,1]),]; index.original = index.original[order(a[,1])]

  index.admiss = 1
  s.max = b[1,2]
  for(i in 2:nrow(b)){
    s.next = b[i, 2]
    if(s.next > s.max){
      s.max = s.next
      index.admiss = c(index.admiss, i)
    }
  }
  d = b[index.admiss,]; index.original = index.original[index.admiss]

  return(index.original)
}


T1 = function(Cost.discrete, EDR.discrete, index.original, N.grid, R.grid, R.discrete, Cost){
  plot(Cost.discrete, EDR.discrete, pch = 4, ylim = c(0.75, 0.95), cex = 2,
       xlim = c(min(Cost.discrete), 300000),
       col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "C", ylab = "EDR", main = "T1", axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)); axis(2)
  index.target = as.integer(1:length(index.original) %in% tail(which(Cost.discrete[index.original] <= Cost), 1)) + 1
  points(Cost.discrete[index.original], EDR.discrete[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])
  abline(v = Cost, col = 4); #lines(d[,1], d[,2], col = 2)

  plot(N.grid, R.grid, pch = 4, col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "N", ylab = "R", axes = F, cex = 2,
  ); box()
  axis(1); axis(2, at = R.discrete, cex.axis = 2, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
  points(N.grid[index.original], R.grid[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])

}


T2 = function(Cost.discrete, EDR.discrete, index.original, N.grid, R.grid, R.discrete, EDR){
  plot(Cost.discrete, EDR.discrete, pch = 4, ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 300000),
       col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "C", ylab = "EDR", main = "T2", cex = 2, axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)); axis(2)
  index.target = as.integer(1:length(index.original) %in% which(EDR.discrete[index.original] > EDR)[1]) + 1
  points(Cost.discrete[index.original], EDR.discrete[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])
  abline(h = EDR, col = 4)
  plot(N.grid, R.grid, pch = 4, col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "N", ylab = "R", axes = F, cex = 2); box()
  axis(1); axis(2, at = R.discrete, cex.axis = 2, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
  points(N.grid[index.original], R.grid[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])

}

T3 = function(Cost.discrete.res, EDR.discrete.res, index.original.res, N.grid.res, R.grid.res, R.discrete, EDR){
  plot(Cost.discrete.res, EDR.discrete.res, pch = 1, cex = 2, ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000),
       col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
       xlab = "C", ylab = "EDR", main = "T3", axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)); axis(2)
  index.target.res = as.integer(1:length(index.original.res) %in% which(EDR.discrete.res[index.original.res] > EDR)[1]) + 1
  index.target.res[length(index.target.res)] = 2
  points(Cost.discrete.res[index.original.res], EDR.discrete.res[index.original.res], pch = 1, col = c(1, 2)[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
  text(Cost.discrete.res[index.original.res][index.target.res == 2],
       EDR.discrete.res[index.original.res][index.target.res == 2], labels = c("B", "A"),
       col = 2, cex = 2, lwd = c(1, 2)[index.target.res], pos = 3)
  plot(N.grid.res, R.grid.res, pch = 1, cex = 2, col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
       xlab = "N", ylab = "R", axes = F, xlim = c(50, 200)); box()
  axis(1); axis(2, at = R.discrete, cex.axis = 2, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
  points(N.grid.res[index.original.res], R.grid.res[index.original.res], pch = 1, col = c(1, 2)[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
  text(N.grid.res[index.original.res][index.target.res == 2],
       R.grid.res[index.original.res][index.target.res == 2], labels = c("B", "A"),
       col = 2, cex = 2, lwd = c(1, 2)[index.target.res], pos = 4)
  abline(v = 80, col = 4)
}

T4 = function(Cost.discrete.res, EDR.discrete.res, index.original.res, Cost.discrete.res.130, EDR.discrete.res.130, index.original.res.130,
              N.grid.res, R.grid.res, N.grid.res.130, R.grid.res.130, R.discrete, EDR){
  plot(Cost.discrete.res[-index.original.res], EDR.discrete.res[-index.original.res], pch = 1,
       ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 1000000),
       col = c("gray55"), cex = 2,
       xlab = "C", ylab = "EDR", main = "T4", axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)); axis(2)
  index.target.res = c(rep(1, length(index.original.res)-1), 2)
  points(Cost.discrete.res[index.original.res], EDR.discrete.res[index.original.res], pch = 1, col = c("black", "red")[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
  points(Cost.discrete.res.130[-index.original.res.130], EDR.discrete.res.130[-index.original.res.130], pch = 4, ylim = c(0.75, 0.95),
         col = c("gray55"), cex = 2)
  index.target.res.130 = rep(1, length(index.original.res.130)); index.target.res.130[16] = 2
  points(Cost.discrete.res.130[index.original.res.130], EDR.discrete.res.130[index.original.res.130], pch = 4, col = c("black", "red")[index.target.res.130],
         cex = c(2, 2)[index.target.res.130], lwd = c(1, 2)[index.target.res.130])

  plot(N.grid.res[-index.original.res], R.grid.res[-index.original.res], pch = 1, col = "gray55", cex = 2,
       xlab = "N", ylab = "R", xlim = c(50, 200), ylim = range(R.grid), axes = F); box()
  axis(1); axis(2, at = R.discrete, cex.axis = 2, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
  points(N.grid.res[index.original.res], R.grid.res[index.original.res], pch = 1, col = c("black", "red")[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])

  points(N.grid.res.130[-index.original.res.130], R.grid.res.130[-index.original.res.130], pch = 4, col = "gray55", cex = 2)
  points(N.grid.res.130[index.original.res.130], R.grid.res.130[index.original.res.130], pch = 4, col = c("black", "red")[index.target.res.130],
         cex = c(2, 2)[index.target.res.130], lwd = c(1, 2)[index.target.res.130])
  abline(v = c(80, 130), col = 4)
}

T5 = function(Cost.discrete, EDR.discrete, R.grid, R.discrete){
  plot(Cost.discrete[6:10], EDR.discrete[6:10], pch = 1:5, col = 2, cex = 2, lwd = 2,
       main = "T5", ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000),
       xlab = "C", ylab = "EDR", axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)

  plot(rep(60, 5), unique(R.grid), pch = 1:5, col = "red", cex = 2, lwd = 2,
       xlab = "N", ylab = "R", xlim = c(50, 200), ylim = range(R.grid), axes = F); box()
  axis(1); axis(2, at = R.discrete, cex.axis = 2, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))

}

