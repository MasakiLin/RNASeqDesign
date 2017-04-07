#----------------------------
# Simulation study
#----------------------------
#----------------------------
# 1.1. simulate true data
#----------------------------
load("Parameter setting estimated from HIP data 20M.RData")
load("Normalized.HIP.data.rdata")

par.setting$lfc.mean.sd = c(0, 0.2)

par.setting$Prop.DE = 0.15
par.setting$lfc = c(-log2(1.4), log2(1.4))

num.repeat = 20
sample.size = c(2, 4, 8, 16)
target.N = c(5, 10, 20, 30, 40, 50, 100)
method = "TwoBeta"
library(snowfall)

for(dispersion in c(2, 5, 10, 20)){

  par.setting$dispersion = dispersion

  #wong05
  sfInit(parallel = TRUE, cpus = 20)
  # sfSource("Simulation subfunction.R")

  set.seed(12345)

  for(i in target.N){
    print(i)
    start = Sys.time()
    assign(paste("True.N.", i, sep = ""), Generate.data.Setting.adjust(i, num.repeat, 25000, par.setting))
    do.call(save, list(paste("True.N.", i, sep = ""),
                       file = paste("Simulated.true.20M.20.repeat.dispersion.", round(par.setting$dispersion, 2),
                                    ".lfc.mean.", par.setting$lfc.mean.sd[1], ".sd.",
                                    par.setting$lfc.mean.sd[2], ".N.", i, ".rdata", sep = "")))
    end = Sys.time()
    end - start
  }

  sfStop()

  #calculate true EDR
  sfInit(parallel = TRUE, cpus = 10)
  sfExportAll()

  for(i in 1:length(target.N)){
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("True.", target.N[i], sep = ""),
           sfLapply(1:num.repeat,function(x){
             Data = get(paste("True.N.", target.N[i], sep = ""))[[x]]
             n = ncol(Data)/2
             status = c(rep("Case", n), rep("Control", n))
             Estimate.true.EDR(Data, status, group.name = c("Control", "Case"),
                               FDR = 0.05, filter = T, filter.level = 5, resample.DE = F)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("True.", target.N[i], sep = ""),
                       file = paste("(Simulation)True.EDR.20M.20.repeat.dispersion.", round(par.setting$dispersion, 2),
                                    ".lfc.mean.", par.setting$lfc.mean.sd[1], ".sd.",
                                    par.setting$lfc.mean.sd[2], ".N.", target.N[i],
                                    ".not.considering.lfc.rdata", sep = "")))
  }
  sfStop()

}

#----------------------------
# 1.2. simulate pilot data
#----------------------------
load("Parameter setting estimated from HIP data 20M.RData")
load("Normalized.HIP.data.rdata")

total.reads = "20M"
par.setting$lfc.mean.sd = c(0, 0.2)
par.setting$Prop.DE = 0.15
par.setting$lfc = c(-log2(1.4), log2(1.4))

for(dispersion in c(2, 5, 10, 20)){

  par.setting$dispersion = dispersion #wong05

  #wong05
  library(snowfall)
  sfInit(parallel = TRUE, cpus = 20)
  # sfSource("Simulation subfunction.R")

  num.repeat = 20
  sample.size = c(2, 4, 8, 16)
  target.N = c(5, 10, 20, 30, 40, 50, 100)
  method = "TwoBeta"

  sfExportAll()
  set.seed(12345)

  for(i in sample.size){
    print(i)
    start = Sys.time()
    assign(paste("simulation.N.", i, sep = ""), Generate.data.Setting.adjust(i, num.repeat, 25000, par.setting))
    do.call(save, list(paste("simulation.N.", i, sep = ""),
                       file = paste("Simulated.data.", total.reads, ".20.repeat.dispersion.", round(par.setting$dispersion, 2),
                                    ".lfc.mean.", par.setting$lfc.mean.sd[1], ".sd.",
                                    par.setting$lfc.mean.sd[2], ".N.", i, ".rdata", sep = "")))
    end = Sys.time()
    end - start
  }

  sfExportAll()

  #predict EDR
  for(i in 1:length(sample.size)){ #16
    # for(i in 5){ #16
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("SeqDesign.", sample.size[i], sep = ""),
           sfLapply(1:num.repeat,function(x){
             Data = get(paste("simulation.N.", sample.size[i], sep = ""))[[x]]
             n = ncol(Data)/2
             status = c(rep("Case", n), rep("Control", n))
             Estimate.EDR.from.pilot(Data, status, group.name = c("Control", "Case"),
                                     FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                     target.N = target.N, method = method, s.lower = 1,
                                     resample.DE = T, fc.cut = 1.4, p.cut = 10^-4,
                                     know.DE = T)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("SeqDesign.", sample.size[i], sep = ""),
                       file = paste("(Simulation)SeqDesign.", total.reads, ".20.repeat.dispersion.", round(par.setting$dispersion, 2)
                                    , ".lfc.mean.", par.setting$lfc.mean.sd[1], ".sd.",
                                    par.setting$lfc.mean.sd[2], ".filter.by.5.2B.automate.CDD.pilot.",
                                    sample.size[i], ".predict.EDR.rdata", sep = "")))
  }

}

#----------------------------
# 1.3. Run on other methods
#----------------------------
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.2.rdata")
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.4.rdata")
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.8.rdata")
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.16.rdata")

dispersion = 5

#----------------------------------------------------
Data.2 = simulation.N.2
Data.4 = simulation.N.4
Data.8 = simulation.N.8
Data.16 = simulation.N.16

sample.size = c(2, 4, 8, 16)
target.N = c(5, 10, 20, 30, 40, 50, 100)

library(snowfall)
sfInit(parallel = TRUE, cpus = 20)
sfSource("Simulation subfunction.R")
sfExportAll()

for(i in 1:4){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Poisson.", i, sep = ""),
         sfSapply(1:length(Data.2),function(x){
           Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           Poisson(Data, status, target.n = target.N, rho = 1.4, DE.prop = .15)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Poisson.", i, sep = ""),
                     file = paste("(Simulation)Poisson.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:4){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("RNASeqPower.", i, sep = ""),
         sfSapply(1:length(Data.2),function(x){
           Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           RNASeqPower(Data, status, target.n = target.N, n.prop = 1,
                       effect = 1.4, alpha = 0.0001)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("RNASeqPower.", i, sep = ""),
                     file = paste("(Simulation)RNASeqPower.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
}

if(1){
  #NB, on cluster 08
  for(i in 1:4){ #only up to 20
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("NB.", i, sep = ""),
           sfSapply(1:length(Data.2),function(x){
             Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
             n = ncol(Data)/2
             status = c(rep("Case", n), rep("Control", n))

             NB.exact(Data, status, target.n = target.N, rho = 1.4, DE.prop = .15)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("NB.", i, sep = ""),
                       file = paste("(Simulation)NB.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
  }
}

#Scotty, on window server
for(i in 1:4){ #only up to 20
  for(j in 1:length(Data.2)){ #repeatment
    print(paste("i", i, "j", j))

    Data = get(paste("Data.", sample.size[i], sep = ""))[[j]]
    n = ncol(Data)/2
    status = c(rep("Case", n), rep("Control", n))

    colnames(Data) = c(paste("Case_", 1:n, sep = ""), paste("Control_", 1:n, sep = ""))

    Data = Data[apply(Data, 1, mean) > 5,]

    Data = rbind(c("Gene", colnames(Data)), cbind(row.names(Data), Data))

    write.table(Data, file = paste("Scotty/Scotty.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#PROPER
for(i in 1:2){ #wong05
  # for(i in 3){ #wong05
  for(i in 3:4){ #wong08
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("PROPER.", i, sep = ""),
           sfSapply(1:length(Data.2),function(x){
             Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
             n = ncol(Data)/2
             status = c(rep("Case", n), rep("Control", n))

             PROPER(Data, status, n.prop = 1, target.n = target.N,
                    effect = log(1.4), DE.prop = 0.15)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("PROPER.", i, sep = ""),
                       file = paste("(Simulation)PROPER.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
  }


#----------------------------
# 1.4. summarize results
#----------------------------
sample.size = c(5, 10, 20, 30, 40, 50, 100)
pilot.size = c(2, 4, 8)
# True EDR-----------
for(i in 1:length(sample.size)) load(paste("TrueEDR/(Simulation)True.EDR.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.", sample.size[i], ".not.considering.lfc.rdata", sep = ""))
True = sapply(1:length(sample.size), function(i){
  x = get(paste("True.", sample.size[i], sep = ""))
  x = x[sapply(x, length) != 1]
  mean(sapply(x, function(y) y$Result[3]))
  # length(x)
}); True
plot(sample.size, True)

# SeqDesign-----------
# load("SeqDesign result/(Simulation)SeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.2.predict.EDR.rdata")
# load("SeqDesign result/(Simulation)SeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.4.predict.EDR.rdata")
# load("SeqDesign result/(Simulation)SeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.8.predict.EDR.rdata")
load("SeqDesign result/(Simulation)SeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.2.predict.EDR.not.considering.lfc.rdata")
load("SeqDesign result/(Simulation)SeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.4.predict.EDR.not.considering.lfc.rdata")
load("SeqDesign result/(Simulation)SeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.8.predict.EDR.not.considering.lfc.rdata")

for(i in c(2, 4, 8)) assign(paste("p.", i, sep = ""),
                            Each.EDR(get(paste("SeqDesign.", i, sep = "")), method = "SeqDesign",
                                     True = True, target.N = sample.size,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))
# Poisson------------
load("Other method result/(Simulation)Poisson.dispersion.5.N.2.rdata")
load("Other method result/(Simulation)Poisson.dispersion.5.N.4.rdata")
load("Other method result/(Simulation)Poisson.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.Poisson.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("Poisson.", i, sep = "")), target.N = sample.size,
                                     method = "Poisson", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# RNASeqPower--------
load("Other method result/(Simulation)RNASeqPower.dispersion.5.N.2.rdata")
load("Other method result/(Simulation)RNASeqPower.dispersion.5.N.4.rdata")
load("Other method result/(Simulation)RNASeqPower.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.RNASeqPower.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("RNASeqPower.", i, sep = "")), target.N = sample.size,
                                     method = "RNASeqPower", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# NB------------
load("Other method result/(Simulation)NB.dispersion.5.N.2.rdata")
load("Other method result/(Simulation)NB.dispersion.5.N.4.rdata")
load("Other method result/(Simulation)NB.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.NB.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("NB.", i, sep = "")), target.N = sample.size,
                                     method = "NB", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# PROPER------------
# load("Other method result/(Simulation)PROPER.dispersion.5.N.2.rdata")
# load("Other method result/(Simulation)PROPER.dispersion.5.N.4.rdata")
# load("Other method result/(Simulation)PROPER.dispersion.5.N.8.rdata")
load("Other method result/(Simulation)PROPER.default.dispersion.5.N.2.rdata")
load("Other method result/(Simulation)PROPER.default.dispersion.5.N.4.rdata")
load("Other method result/(Simulation)PROPER.default.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.PROPER.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("PROPER.", i, sep = "")), target.N = sample.size,
                                     method = "PROPER", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))
# Scotty----------
for(i in c(1, 2, 3)){
  tmp=read.table(paste("Other method result/Simulation_Scotty_dispersion_5_N_", pilot.size[i], ".txt", sep = ""), sep = ",")
  assign(paste("Scotty.", i, sep = ""), {
    x = sapply(1:10, function(j) tmp[(1+149*(j-1)):(149*(j-1)+149),][c(sample.size)-1,10]/100)
  })
}

for(i in c(1, 2, 3)) assign(paste("p.Scotty.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("Scotty.", i, sep = "")), target.N = sample.size,
                                     method = "Scotty", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# plot-----------------------------------------------------------------------------------------------
# pdf("Simulation dispersion 5.pdf", width = 18, height = 9)
# pdf("Simulation dispersion 5 SeqDesign not considering lfc PROPER default.pdf", width = 18, height = 9)
pdf("Simulation dispersion 5 SeqDesign not considering lfc PROPER default.pdf", width = 18, height = 9)
multiplot(p.Poisson.2[[1]], p.RNASeqPower.2[[1]], p.NB.2[[1]], p.Scotty.2[[1]], p.PROPER.2[[1]], p.2[[1]],
          p.Poisson.4[[1]], p.RNASeqPower.4[[1]], p.NB.4[[1]], p.Scotty.4[[1]], p.PROPER.4[[1]], p.4[[1]],
          p.Poisson.8[[1]], p.RNASeqPower.8[[1]], p.NB.8[[1]], p.Scotty.8[[1]], p.PROPER.8[[1]], p.8[[1]],

          cols=6,layout=matrix(1:18,ncol=6,byrow=T))

dev.off()

# MSE table
MSE.EDR=matrix(c(p.Poisson.2[[2]],p.Poisson.4[[2]],p.Poisson.8[[2]],
                 p.RNASeqPower.2[[2]],p.RNASeqPower.4[[2]],p.RNASeqPower.8[[2]],
                 p.NB.2[[2]],p.NB.4[[2]],p.NB.8[[2]],
                 p.Scotty.2[[2]],p.Scotty.4[[2]],p.Scotty.8[[2]],
                 p.PROPER.2[[2]],p.PROPER.4[[2]],p.PROPER.8[[2]],
                 p.2[[2]], p.4[[2]], p.8[[2]]),nrow=6,byrow=T)

MSE.N=matrix(c(p.Poisson.2[[3]],p.Poisson.4[[3]],p.Poisson.8[[3]],
               p.RNASeqPower.2[[3]],p.RNASeqPower.4[[3]],p.RNASeqPower.8[[3]],
               p.NB.2[[3]],p.NB.4[[3]],p.NB.8[[3]],
               p.Scotty.2[[3]],p.Scotty.4[[3]],p.Scotty.8[[3]],
               p.PROPER.2[[3]],p.PROPER.4[[3]],p.PROPER.8[[3]],
               p.2[[3]], p.4[[3]], p.8[[3]]),nrow=6,byrow=T)


#---------------------------------------
# 1.5. Wald test vs exact test and LRT
#---------------------------------------


#----------------------------
# 2. Real data applications
#----------------------------
#---------------------------------------
# 2.1. HIV mice data
#---------------------------------------

num.resample = 10

load("HIP.n.2.data.rdata")
load("HIP.n.4.data.rdata")
load("HIP.n.6.data.rdata")
load("HIP.n.8.new50.data.rdata")
load("HIP.n.10.data.rdata")
load("HIP.n.12.data.rdata")

sample.size = c(2, 4, 6, 8, 10, 12)

#RNASeqDesign
method = "TwoBeta"

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfSource("Simulation subfunction.R")
sfExportAll()

for(i in 1:5){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("HIP.n.", sample.size[i], sep = ""),
         sfLapply(1:num.resample,function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))
           Estimate.EDR.from.pilot(Data, status, group.name = c("Control", "Case"),
                                   FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                   target.N = c(n:12,20,30,40,50,100),
                                   method = method, s.lower = 1,
                                   resample.DE = T, fc.cut = 1.4, p.cut = 10^-4)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("HIP.n.", sample.size[i], sep = ""),
                     # file = paste("(HIP)SeqDesign.filter.by.5.2B.automate.decide.if.CDD.pilot.", sample.size[i], ".predict.EDR.rdata", sep = "")))
                     # file = paste("(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.", sample.size[i], ".predict.EDR.rdata", sep = "")))
                     file = paste("(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.", sample.size[i], ".predict.EDR.not.considering.lfc.rdata", sep = "")))
}

i=6
t1=Sys.time()
Data = get(paste("HIP.n.", sample.size[i], sep = ""))
n = ncol(Data)/2
status = c(rep("Case", n), rep("Control", n))
HIP.n.12 = Estimate.EDR.from.pilot(Data, status, group.name = c("Control", "Case"),
                                   FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                   target.N = c(n:12,20,30,40,50,100), method = method, s.lower = 1,
                                   resample.DE = T, fc.cut = 1.4, p.cut = 10^-4)
t2=Sys.time()
print(t2-t1)
do.call(save, list(paste("HIP.n.", sample.size[i], sep = ""),
                   file = paste("(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.", sample.size[i], ".predict.EDR.not.considering.lfc.rdata", sep = "")))

#other methods
load("HIP.n.2.data.rdata")
load("HIP.n.4.data.rdata")
load("HIP.n.6.data.rdata")
load("HIP.n.8.new50.data.rdata")
load("HIP.n.10.data.rdata")
load("HIP.n.12.data.rdata")
HIP.n.12 = list(HIP.n.12)
#----------------------------------------------------
sample.size = c(2, 4, 6, 8, 10, 12)

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfSource("Simulation subfunction.R")
sfExportAll()

#Poisson
for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Poisson.", i, sep = ""),
         sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           Poisson(Data, status, target.n = c(n:12,20,30,40,50,100), rho = 1.4, DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Poisson.", i, sep = ""),
                     file = paste("(HIP)Poisson.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("RNASeqPower.", i, sep = ""),
         sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           RNASeqPower(Data, status, target.n = c(n:12,20,30,40,50,100), n.prop = 1,
                       effect = 1.4, alpha = 0.0001)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("RNASeqPower.", i, sep = ""),
                     file = paste("(HIP)RNASeqPower.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
}

if(1){
  #NB, on cluster 08
  for(i in 1:6){ #only up to 20
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("NB.", i, sep = ""),
           # sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           sapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
             Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
             n = ncol(Data)/2
             status = c(rep("Case", n), rep("Control", n))

             NB.exact(Data, status, target.n = c(n:12,20,30,40,50,100), rho = 1.4, DE.prop = 0.1)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("NB.", i, sep = ""),
                       file = paste("(HIP)NB.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
  }
}

#Scotty, on window server
for(i in 1:6){ #only up to 20
  for(j in 1:length(get(paste("HIP.n.", sample.size[i], sep = "")))){ #repeatment
    print(paste("i", i, "j", j))

    Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[j]]
    n = ncol(Data)/2
    status = c(rep("Case", n), rep("Control", n))

    colnames(Data) = c(paste("Case_", 1:n, sep = ""), paste("Control_", 1:n, sep = ""))

    Data = Data[apply(Data, 1, mean) > 5,]

    Data = rbind(c("Gene", colnames(Data)), cbind(row.names(Data), Data))

    write.table(Data, file = paste("Scotty/Scotty.fold.change.1.4.N.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#PROPER
# for(i in 1:3){ #wong05
# for(i in 3){ #wong05
for(i in 4:5){ #wong08
  i = 5
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("PROPER.", i, sep = ""),
         sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           PROPER(Data, status, n.prop = 1, target.n = c(n:12,20,30,40,50,100),
                  effect = log(1.4), DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("PROPER.", i, sep = ""),
                     file = paste("(HIP)PROPER.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
}

i=6
print(i)
t1=Sys.time()

assign(paste("PROPER.", i, sep = ""),
       sapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
         Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
         n = ncol(Data)/2
         status = c(rep("Case", n), rep("Control", n))

         PROPER(Data, status, n.prop = 1, target.n = c(n:12,20,30,40,50,100),
                effect = log(1.4), DE.prop = 0.1)
       }))
t2=Sys.time()
print(t2-t1)
do.call(save, list(paste("PROPER.", i, sep = ""),
                   file = paste("(HIP)PROPER.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))

#---------------------------------------
# 2.2. TCGA ER positive vs netative
#---------------------------------------
#DE
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))
source("Real data analysis TCGA subfunction.R")

repeat.num = 30
# sample.size = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 105, 110, 115, 120)
sample.size = c(2:20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130:150)
set.seed(12345)
Indicator.ER.pred = lapply(1:repeat.num,function(x) Sampling.function.2(sample.size, Data = TCGA.ER$Data, Status = TCGA.ER$status, 10, sample.ratio = 3))

for(i in 1:length(sample.size)) assign(paste("Indicator.ER.pred.", i, sep = ""),
                                       lapply(Indicator.ER.pred, function(x) list(x[[1]][[i]],x[[2]])))

a = paste("Indicator.ER.pred.", 1:length(sample.size), sep = "")
do.call(save, list(list = a, file = "TCGA.ER.indicator.N2-20.30.40.50.60.70.80.90.100.120.130-150.RData"))

# run-------------------------------------------
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))
load("TCGA.ER.indicator.N2-20.30.40.50.60.70.80.90.100.120.130-150.RData")
repeat.num = 30
sample.size = c(2:20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130:150)

library(snowfall)
sfInit(parallel = TRUE, cpus = 30)
sfSource("Real data analysis TCGA subfunction.R")

sfExportAll()

for(i in 1:length(sample.size)){
  print(i); sfExport("i")
  t1=Sys.time()
  assign(paste("ER.DE.", i, sep = ""),
         sfLapply(1:repeat.num,function(x){
           Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
           Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.ER$status[Ind[[1]]]
           DE(Data, status, group.name = c("Negative", "Positive"), method = "GLM", FDR = 0.05)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("ER.DE.", i, sep = ""), file = paste("ER/DE/TCGA.ER.N", sample.size[i], ".DE.rdata", sep = "")))
}

#RNASeqDesign
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))

if(0){
  source("Real data analysis TCGA subfunction.R")

  sample.size = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

  set.seed(12345)
  Indicator.ER.pred = lapply(1:10,function(x) Sampling.function.2(sample.size, Data = TCGA.ER$Data, Status = TCGA.ER$status, filter.cutoff = 10, sample.ratio = 3))

  for(i in 1:length(sample.size)) assign(paste("Indicator.ER.pred.", i, sep = ""),
                                         lapply(Indicator.ER.pred, function(x) list(x[[1]][[i]],x[[2]])))

  save(Indicator.ER.pred.1, Indicator.ER.pred.2, Indicator.ER.pred.3,
       Indicator.ER.pred.4, Indicator.ER.pred.5, Indicator.ER.pred.6,
       Indicator.ER.pred.7, Indicator.ER.pred.8, Indicator.ER.pred.9,
       Indicator.ER.pred.10, Indicator.ER.pred.11, Indicator.ER.pred.12,
       Indicator.ER.pred.13, Indicator.ER.pred.14,
       file = "TCGA.ER.indicator.UnBalance.3to1.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

}

load("TCGA.ER.indicator.UnBalance.3to1.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

sample.size = target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
method = "TwoBeta"

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
# sfSource("Real data analysis TCGA subfunction.R")
sfSource("Simulation subfunction.R")
sfExportAll()
# for(i in 1:length(sample.size)) sfExport(paste("Indicator.ER.pred.", i, sep = ""))

for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("ER.Predict.", i, ".2B", sep = ""),
         sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
           Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.ER$status[Ind[[1]]]; print(status)
           Estimate.EDR.from.pilot(Data, status, group.name = c("Negative", "Positive"),
                                   FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                   target.N = target.N, method = method, s.lower = 1,
                                   resample.DE = T, fc.cut = 1.4, p.cut = 10^-4)
         }))
  t2=Sys.time()
  print(t2-t1)
  # do.call(save, list(paste("ER.Predict.", i, ".2B", sep = ""), file = paste("ER/TCGA.ER.N", target.N[i], ".predict.EDR.2B.automate.CDD.rdata", sep = "")))
  # do.call(save, list(paste("ER.Predict.", i, ".2B", sep = ""), file = paste("ER/TCGA.ER.N", target.N[i], ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.rdata", sep = "")))
  do.call(save, list(paste("ER.Predict.", i, ".2B", sep = ""), file = paste("ER/TCGA.ER.N", target.N[i], ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.not.considering.lfc.rdata", sep = "")))
}

#Other methods
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))

load("TCGA.ER.indicator.UnBalance.3to1.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

sample.size = target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfSource("Simulation subfunction.R")
sfExportAll()

for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("ER.Poisson.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
           Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.ER$status[Ind[[1]]]

           Poisson(Data, status, target.n = target.N, rho = 1.4, DE.prop = .3)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("ER.Poisson.", i, sep = ""), file = paste("ER/TCGA.ER.Poisson.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("ER.RNASeqPower.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
           Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.ER$status[Ind[[1]]]

           RNASeqPower(Data, status, target.N, n.prop = 3, effect = 1.4, alpha = 0.0001)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("ER.RNASeqPower.", i, sep = ""), file = paste("ER/TCGA.ER.RNASeqPower.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#NB, on cluster 08
if(0){
  for(i in 1:6){ #only up to 20
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("ER.NB.", i, sep = ""),
           sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
             Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
             Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
             status = TCGA.ER$status[Ind[[1]]]

             NB.exact(Data, status, target.n = target.N, rho = 1.4, DE.prop = 0.3)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("ER.NB.", i, sep = ""), file = paste("ER/TCGA.ER.NB.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
  }
}
#PROPER
for(i in 1:6){ #only up to 20
  # for(i in 3:4){ #only up to 20
  for(i in 5:6){ #only up to 20
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("ER.PROPER.", i, sep = ""),
           sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
             Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
             Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
             status = TCGA.ER$status[Ind[[1]]]

             PROPER(Data, status, n.prop = 3, target.n = target.N, effect = log(1.4), DE.prop = 0.3)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("ER.PROPER.", i, sep = ""), file = paste("ER/TCGA.ER.PROPER.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
  }

  #Scotty, on window server
  for(i in 1:6){ #only up to 20
    print(paste("i", i))
    for(j in 1:10){ #repeatment
      print(paste("j", j))
      Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[j]]
      Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
      status = TCGA.ER$status[Ind[[1]]]
      n = sample.size[i]

      colnames(Data) = c(paste("Early_", 1:n, sep = ""), paste("Late_", 1:n, sep = ""))

      Data = Data[apply(Data, 1, mean) > 5,]

      Data=rbind(c("Gene", colnames(Data)), cbind(row.names(Data), Data))

      write.table(Data, file = paste("Scotty/ER.Data.", sample.size[i], ".", j, ".txt", sep = ""),
                  row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    }
  }




















#---------------------------------------
# 2.3. TCGA early stage vs late stage
#---------------------------------------
#DE
TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))
source("Real data analysis TCGA subfunction.R")

repeat.num = 30
# sample.size = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 105, 110, 115, 120)
sample.size = c(2:20, 30, 40, 50, 60, 70, 80:120)
set.seed(12345)
Indicator.Stage.pred = lapply(1:repeat.num,function(x) Sampling.function.2(sample.size, Data = TCGA.Stage$Data, Status = TCGA.Stage$status, 10))

for(i in 1:length(sample.size)) assign(paste("Indicator.Stage.pred.", i, sep = ""),
                                       lapply(Indicator.Stage.pred, function(x) list(x[[1]][[i]],x[[2]])))

a = paste("Indicator.Stage.pred.", 1:length(sample.size), sep = "")
# do.call(save, list(list = a, file = "TCGA.Stage.indicator.N10.20.30.40.50.60.70.80.90.100-120.RData"))
# do.call(save, list(list = a, file = "TCGA.Stage.indicator.N10.20.30.40.50.60.70.80.90.100.105.110.115.120.RData"))
# do.call(save, list(list = a, file = "TCGA.Stage.indicator.N10.20.30.40.50.60.70.80-120.RData"))
do.call(save, list(list = a, file = "TCGA.Stage.indicator.N2-20.30.40.50.60.70.80-120.RData"))

# run-------------------------------------------
TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))
# load("TCGA.Stage.indicator.N10.20.30.40.50.60.70.80.90.100.105.110.115.120.RData")
# load("TCGA.Stage.indicator.N10.20.30.40.50.60.70.80-120.RData")
load("TCGA.Stage.indicator.N2-20.30.40.50.60.70.80-120.RData")
repeat.num = 30
# sample.size = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 105, 110, 115, 120)
# sample.size = c(10, 20, 30, 40, 50, 60, 70, 80:120)
sample.size = c(2:20, 30, 40, 50, 60, 70, 80:120)

library(snowfall)
sfInit(parallel = TRUE, cpus = 30)
sfSource("Real data analysis TCGA subfunction.R")

sfExportAll()

for(i in 1:length(sample.size)){
  # result = sfLapply(1:30, function(i){
  print(i); sfExport("i")
  t1=Sys.time()
  assign(paste("Stage.DE.", i, sep = ""),
         sfLapply(1:repeat.num,function(x){
           # lapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]
           # DE(Data, status, group.name = c("Early", "Late"), method = "edgeR", FDR = 0.05)
           DE(Data, status, group.name = c("Early", "Late"), method = "GLM", FDR = 0.05)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.DE.", i, sep = ""), file = paste("Stage/DE/TCGA.Stage.N", sample.size[i], ".DE.rdata", sep = "")))
  # do.call(save, list(paste("Stage.DE.", i, sep = ""), file = paste("Stage/DE/TCGA.Stage.N", sample.size[i], ".DE.edgeR.rdata", sep = "")))
  # return(t2-t1)
  # })
}

#RNASeqDesign
TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))
source("Real data analysis TCGA subfunction.R")

if(0){
  set.seed(12345)
  Indicator.Stage.pred = lapply(1:10,function(x) Sampling.function.2(c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), Data = TCGA.Stage$Data, Status = TCGA.Stage$status, 10))

  for(i in 1:length(sample.size)) assign(paste("Indicator.Stage.pred.", i, sep = ""),
                                         lapply(Indicator.Stage.pred, function(x) list(x[[1]][[i]],x[[2]])))

  save(Indicator.Stage.pred.1, Indicator.Stage.pred.2, Indicator.Stage.pred.3,
       Indicator.Stage.pred.4, Indicator.Stage.pred.5, Indicator.Stage.pred.6,
       Indicator.Stage.pred.7, Indicator.Stage.pred.8, Indicator.Stage.pred.9,
       Indicator.Stage.pred.10, Indicator.Stage.pred.11, Indicator.Stage.pred.12,
       Indicator.Stage.pred.13, Indicator.Stage.pred.14,
       file = "TCGA.Stage.indicator.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

  TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))
}

load("TCGA.Stage.indicator.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

sample.size = target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
method = "TwoBeta"

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
# sfSource("Real data analysis TCGA subfunction.R")
sfSource("Simulation subfunction.R")

sfExportAll()

for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.Predict.", i, ".2B", sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]
           Estimate.EDR.from.pilot(Data, status, group.name = c("Early", "Late"),
                                   FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                   target.N = target.N, method = method, s.lower = 1,
                                   resample.DE = T, fc.cut = 1.4, p.cut = 10^-4)
         }))
  t2=Sys.time()
  print(t2-t1)
  # do.call(save, list(paste("Stage.Predict.", i, ".2B", sep = ""), file = paste("Stage/TCGA.Stage.N", sample.size[i], ".predict.EDR.2B.rdata", sep = "")))
  # do.call(save, list(paste("Stage.Predict.", i, ".2B", sep = ""), file = paste("Stage/TCGA.Stage.N", sample.size[i], ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.rdata", sep = "")))
  do.call(save, list(paste("Stage.Predict.", i, ".2B", sep = ""), file = paste("Stage/TCGA.Stage.N", sample.size[i], ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.not.considering.lfc.rdata", sep = "")))
}


#Other methods
TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))

load("TCGA.Stage.indicator.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

sample.size = target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfSource("Simulation subfunction.R")
sfExportAll()
# for(i in 1:length(sample.size)) sfExport(paste("Indicator.Stage.pred.", i, sep = ""))

for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.Poisson.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]

           Poisson(Data, status, target.n = target.N, rho = 1.4, DE.prop = .1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.Poisson.", i, sep = ""), file = paste("Stage/TCGA.Stage.Poisson.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.RNASeqPower.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]

           RNASeqPower(Data, status, target.N, n.prop = 1, effect = 1.4, alpha = 0.0001)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.RNASeqPower.", i, sep = ""), file = paste("Stage/TCGA.Stage.RNASeqPower.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#NB
if(0){
  for(i in 1:6){ #only up to 20
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("Stage.NB.", i, sep = ""),
           sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
             Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
             Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
             status = TCGA.Stage$status[Ind[[1]]]

             NB.exact(Data, status, target.n = target.N, rho = 1.4, DE.prop = 0.1)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("Stage.NB.", i, sep = ""), file = paste("Stage/TCGA.Stage.NB.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
  }
}
#PROPER
for(i in 1:6){ #only up to 20
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.PROPER.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]

           PROPER(Data, status, n.prop = 1, target.n = target.N, effect = log(1.4), DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.PROPER.", i, sep = ""), file = paste("Stage/TCGA.Stage.PROPER.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#Scotty
for(i in 1:6){ #only up to 20
  print(paste("i", i))
  for(j in 1:10){ #repeatment
    print(paste("j", j))
    Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[j]]
    Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
    status = TCGA.Stage$status[Ind[[1]]]
    n = sample.size[i]

    colnames(Data) = c(paste("Early_", 1:n, sep = ""), paste("Late_", 1:n, sep = ""))

    Data = Data[apply(Data, 1, mean) > 5,]

    Data=rbind(c("Gene", colnames(Data)), cbind(row.names(Data), Data))

    write.table(Data, file = paste("Scotty/Stage.Data.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#-----------------------------------------------
# 2.4. Summarize real data application results
#-----------------------------------------------
source("../Simulation/Simulation subfunction.R")
# HIP-----------------
sample.size = c(2, 4, 6, 8, 10, 20) #for predicted EDR
target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
# SeqDesign-----------
# load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.2.predict.EDR.rdata")
# load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.4.predict.EDR.rdata")
# load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.10.predict.EDR.rdata")
# load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.12.predict.EDR.rdata")
load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.2.predict.EDR.not.considering.lfc.rdata")
load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.4.predict.EDR.not.considering.lfc.rdata")
load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.10.predict.EDR.not.considering.lfc.rdata")
load("Mouse HIV/SeqDesign result/(HIP)SeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.12.predict.EDR.not.considering.lfc.rdata")
True.SeqDEsign=apply(sapply(HIP.n.12[[1]], function(x) x[3,]), 1, median)

for(i in c(2, 4, 10)) assign(paste("p.", i, sep = ""),
                             Each.EDR(get(paste("HIP.n.", i, sep = "")), method = "SeqDesign",
                                      target.N = c(i:12,20,30,40,50,100),
                                      target.N.true = c(12, 20, 30, 40, 50, 100),
                                      True = True.SeqDEsign,
                                      variation.True = F, sd.True = NULL))
# Poisson------------
load("Mouse HIV/Other method results/(HIP)Poisson.fold.change.1.4.N.2.rdata")
load("Mouse HIV/Other method results/(HIP)Poisson.fold.change.1.4.N.4.rdata")
load("Mouse HIV/Other method results/(HIP)Poisson.fold.change.1.4.N.10.rdata")
load("Mouse HIV/Other method results/(HIP)Poisson.fold.change.1.4.N.12.rdata")
for(i in c(1, 2, 5)) assign(paste("p.Poisson.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("Poisson.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     method = "Poisson", True = Poisson.6, target.N.true = c(12, 20, 30, 40, 50, 100),
                                     variation.True = F, sd.True = NULL))

# RNASeqPower--------
load("Mouse HIV/Other method results/(HIP)RNASeqPower.fold.change.1.4.N.2.rdata")
load("Mouse HIV/Other method results/(HIP)RNASeqPower.fold.change.1.4.N.4.rdata")
load("Mouse HIV/Other method results/(HIP)RNASeqPower.fold.change.1.4.N.10.rdata")
load("Mouse HIV/Other method results/(HIP)RNASeqPower.fold.change.1.4.N.12.rdata")

for(i in c(1, 2, 5)) assign(paste("p.RNASeqPower.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("RNASeqPower.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "RNASeqPower", True = RNASeqPower.6,
                                     variation.True = F, sd.True = NULL))

# NB------------
load("Mouse HIV/Other method results/(HIP)NB.fold.change.1.4.N.2.rdata")
load("Mouse HIV/Other method results/(HIP)NB.fold.change.1.4.N.4.rdata")
load("Mouse HIV/Other method results/(HIP)NB.fold.change.1.4.N.10.rdata")
load("Mouse HIV/Other method results/(HIP)NB.fold.change.1.4.N.12.rdata")
for(i in c(1, 2, 5)) assign(paste("p.NB.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("NB.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "NB", True = NB.6,
                                     variation.True = F, sd.True = NULL))

# PROPER------------
# load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.2.rdata")
# load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.4.rdata")
# load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.10.rdata")
# load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.12.rdata")
load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.2.default.rdata")
load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.4.default.rdata")
load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.10.default.rdata")
load("Mouse HIV/Other method results/(HIP)PROPER.fold.change.1.4.N.12.default.rdata")
for(i in c(1, 2, 5)) assign(paste("p.PROPER.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("PROPER.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "PROPER", True = PROPER.6,
                                     variation.True = F, sd.True = NULL))
# Scotty----------
for(i in c(1, 2, 5)){
  tmp=read.table(paste("Mouse HIV/Scotty result/scotty_fold_change_1.4_result_n_", sample.size[i], ".txt", sep = ""), sep = ",")
  assign(paste("Scotty.", i, sep = ""), {
    x = sapply(1:10, function(j) tmp[(1+149*(j-1)):(149*(j-1)+149),][c(sample.size[i]:12,20,30,40,50,100)-1,10]/100)
  })
}
tmp=read.table(paste("Mouse HIV/Scotty result/scotty_fold_change_1.4_result_n_12.txt", sep = ""), sep = ",")
Scotty.6 = tmp[c(12,20,30,40,50,100)-1,10]/100

for(i in c(1, 2, 5)) assign(paste("p.Scotty.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("Scotty.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "Scotty", True = Scotty.6,
                                     variation.True = F, sd.True = NULL))


# TCGA ER-----------------
# load("TCGA/ER.True.EDR.6000.cut.RData") #True.EDR, True.CI
load("TCGA/ER.True.EDR.RData") #True.EDR, True.CI
target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
target.N.true = c(20, 30, 40, 50, 60, 70, 80, 90, 100)
# SeqDesign-----
# for(i in sample.size[1:6]) load(paste("TCGA/ER/TCGA.ER.N", i, ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.rdata", sep = ""))
for(i in sample.size[1:6]) load(paste("TCGA/ER/TCGA.ER.N", i, ".predict.EDR.2B.rdata", sep = ""))

for(i in c(2, 5, 6)) assign(paste("p.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.Predict.", i, ".2B", sep = "")),
                                     method = "SeqDesign", target.N = target.N,
                                     target.N.true = target.N.true,
                                     pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Poisson-------
for(i in sample.size[1:6]) load(paste("TCGA/ER/Result other methods/TCGA.ER.Poisson.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.Poisson.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.Poisson.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "Poisson", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# RNASeqPower-------
for(i in sample.size[1:6]) load(paste("TCGA/ER/Result other methods/TCGA.ER.RNASeqPower.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.RNASeqPower.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.RNASeqPower.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "RNASeqPower", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# NB-------
for(i in sample.size[1:6]) load(paste("TCGA/ER/Result other methods/TCGA.ER.NB.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.NB.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.NB.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "NB", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))
# PROPER-------
# for(i in sample.size[1:6]) load(paste("TCGA/ER/Result other methods/TCGA.ER.PROPER.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in sample.size[1:6]) load(paste("TCGA/ER/Result other methods/TCGA.ER.PROPER.fold.change.1.4.N", i, ".default.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.PROPER.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.PROPER.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "PROPER", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Scotty------
for(i in 1:6){
  tmp=read.table(paste("TCGA/Scotty result/er_result_fc_1.4_n_", sample.size[i], ".txt", sep = ""),sep = ",")
  assign(paste("ER.Scotty.", i, sep = ""), {
    x = sapply(1:10, function(j) tmp[(1+149*(j-1)):(149*(j-1)+149),10]/100)
    x = x[target.N-1,]
  })
}

for(i in c(2, 5, 6)) assign(paste("p.Scotty.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.Scotty.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "Scotty", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# TCGA Stage-----------------
# load("TCGA/Stage.True.EDR.1400.cut.RData") #True.EDR, True.CI
load("TCGA/Stage.True.EDR.RData") #True.EDR, True.CI
# SeqDesign-----
# for(i in sample.size[1:6]) load(paste("TCGA/Stage/TCGA.Stage.N", i, ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.rdata", sep = ""))
for(i in sample.size[1:6]) load(paste("TCGA/Stage/TCGA.Stage.N", i, ".predict.EDR.2B.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.Predict.", i, ".2B", sep = "")),
                                     method = "SeqDesign",
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Poisson-------
for(i in sample.size[1:6]) load(paste("TCGA/Stage/Result other methods/TCGA.Stage.Poisson.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.Poisson.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.Poisson.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "Poisson", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# RNASeqPower-------
for(i in sample.size[1:6]) load(paste("TCGA/Stage/Result other methods/TCGA.Stage.RNASeqPower.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.RNASeqPower.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.RNASeqPower.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "RNASeqPower", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# NB-------
for(i in sample.size[1:6]) load(paste("TCGA/Stage/Result other methods/TCGA.Stage.NB.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.NB.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.NB.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "NB", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))
# PROPER-------
# for(i in sample.size[1:6]) load(paste("TCGA/Stage/Result other methods/TCGA.Stage.PROPER.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in sample.size[1:6]) load(paste("TCGA/Stage/Result other methods/TCGA.Stage.PROPER.fold.change.1.4.N", i, ".default.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.PROPER.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.PROPER.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "PROPER", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Scotty------
for(i in 1:6){
  tmp=read.table(paste("TCGA/Scotty result/stage_result_fc_1.4_n_", sample.size[i], ".txt", sep = ""),sep = ",")
  assign(paste("Stage.Scotty.", i, sep = ""), {
    x = sapply(1:10, function(j) tmp[(1+149*(j-1)):(149*(j-1)+149),10]/100)
    x = x[target.N-1,]
  })
}

for(i in c(2, 5, 6)) assign(paste("p.Scotty.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.Scotty.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "Scotty", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# plot-----------------------------------------------------------------------------------------------
# pdf("RealDataEDR new cut.pdf", width = 18, height = 27)
pdf("RealDataEDR PROPER default SeqDesign not considering lfc.pdf", width = 18, height = 27)
multiplot(p.Poisson.2, p.RNASeqPower.2, p.NB.2, p.Scotty.2, p.PROPER.2, p.2,
          p.Poisson.4, p.RNASeqPower.4, p.NB.4, p.Scotty.4, p.PROPER.4, p.4,
          p.Poisson.10, p.RNASeqPower.10, p.NB.10, p.Scotty.10, p.PROPER.10, p.10,

          p.Poisson.4.ER, p.RNASeqPower.4.ER, p.NB.4.ER, p.Scotty.4.ER, p.PROPER.4.ER, p.4.ER,
          p.Poisson.10.ER, p.RNASeqPower.10.ER, p.NB.10.ER, p.Scotty.10.ER, p.PROPER.10.ER, p.10.ER,
          p.Poisson.20.ER, p.RNASeqPower.20.ER, p.NB.20.ER, p.Scotty.20.ER, p.PROPER.20.ER, p.20.ER,

          p.Poisson.4.Stage, p.RNASeqPower.4.Stage, p.NB.4.Stage, p.Scotty.4.Stage, p.PROPER.4.Stage, p.4.Stage,
          p.Poisson.10.Stage, p.RNASeqPower.10.Stage, p.NB.10.Stage, p.Scotty.10.Stage, p.PROPER.10.Stage, p.10.Stage,
          p.Poisson.20.Stage, p.RNASeqPower.20.Stage, p.NB.20.Stage, p.Scotty.20.Stage, p.PROPER.20.Stage, p.20.Stage,
          cols=6,layout=matrix(1:54,ncol=6,byrow=T))

dev.off()


#---------------------------------------
# 3. Illustration of T1-T5
#---------------------------------------
# parameter ---------------------------------------------------------------
# samplesize.predict = c(5,10,20,30,40,50,100,150)
samplesize.predict = c(10,20,30,40,50,80,100,120,150)
read = seq(10^6, 6*10^6, 10^6)
# read = seq(10^6, 3*10^6, length = 6)
# read = seq(10^5, 10^6, length = 6)
# read.predict = seq(10^6, 12*10^6, 10^6)
# read.predict = c(10*10^6, 12*10^6, 15*10^6, 20*10^6, 30*10^6, 60*10^6)
read.predict = seq(10^7, 10^8, 10^7)
# n.repeat = 20
n.repeat = 20

Cost.N=500*2 #cost per sample
Cost.R=(25/10^6)*2 #cost per read
# Cost=80000 #total budget
# cost.R.function = function(N) (Cost/(N-1)-Cost.N)/Cost.R
cost.R.function = function(N) (Cost/(N)-Cost.N)/Cost.R
# Power = 0.8
power.R.function = function(N) ((Power-1+true.para[1]*N^(-true.para[2]))/(-true.para[3]))^(-1/true.para[4])
#========================================
# if(0){
# curve.fitting = function(n, r, EDR){
#   lkh<-function(x){
#     return(sum((-x[4]*r^(-x[5])-x[2]*n^(-x[3])+x[1]-EDR)^2))
#   }
#   outOptim<-DEoptim(lkh,lower=c(0, 0, 0, 0, 0),upper=c(5, 5, 5, 5, 5), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
#   return(outOptim$optim$bestmem)
# }
#
# EDR.curve = function(n, r, x) x[1]-x[2]*n^(-x[3])-x[4]*r^(-x[5])
#
# Optim.Function <- function(a,b,Cost, EDR.curve.para){
#   fBUM <- function(z) {
#     N=z
#     Target.Func=-(EDR.curve.para[1]-EDR.curve.para[2]*N^(-EDR.curve.para[3])-EDR.curve.para[4]*((Cost/N-a)/b)^(-EDR.curve.para[5]))
#     return(Target.Func)
#   }
#   pars = optim(6, fBUM,method= "L-BFGS-B",upper=Cost/a,lower=5)$par
#   N.opt=round(pars)
#   R=(Cost/N.opt-a)/b
#   Target.Func=(EDR.curve.para[1]-EDR.curve.para[2]*N.opt^(-EDR.curve.para[3])-EDR.curve.para[4]*((Cost/N.opt-a)/b)^(-EDR.curve.para[5]))
#   Cost.opt=N.opt*(a+b*R)
#
#   return(c(N.star=N.opt,R.star=R,EDR.star=Target.Func,Cost.star=Cost.opt))
# }
# }

#for discrete plot

N.discrete = seq(50, 200, 10)
R.discrete = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid = rep(N.discrete, each = length(R.discrete))
R.grid = rep(R.discrete, times = length(N.discrete))

# par(mfrow = c(1, 2))
# plot(N.grid, R.grid, pch = 3, col = "gray")
Cost.discrete = apply(cbind(N.grid, R.grid), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete = apply(cbind(N.grid, R.grid), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete[EDR.discrete < 0] = 0

# plot(Cost.discrete, EDR.discrete, pch = 3)

if(0){
  index = seq(1, 546, 6)
  plot(Cost.discrete[index], EDR.discrete[index], pch = 3)
  points(Cost.discrete[index+1], EDR.discrete[index+1], pch = 2)
}

#find admissible points
a = cbind(Cost.discrete, EDR.discrete); index.original = 1:nrow(a)
b = a[order(a[,1]),]; index.original = index.original[order(a[,1])]
# b = tapply(a[,2], a[,1], function(x) max(x))
# b = cbind(as.numeric(names(b)), b)

# b[b[,2] < 0,2] = 0
# index.admiss = diff(c(-1, b[,2], 1)) > 0; index.admiss = c(index.admiss[-length(index.admiss)])

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

#Restricted analysis----N=100
N.discrete.res = seq(50, 80, 10)
R.discrete.res = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid.res = rep(N.discrete.res, each = length(R.discrete.res))
R.grid.res = rep(R.discrete.res, times = length(N.discrete.res))

# par(mfrow = c(1, 2))
# plot(N.grid, R.grid, pch = 3, col = "gray")
Cost.discrete.res = apply(cbind(N.grid.res, R.grid.res), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res = apply(cbind(N.grid.res, R.grid.res), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res[EDR.discrete.res < 0] = 0

#find admissible points
a.res = cbind(Cost.discrete.res, EDR.discrete.res); index.original.res = 1:nrow(a.res)
b.res = a.res[order(a.res[,1]),]; index.original.res = index.original.res[order(a.res[,1])]
# b = tapply(a[,2], a[,1], function(x) max(x))
# b = cbind(as.numeric(names(b)), b)

# b[b[,2] < 0,2] = 0
# index.admiss = diff(c(-1, b[,2], 1)) > 0; index.admiss = c(index.admiss[-length(index.admiss)])

index.admiss.res = 1
s.max.res = b.res[1,2]
for(i in 2:nrow(b.res)){
  s.next.res = b.res[i, 2]
  if(s.next.res > s.max.res){
    s.max.res = s.next.res
    index.admiss.res = c(index.admiss.res, i)
  }
}
d.res = b.res[index.admiss.res,]; index.original.res = index.original.res[index.admiss.res]

#Restricted analysis N = 130------
N.discrete.res.130 = seq(50, 130, 10)
R.discrete.res.130 = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid.res.130 = rep(N.discrete.res.130, each = length(R.discrete.res.130))
R.grid.res.130 = rep(R.discrete.res.130, times = length(N.discrete.res.130))

# par(mfrow = c(1, 2))
# plot(N.grid, R.grid, pch = 3, col = "gray")
Cost.discrete.res.130 = apply(cbind(N.grid.res.130, R.grid.res.130), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.130 = apply(cbind(N.grid.res.130, R.grid.res.130), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.130[EDR.discrete.res.130 < 0] = 0

#find admissible points
a.res.130 = cbind(Cost.discrete.res.130, EDR.discrete.res.130); index.original.res.130 = 1:nrow(a.res.130)
b.res.130 = a.res.130[order(a.res.130[,1]),]; index.original.res.130 = index.original.res.130[order(a.res.130[,1])]
# b = tapply(a[,2], a[,1], function(x) max(x))
# b = cbind(as.numeric(names(b)), b)

# b[b[,2] < 0,2] = 0
# index.admiss = diff(c(-1, b[,2], 1)) > 0; index.admiss = c(index.admiss[-length(index.admiss)])

index.admiss.res.130 = 1
s.max.res.130 = b.res.130[1,2]
for(i in 2:nrow(b.res.130)){
  s.next.res.130 = b.res.130[i, 2]
  if(s.next.res.130 > s.max.res.130){
    s.max.res.130 = s.next.res.130
    index.admiss.res.130 = c(index.admiss.res.130, i)
  }
}
d.res.130 = b.res.130[index.admiss.res.130,]; index.original.res.130 = index.original.res.130[index.admiss.res.130]

#Restricted analysis N=100, R = 1/4------
N.discrete.res.R1 = seq(50, 100, 10)
R.discrete.res.R1 = 6*10^7/c(4)
N.grid.res.R1 = rep(N.discrete.res.R1, each = length(R.discrete.res.R1))
R.grid.res.R1 = rep(R.discrete.res.R1, times = length(N.discrete.res.R1))

Cost.discrete.res.R1 = apply(cbind(N.grid.res.R1, R.grid.res.R1), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.R1 = apply(cbind(N.grid.res.R1, R.grid.res.R1), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.R1[EDR.discrete.res.R1 < 0] = 0

#find admissible points
a.res.R1 = cbind(Cost.discrete.res.R1, EDR.discrete.res.R1); index.original.res.R1 = 1:nrow(a.res.R1)
b.res.R1 = a.res.R1[order(a.res.R1[,1]),]; index.original.res.R1 = index.original.res.R1[order(a.res.R1[,1])]

index.admiss.res.R1 = 1
s.max.res.R1 = b.res.R1[1,2]
for(i in 2:nrow(b.res.R1)){
  s.next.res.R1 = b.res.R1[i, 2]
  if(s.next.res.R1 > s.max.res.R1){
    s.max.res.R1 = s.next.res.R1
    index.admiss.res.R1 = c(index.admiss.res.R1, i)
  }
}
d.res.R1 = b.res.R1[index.admiss.res.R1,]; index.original.res.R1 = index.original.res.R1[index.admiss.res.R1]

#Restricted analysis N=100, R = 1/2------
N.discrete.res.R2 = seq(50, 100, 10)
R.discrete.res.R2 = 6*10^7/c(4, 2)
N.grid.res.R2 = rep(N.discrete.res.R2, each = length(R.discrete.res.R2))
R.grid.res.R2 = rep(R.discrete.res.R2, times = length(N.discrete.res.R2))

Cost.discrete.res.R2 = apply(cbind(N.grid.res.R2, R.grid.res.R2), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.R2 = apply(cbind(N.grid.res.R2, R.grid.res.R2), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.R2[EDR.discrete.res.R2 < 0] = 0

#find admissible points
a.res.R2 = cbind(Cost.discrete.res.R2, EDR.discrete.res.R2); index.original.res.R2 = 1:nrow(a.res.R2)
b.res.R2 = a.res.R2[order(a.res.R2[,1]),]; index.original.res.R2 = index.original.res.R2[order(a.res.R2[,1])]

index.admiss.res.R2 = 1
s.max.res.R2 = b.res.R2[1,2]
for(i in 2:nrow(b.res.R2)){
  s.next.res.R2 = b.res.R2[i, 2]
  if(s.next.res.R2 > s.max.res.R2){
    s.max.res.R2 = s.next.res.R2
    index.admiss.res.R2 = c(index.admiss.res.R2, i)
  }
}
d.res.R2 = b.res.R2[index.admiss.res.R2,]; index.original.res.R2 = index.original.res.R2[index.admiss.res.R2]



c(4, 2, 1, 2/3, 1/2)
#Restricted analysis N=100, R = 1------
N.discrete.res.R3 = seq(50, 100, 10)
R.discrete.res.R3 = 6*10^7/c(4, 2, 1)
N.grid.res.R3 = rep(N.discrete.res.R3, each = length(R.discrete.res.R3))
R.grid.res.R3 = rep(R.discrete.res.R3, times = length(N.discrete.res.R3))

Cost.discrete.res.R3 = apply(cbind(N.grid.res.R3, R.grid.res.R3), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.R3 = apply(cbind(N.grid.res.R3, R.grid.res.R3), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.R3[EDR.discrete.res.R3 < 0] = 0

#find admissible points
a.res.R3 = cbind(Cost.discrete.res.R3, EDR.discrete.res.R3); index.original.res.R3 = 1:nrow(a.res.R3)
b.res.R3 = a.res.R3[order(a.res.R3[,1]),]; index.original.res.R3 = index.original.res.R3[order(a.res.R3[,1])]

index.admiss.res.R3 = 1
s.max.res.R3 = b.res.R3[1,2]
for(i in 2:nrow(b.res.R3)){
  s.next.res.R3 = b.res.R3[i, 2]
  if(s.next.res.R3 > s.max.res.R3){
    s.max.res.R3 = s.next.res.R3
    index.admiss.res.R3 = c(index.admiss.res.R3, i)
  }
}
d.res.R3 = b.res.R3[index.admiss.res.R3,]; index.original.res.R3 = index.original.res.R3[index.admiss.res.R3]


c(4, 2, 1, 2/3, 1/2)
#Restricted analysis N=100, R = 2/3------
N.discrete.res.R4 = seq(50, 100, 10)
R.discrete.res.R4 = 6*10^7/c(4, 2, 1, 2/3)
N.grid.res.R4 = rep(N.discrete.res.R4, each = length(R.discrete.res.R4))
R.grid.res.R4 = rep(R.discrete.res.R4, times = length(N.discrete.res.R4))

Cost.discrete.res.R4 = apply(cbind(N.grid.res.R4, R.grid.res.R4), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.R4 = apply(cbind(N.grid.res.R4, R.grid.res.R4), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.R4[EDR.discrete.res.R4 < 0] = 0

#find admissible points
a.res.R4 = cbind(Cost.discrete.res.R4, EDR.discrete.res.R4); index.original.res.R4 = 1:nrow(a.res.R4)
b.res.R4 = a.res.R4[order(a.res.R4[,1]),]; index.original.res.R4 = index.original.res.R4[order(a.res.R4[,1])]

index.admiss.res.R4 = 1
s.max.res.R4 = b.res.R4[1,2]
for(i in 2:nrow(b.res.R4)){
  s.next.res.R4 = b.res.R4[i, 2]
  if(s.next.res.R4 > s.max.res.R4){
    s.max.res.R4 = s.next.res.R4
    index.admiss.res.R4 = c(index.admiss.res.R4, i)
  }
}
d.res.R4 = b.res.R4[index.admiss.res.R4,]; index.original.res.R4 = index.original.res.R4[index.admiss.res.R4]

c(4, 2, 1, 2/3, 1/2)
#Restricted analysis N=100, R = 1/2------
N.discrete.res.R5 = seq(50, 100, 10)
R.discrete.res.R5 = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid.res.R5 = rep(N.discrete.res.R5, each = length(R.discrete.res.R5))
R.grid.res.R5 = rep(R.discrete.res.R5, times = length(N.discrete.res.R5))

Cost.discrete.res.R5 = apply(cbind(N.grid.res.R5, R.grid.res.R5), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.R5 = apply(cbind(N.grid.res.R5, R.grid.res.R5), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.R5[EDR.discrete.res.R5 < 0] = 0

#find admissible points
a.res.R5 = cbind(Cost.discrete.res.R5, EDR.discrete.res.R5); index.original.res.R5 = 1:nrow(a.res.R5)
b.res.R5 = a.res.R5[order(a.res.R5[,1]),]; index.original.res.R5 = index.original.res.R5[order(a.res.R5[,1])]

index.admiss.res.R5 = 1
s.max.res.R5 = b.res.R5[1,2]
for(i in 2:nrow(b.res.R5)){
  s.next.res.R5 = b.res.R5[i, 2]
  if(s.next.res.R5 > s.max.res.R5){
    s.max.res.R5 = s.next.res.R5
    index.admiss.res.R5 = c(index.admiss.res.R5, i)
  }
}
d.res.R5 = b.res.R5[index.admiss.res.R5,]; index.original.res.R5 = index.original.res.R5[index.admiss.res.R5]










#Admissible points----
pdf("Admissible.pdf", width = 6, height = 8)
par(mgp = c(2, .5, 0))

index.target = as.integer(1:length(index.original) %in% which(EDR.discrete[index.original] > EDR)[1]) + 1
# ii = as.integer(1:length(N.grid) %in% index.original[which(index.target == 2)]) + 1
ii = as.integer(1:length(N.grid) %in% index.original) + 1
plot(Cost.discrete, EDR.discrete, pch = 4,
     col = c("darkgray", "white")[ii], lwd = 1,
     xlab = "Cost", ylab = "EDR", main = "Illustration for admissible design", ylim = c(0.8, 0.95), xlim = c(200000, 400000)
)
# text(Cost.discrete[index.original], EDR.discrete[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
points(Cost.discrete[index.original], EDR.discrete[index.original], pch = 4, col = 2, lwd = 3)
# text(Cost.discrete, EDR.discrete, labels = 1:length(Cost.discrete))
#37, 66, 42
text(Cost.discrete[c(37, 66, 42)], EDR.discrete[c(37, 66, 42)], label = c("A", "B", "C"), pos = 4, offset = 0.3)

# abline(v = Cost.discrete[index.original[which(index.target == 2)]], h = EDR.discrete[index.original[which(index.target == 2)]], lty = 2)
abline(v = Cost.discrete[37], h = EDR.discrete[37], lty = 2)
# v = Cost.discrete[index.original[which(index.target == 2)]]
# h = EDR.discrete[index.original[which(index.target == 2)]]
v = Cost.discrete[37]
h = EDR.discrete[37]
polygon(x = c(0, 0, v, v), y = c(2, h, h, 2), density = 5, border = NA, col = 2)
legend("topright", c("Feasible design", "Admissible design"), pch = 4, col = c("gray", "red"), pt.lwd = c(1,3))
dev.off()


pdf("(Q1-Q4) Cost benefit plot discrete version.pdf", width = 15, height = 8)
# pdf("(Q1) Cost benefit plot discrete version.pdf", width = 12, height = 8)
par(mfcol = c(2, 5), mgp = c(2.8, 1.2, 0), mar = c(4, 5, 3, 1) + .1, cex.axis = 1.8, cex.lab = 2, cex.main = 3)
# par(mgp = c(2, .5, 0), mar = c(3, 4, 2, 1) + .1)

# cex.a = .8
cex.a = 2
#T1-----
# pdf("(Q2) Cost benefit plot discrete version.pdf", width = 12, height = 8)
# par(mfrow = c(1, 2), mgp = c(2, .5, 0))
Cost = 200000
plot(Cost.discrete, EDR.discrete, pch = 4, ylim = c(0.75, 0.95), cex = 2,
     #      xlim = c(min(Cost.discrete), 600000),
     xlim = c(min(Cost.discrete), 300000),
     col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
     #      xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T1")
     xlab = "C", ylab = "EDR", main = "T1", axes = F); box()
axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)
index.target = as.integer(1:length(index.original) %in% tail(which(Cost.discrete[index.original] <= Cost), 1)) + 1
# text(Cost.discrete[index.original], EDR.discrete[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
points(Cost.discrete[index.original], EDR.discrete[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])
abline(v = Cost, col = 4); #lines(d[,1], d[,2], col = 2)
# legend("bottomright", c("Feasible design", "Admissible design", "Optimal design"), pch = 4, col = c("gray55", "black", "red"))
plot(N.grid, R.grid, pch = 4, col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
     #      xlab = "Sample size (N)", ylab = "Sequencing depth (R)", axes = F); box()
     xlab = "N", ylab = "R", axes = F, cex = 2,
); box()
axis(1); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
# text(N.grid[index.original], R.grid[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
points(N.grid[index.original], R.grid[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])
# curve(cost.R.function, col = 4, add = T)
# Cost = Cost.discrete[index.original][which(EDR.discrete[index.original] > EDR)[1]]; curve(cost.R.function, col = "green", add = T)
# dev.off()


#T2-----
# h = layout(rbind(matrix(1:8, 2), c(9, 9, 9, 9)), heights = c(1, 1, .3)); #layout.show(h)
EDR = 0.85
plot(Cost.discrete, EDR.discrete, pch = 4, ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 300000),
     #      col = c("gray55", "red")[as.integer(1:length(N.grid) %in% index.original) + 1],
     col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
     #      xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T2")
     xlab = "C", ylab = "EDR", main = "T2", cex = 2, axes = F); box()
axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)
index.target = as.integer(1:length(index.original) %in% which(EDR.discrete[index.original] > EDR)[1]) + 1
# text(Cost.discrete[index.original], EDR.discrete[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
points(Cost.discrete[index.original], EDR.discrete[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])
abline(h = EDR, col = 4); #lines(d[,1], d[,2], col = 2)
# legend("bottomright", c("Feasible design", "Admissible design", "Optimal design"), pch = 4, col = c("gray55", "black", "red"))
plot(N.grid, R.grid, pch = 4, col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
     #      xlab = "Sample size (N)", ylab = "Sequencing depth (R)", axes = F); box()
     xlab = "N", ylab = "R", axes = F, cex = 2); box()
axis(1); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
# text(N.grid[index.original], R.grid[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
points(N.grid[index.original], R.grid[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(2, 2)[index.target], lwd = c(1, 2)[index.target])
# Power = EDR; curve(power.R.function, col = 4, add = T)
# Cost = Cost.discrete[index.original][which(EDR.discrete[index.original] > EDR)[1]]; curve(cost.R.function, col = "green", add = T)
# dev.off()


#T3-----
# pdf("(Q3) Cost benefit plot discrete version.pdf", width = 12, height = 8)
# par(mfrow = c(1, 2), mgp = c(2, .5, 0))
EDR = 0.85
plot(Cost.discrete.res, EDR.discrete.res, pch = 1, cex = 2, ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000),
     #      col = c("gray55", "red")[as.integer(1:length(N.grid) %in% index.original) + 1],
     col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
     #      xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T3")
     xlab = "C", ylab = "EDR", main = "T3", axes = F); box()
axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)
index.target.res = as.integer(1:length(index.original.res) %in% which(EDR.discrete.res[index.original.res] > EDR)[1]) + 1
index.target.res[length(index.target.res)] = 2
# text(Cost.discrete.res[index.original.res], EDR.discrete.res[index.original.res], labels = 1:length(index.original.res), col = c(1, 2)[index.target.res], cex = c(1, 1)[index.target.res])
points(Cost.discrete.res[index.original.res], EDR.discrete.res[index.original.res], pch = 1, col = c(1, 2)[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
text(Cost.discrete.res[index.original.res][index.target.res == 2],
     EDR.discrete.res[index.original.res][index.target.res == 2], labels = c("B", "A"),
     col = 2, cex = 2, lwd = c(1, 2)[index.target.res], pos = 3)
# abline(h = EDR, col = 4); #lines(d.res[,1], d.res[,2], col = 2)
# legend("bottomright", c("Feasible design", "Admissible design", "Optimal design"), pch = 4, col = c("gray55", "black", "red"))
plot(N.grid.res, R.grid.res, pch = 1, cex = 2, col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
     #      xlab = "Sample size (N)", ylab = "Sequencing depth (R)", axes = F, xlim = c(50, 200)); box()
     xlab = "N", ylab = "R", axes = F, xlim = c(50, 200)); box()
axis(1); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
# text(N.grid.res[index.original.res], R.grid.res[index.original.res], labels = 1:length(index.original.res), col = c(1, 2)[index.target.res], cex = c(1, 1)[index.target.res])
points(N.grid.res[index.original.res], R.grid.res[index.original.res], pch = 1, col = c(1, 2)[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
text(N.grid.res[index.original.res][index.target.res == 2],
     R.grid.res[index.original.res][index.target.res == 2], labels = c("B", "A"),
     col = 2, cex = 2, lwd = c(1, 2)[index.target.res], pos = 4)
abline(v = 80, col = 4)
# Power = EDR; curve(power.R.function, col = 4, add = T)
# Cost = Cost.discrete[index.original][which(EDR.discrete[index.original] > EDR)[1]]; curve(cost.R.function, col = "green", add = T)

# dev.off()



#T4-----
# pdf("(Q4) Cost benefit plot discrete version.pdf", width = 12, height = 8)
# par(mfrow = c(1, 2), mgp = c(2, .5, 0))
EDR = 0.85
if(0){
  plot(Cost.discrete, EDR.discrete, pch = 4, ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000),
       #      col = c("gray55", "red")[as.integer(1:length(N.grid) %in% index.original) + 1],
       col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       #      xlab = "Cost", ylab = "EDR", main = "Q4")
       xlab = "C", ylab = "EDR", main = "T4", axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)
  index.target = as.integer(1:length(index.original) %in% which(EDR.discrete[index.original] > EDR)[1]) + 1
  # text(Cost.discrete[index.original], EDR.discrete[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
  points(Cost.discrete[index.original], EDR.discrete[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(1, 1)[index.target], lwd = c(1, 2)[index.target])
  abline(h = EDR, col = 4); #lines(d[,1], d[,2], col = 2)
  # legend("bottomright", c("Feasible design", "Admissible design", "Admissible design (restricted)", "Optimal design", "Optimal design (restricted)"), pch = 4, col = c("gray55", "black", "purple", "red", "orange"))
  legend("bottomright", c("Feasible design", "Feasible design (restricted)", "Admissible design", "Admissible design (restricted)", "Optimal design", "Optimal design (restricted)"),
         pch = c(4, 1, 4, 1, 4, 1), col = c("gray55", "gray55", "black", "black", "red", "red"), cex = 1)
}

plot(Cost.discrete.res[-index.original.res], EDR.discrete.res[-index.original.res], pch = 1,
     ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 1000000),
     #      col = c("gray55", "red")[as.integer(1:length(N.grid) %in% index.original) + 1],
     #        col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
     col = c("gray55"), cex = 2,
     #        xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T4")
     xlab = "C", ylab = "EDR", main = "T4", axes = F); box()
axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)
# index.target.res = as.integer(1:length(index.original.res) %in% which(EDR.discrete.res[index.original.res] > EDR)[1]) + 1
index.target.res = c(rep(1, length(index.original.res)-1), 2)
# text(Cost.discrete.res[index.original.res], EDR.discrete.res[index.original.res], labels = 1:length(index.original.res), col = c("purple", "orange")[index.target.res], cex = c(1, 1)[index.target.res])
points(Cost.discrete.res[index.original.res], EDR.discrete.res[index.original.res], pch = 1, col = c("black", "red")[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
# lines(d.res[,1], d.res[,2], col = 2)
# legend("bottomright", c("Feasible design", "Feasible design (restricted)", "Admissible design", "Admissible design (restricted)", "Optimal design", "Optimal design (restricted)"),
#        pch = c(4, 1, 4, 1, 4, 1), col = c("gray55", "gray55", "black", "black", "red", "red"), cex = 1)


# arrows(x0 = Cost.discrete.res[index.original.res][12], y0 = EDR.discrete.res[index.original.res][12],
#        x1 = Cost.discrete[index.original][16]*1.1, y1 = EDR.discrete[index.original][16]*0.999, col = 1,
#        length = 0.1, lwd = 3)
#restricted N130
points(Cost.discrete.res.130[-index.original.res.130], EDR.discrete.res.130[-index.original.res.130], pch = 4, ylim = c(0.75, 0.95),
       #      col = c("gray55", "red")[as.integer(1:length(N.grid) %in% index.original) + 1],
       #        col = c("gray55", "white")[as.integer(1:length(N.grid.res.130) %in% index.original.res.130) + 1],
       col = c("gray55"), cex = 2)
# index.target.res.130 = as.integer(1:length(index.original.res.130) %in% which(EDR.discrete.res.130[index.original.res.130] > EDR)[1]) + 1
# index.target.res.130 = c(rep(1, length(index.original.res.130)-1), 2)
index.target.res.130 = rep(1, length(index.original.res.130)); index.target.res.130[16] = 2
# text(Cost.discrete.res.130[index.original.res.130], EDR.discrete.res.130[index.original.res.130], labels = 1:length(index.original.res.130), col = c("purple", "orange")[index.target.res.130], cex = c(1, 1)[index.target.res.130])
points(Cost.discrete.res.130[index.original.res.130], EDR.discrete.res.130[index.original.res.130], pch = 4, col = c("black", "red")[index.target.res.130],
       cex = c(2, 2)[index.target.res.130], lwd = c(1, 2)[index.target.res.130])

112000; 0.88
104000; 0.92

#second figure
if(0){
  plot(N.grid, R.grid, pch = 4, col = c("gray55", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "Sample size (N)", ylab = "Sequencing depth (R)", axes = F, cex = 1); box()
  axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
  # text(N.grid[index.original], R.grid[index.original], labels = 1:length(index.original), col = c(1, 2)[index.target], cex = c(1, 1)[index.target])
  points(N.grid[index.original], R.grid[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(1, 1)[index.target], lwd = c(1, 2)[index.target])
  Power = EDR; curve(power.R.function, col = 4, add = T)
}
#restricted N100
# points(N.grid.res, R.grid.res, pch = 1, col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
plot(N.grid.res[-index.original.res], R.grid.res[-index.original.res], pch = 1, col = "gray55", cex = 2,
     #        xlab = "Sample size (N)", ylab = "Sequencing depth (R)", xlim = c(50, 200), ylim = range(R.grid), axes = F); box()
     xlab = "N", ylab = "R", xlim = c(50, 200), ylim = range(R.grid), axes = F); box()
axis(1); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
# text(N.grid.res[index.original.res], R.grid.res[index.original.res], labels = 1:length(index.original.res), col = c("purple", "orange")[index.target.res], cex = c(1, 1)[index.target.res])
points(N.grid.res[index.original.res], R.grid.res[index.original.res], pch = 1, col = c("black", "red")[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
# text(N.grid.res[index.original.res], R.grid.res[index.original.res], label = 1:length(index.original.res), col = c("black", "red")[index.target.res], cex = c(2, 2)[index.target.res], lwd = c(1, 2)[index.target.res])
# Cost = Cost.discrete[index.original][which(EDR.discrete[index.original] > EDR)[1]]; curve(cost.R.function, col = "green", add = T)

#restricted N130
# points(N.grid.res, R.grid.res, pch = 1, col = c("gray55", "white")[as.integer(1:length(N.grid.res) %in% index.original.res) + 1],
points(N.grid.res.130[-index.original.res.130], R.grid.res.130[-index.original.res.130], pch = 4, col = "gray55", cex = 2)
# text(N.grid.res.130[index.original.res.130], R.grid.res.130[index.original.res.130], labels = 1:length(index.original.res.130), col = c("purple", "orange")[index.target.res.130], cex = c(1, 1)[index.target.res.130])
points(N.grid.res.130[index.original.res.130], R.grid.res.130[index.original.res.130], pch = 4, col = c("black", "red")[index.target.res.130],
       cex = c(2, 2)[index.target.res.130], lwd = c(1, 2)[index.target.res.130])
# Cost = Cost.discrete[index.original][which(EDR.discrete[index.original] > EDR)[1]]; curve(cost.R.function, col = "green", add = T)
abline(v = c(80, 130), col = 4)

#legend
if(0){
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes = F, xlab = "", ylab = ""); box()
  legend("center", c("Feasible design", "Feasible design (restricted)", "Admissible design", "Admissible design (restricted)", "Optimal design", "Optimal design (restricted)"),
         pch = c(4, 1, 4, 1, 4, 1), col = c("gray55", "gray55", "black", "black", "red", "red"), horiz = T, bty = "n", cex = 2)
}


#T5-----
#R = 1/4
if(0){
  index.target.res.R1 = as.integer(1:length(index.original.res.R1) %in% which(EDR.discrete.res.R1[index.original.res.R1] > EDR)[1]) + 1
  index.target.res.R2 = as.integer(1:length(index.original.res.R2) %in% which(EDR.discrete.res.R2[index.original.res.R2] > EDR)[1]) + 1
  index.target.res.R3 = as.integer(1:length(index.original.res.R3) %in% which(EDR.discrete.res.R3[index.original.res.R3] > EDR)[1]) + 1
  index.target.res.R4 = as.integer(1:length(index.original.res.R4) %in% which(EDR.discrete.res.R4[index.original.res.R4] > EDR)[1]) + 1
  index.target.res.R5 = as.integer(1:length(index.original.res.R5) %in% which(EDR.discrete.res.R5[index.original.res.R5] > EDR)[1]) + 1

  plot(Cost.discrete, EDR.discrete, pch = 4, ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000),
       col = c("white", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T1")
  text(Cost.discrete, EDR.discrete, labels = 1:length(Cost.discrete))

  plot(N.grid, R.grid, pch = 4, col = c("white", "white")[as.integer(1:length(N.grid) %in% index.original) + 1],
       xlab = "Sample size (N)", ylab = "Sequencing depth (R)", axes = F); box()
  axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))
  text(N.grid, R.grid, labels = 1:length(N.grid))
  # points(N.grid[index.original], R.grid[index.original], pch = 4, col = c(1, 2)[index.target], cex = c(1, 1)[index.target], lwd = c(1, 2)[index.target])
}

plot(Cost.discrete[6:10], EDR.discrete[6:10], pch = 1:5, col = 2, cex = 2, lwd = 2,
     main = "T5", ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000),
     #      xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)"
     xlab = "C", ylab = "EDR", axes = F); box()
axis(1, at = axTicks(1), labels = axTicks(1)/5); axis(2)
if(0){
  points(Cost.discrete.res.R2[index.original.res.R2][which(index.target.res.R2 == 2)],
         EDR.discrete.res.R2[index.original.res.R2][which(index.target.res.R2 == 2)], pch = 2, col = 2, cex = 2, lwd = 2,
         xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T5", ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000)
  )
  points(Cost.discrete.res.R3[index.original.res.R3][which(index.target.res.R3 == 2)],
         EDR.discrete.res.R3[index.original.res.R3][which(index.target.res.R3 == 2)], pch = 3, col = 2, cex = 2, lwd = 2,
         xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T5", ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000)
  )
  points(Cost.discrete.res.R4[index.original.res.R4][which(index.target.res.R4 == 2)],
         EDR.discrete.res.R4[index.original.res.R4][which(index.target.res.R4 == 2)], pch = 4, col = 2, cex = 2, lwd = 2,
         xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T5", ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000)
  )
  points(Cost.discrete.res.R5[index.original.res.R5][which(index.target.res.R5 == 2)],
         EDR.discrete.res.R5[index.original.res.R5][which(index.target.res.R5 == 2)], pch = 5, col = 2, cex = 2, lwd = 2,
         xlab = "Cost (C)", ylab = "Expected discovery rate (EDR)", main = "T5", ylim = c(0.75, 0.95), xlim = c(min(Cost.discrete), 600000)
  )
}

plot(rep(60, 5), unique(R.grid), pch = 1:5, col = "red", cex = 2, lwd = 2,
     #      xlab = "Sample size (N)", ylab = "Sequencing depth (R)", xlim = c(50, 200), ylim = range(R.grid), axes = F); box()
     xlab = "N", ylab = "R", xlim = c(50, 200), ylim = range(R.grid), axes = F); box()
axis(1); axis(2, at = R.discrete, cex.axis = cex.a, labels = c("1/4 Lane", "1/2 Lane", "1 Lane", "1.5 Lane", "2 Lane"))

if(0){
  points(N.grid.res.R2[index.original.res.R2][which(index.target.res.R2 == 2)], R.grid.res.R2[index.original.res.R2][which(index.target.res.R2 == 2)],
         pch = 2, col = "red", cex = 2, lwd = 2,); box()

  points(N.grid.res.R3[index.original.res.R3][which(index.target.res.R3 == 2)], R.grid.res.R3[index.original.res.R3][which(index.target.res.R3 == 2)],
         pch = 3, col = "red", cex = 2, lwd = 2,
         xlab = "Sample size (N)", ylab = "Sequencing depth (R)", xlim = c(50, 200), ylim = range(R.grid)); box()

  points(N.grid.res.R4[index.original.res.R4][which(index.target.res.R4 == 2)], R.grid.res.R4[index.original.res.R4][which(index.target.res.R4 == 2)],
         pch = 4, col = "red", cex = 2, lwd = 2,
         xlab = "Sample size (N)", ylab = "Sequencing depth (R)", xlim = c(50, 200), ylim = range(R.grid)); box()

  points(N.grid.res.R5[index.original.res.R5][which(index.target.res.R5 == 2)], R.grid.res.R5[index.original.res.R5][which(index.target.res.R5 == 2)], pch = 5, col = "red", cex = 2, lwd = 2,
         xlab = "Sample size (N)", ylab = "Sequencing depth (R)", xlim = c(50, 200), ylim = range(R.grid)); box()
}
dev.off()









# function ----------------------------------------------------------------
curve.fitting = function(n, r, EDR){
  lkh<-function(x){
    return(sum((-x[3]*r^(-x[4])-x[1]*n^(-x[2])+1-EDR)^2))
  }
  #   outOptim<-DEoptim(lkh,lower=c(0, 0, 0, 0),upper=c(100, 100, 10^6, 100), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
  outOptim<-DEoptim(lkh,lower=c(0, 0, 0, 0),upper=c(100, 100, 10^8, 100), DEoptim.control(trace = F, itermax = 500))###stochastic fitting
  return(outOptim$optim$bestmem)
}

Cost.N.R = function(N, R) N*(Cost.N + Cost.R*R)

EDR.curve = function(n, r, x) 1-x[1]*n^(-x[2])-x[3]*r^(-x[4])

Optim.Function <- function(a,b,Cost, EDR.curve.para){
  fBUM <- function(z) {
    N=z
    Target.Func=-(1-EDR.curve.para[1]*N^(-EDR.curve.para[2])-EDR.curve.para[3]*((Cost/N-a)/b)^(-EDR.curve.para[4]))
    return(Target.Func)
  }
  pars = optim(6, fBUM,method= "L-BFGS-B",upper=Cost/(a+b*10^6),lower=5, control=list(trace = F))$par
  N.opt=round(pars)
  R=(Cost/N.opt-a)/b
  Target.Func=(1-EDR.curve.para[1]*N.opt^(-EDR.curve.para[2])-EDR.curve.para[3]*((Cost/N.opt-a)/b)^(-EDR.curve.para[4]))
  Cost.opt=N.opt*(a+b*R)

  return(c(N.star=N.opt,R.star=R,EDR.star=Target.Func,Cost.star=Cost.opt))
}

Optim.Function.restrict.N <- function(a,b,Cost, EDR.curve.para, N.restrict){
  fBUM <- function(z) {
    N=z
    Target.Func=-(1-EDR.curve.para[1]*N^(-EDR.curve.para[2])-EDR.curve.para[3]*((Cost/N-a)/b)^(-EDR.curve.para[4]))
    return(Target.Func)
  }
  pars = optim(6, fBUM,method= "L-BFGS-B",upper=min(N.restrict, Cost/(a+b*10^6)),lower=5, control=list(trace = F))$par
  N.opt=round(pars)
  R=(Cost/N.opt-a)/b
  Target.Func=(1-EDR.curve.para[1]*N.opt^(-EDR.curve.para[2])-EDR.curve.para[3]*((Cost/N.opt-a)/b)^(-EDR.curve.para[4]))
  Cost.opt=N.opt*(a+b*R)

  return(c(N.star=N.opt,R.star=R,EDR.star=Target.Func,Cost.star=Cost.opt))
}

Optim.Function.2 <- function(a,b,Power, EDR.curve.para){
  fBUM <- function(z) {
    N=z
    R = ((Power-1+EDR.curve.para[1]*N^(-EDR.curve.para[2]))/(-EDR.curve.para[3]))^(-1/EDR.curve.para[4])
    Target.Func=N*(a+b*R)
    return(Target.Func)
  }
  pars = optim(6, fBUM,method= "L-BFGS-B",upper=100,lower=((1-Power)/EDR.curve.para[1])^(-1/EDR.curve.para[2])+1, control=list(trace = F))$par #minimize cost
  N.opt=round(pars)
  R=((Power-1+EDR.curve.para[1]*N.opt^(-EDR.curve.para[2]))/(-EDR.curve.para[3]))^(-1/EDR.curve.para[4])
  Target.Func=N.opt*(a+b*R) #cost
  Power.opt=(1-EDR.curve.para[1]*N.opt^(-EDR.curve.para[2])-EDR.curve.para[3]*R^(-EDR.curve.para[4])) #power

  return(c(N.star=N.opt,R.star=R,EDR.star=Power.opt,Cost.star=Target.Func))
}

# setwd(paste("/Users/Masaki/Dropbox/Research/Power calculation/New4/Dispersion", dispersion, "/", sep = ""))
#Parametric bootstrap
library(snowfall)
sfLibrary(DEoptim)
sfInit(parallel=T, cpus=3)
# for(fold.change in 1:4){
#=========================================================#
surface.para = lapply(1:4, function(fold.change){
  EDR.PS.2 = lapply(samplesize.pilot, function(j){
    EDR = lapply(1:n.repeat, function(i){ #for each repeat
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".Storey.RData",sep=""))
      load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.RData",sep=""))

      EDR = t(matrix(apply(sapply(Pilot.power$Result, function(x) sapply(x, function(y) y[3,])), 1, mean), nrow = length(samplesize.predict)))
      rownames(EDR) = read; colnames(EDR) = paste("Predict", samplesize.predict, sep = "")
      return(EDR)
    })
    names(EDR) = paste("Repeat", 1:n.repeat, sep = "")
    return(EDR)
  })
  names(EDR.PS.2) = paste("Pilot", samplesize.pilot, sep = "")

  EDR.PS.2 = lapply(EDR.PS.2, function(x){
    y = matrix(rowMeans(sapply(x, function(y) y)), nrow = length(read))
    rownames(y) = read
    colnames(y) = paste("Predict", samplesize.predict, sep = "")
    return(y)
  })
  #==================================================#
  #True
  #==================================================#
  EDR.T = t(sapply(1:length(read.predict), function(r){
    load(paste("True.Dispersion.50.Depth.", r, "M.rdata", sep = ""))
    EDR.T = sapply(True.Dispersion.50[[fold.change]], function(x) mean(x[,1]))
    EDR.T = EDR.T[-1]
    names(EDR.T) = samplesize.predict
    return(EDR.T)
  }))

  rownames(EDR.T) = read.predict
  #==================================================#
  #Surface fitting
  #==================================================#

  #==================================================#
  #True
  #==================================================#
  #   NR.grid = expand.grid(seq(10, 100, length = 500), seq(10^5, 2*10^6, length = 500))
  N = rep(samplesize.predict, each = length(read))
  R = rep(read, length(samplesize.predict))

  predict.N = rep(samplesize.predict, each = length(read.predict))
  predict.R = rep(read.predict, length(samplesize.predict))

  pilot.surface.para = sapply(EDR.PS.2, function(y) EDR.curve.para = curve.fitting(n = N, r = R, EDR = y))

  true.para = curve.fitting(n = predict.N, r = predict.R, EDR = EDR.T)

  return(list(pilot = pilot.surface.para, true = true.para))
})
names(surface.para) = paste("Fold.change.", 1:4, sep = "")
# save(surface.para, file = paste("(Dispersion", dispersion, ")Surface.parameter.constraint.RData", sep = ""))
# save(surface.para, file = paste("(Dispersion", dispersion, ")Surface.parameter.constraint.Storey.RData", sep = ""))
save(surface.para, file = paste("(Dispersion", dispersion, ")Surface.parameter.constraint.CDD.RData", sep = ""))

#=========================================================#
sfExportAll()
# for(fold.change in 1:4){
# EDR.result = sfLapply(1:4, function(fold.change){
EDR.result.CDD.cost = lapply(1:4, function(fold.change){
  print(fold.change)
  # EDR.result.Storey = lapply(1:4, function(fold.change){
  # EDR.result.down = lapply(1:4, function(fold.change){
  EDR.PS.2 = lapply(samplesize.pilot, function(j){
    EDR = lapply(1:n.repeat, function(i){ #for each repeat
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".Storey.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".Storey.new.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".Storey.3M.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".Storey.1M.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.3M.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.1M.RData",sep=""))
      load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.new.RData",sep=""))
      EDR = t(matrix(apply(sapply(Pilot.power$Result, function(x) sapply(x, function(y) y[3,])), 1, mean), nrow = length(samplesize.predict)))
      rownames(EDR) = read; colnames(EDR) = paste("Predict", samplesize.predict, sep = "")
      return(EDR)
    })
    names(EDR) = paste("Repeat", 1:n.repeat, sep = "")
    return(EDR)
  })
  names(EDR.PS.2) = paste("Pilot", samplesize.pilot, sep = "")

  #==================================================#
  #True
  #==================================================#
  EDR.T = t(sapply(1:length(read.predict), function(r){
    load(paste("True.Dispersion.50.Depth.", r, "M.rdata", sep = ""))
    EDR.T = sapply(True.Dispersion.50[[fold.change]], function(x) mean(x[,1]))
    EDR.T = EDR.T[-1]
    names(EDR.T) = samplesize.predict
    return(EDR.T)
  }))

  rownames(EDR.T) = read.predict
  #==================================================#
  #Surface fitting
  #==================================================#
  #pilot
  N = rep(samplesize.predict, each = length(read))
  R = rep(read, length(samplesize.predict))

  #true
  predict.N = rep(samplesize.predict, each = length(read.predict))
  predict.R = rep(read.predict, length(samplesize.predict))

  true.para = curve.fitting(n = predict.N, r = predict.R, EDR = EDR.T)

  #pilot optim
  NR.optim.pilot = lapply(EDR.PS.2, function(x){ #for each pilot data sample size
    #   NR.optim.pilot = lapply(1:length(EDR.PS.2), function(xx){ #for each pilot data sample size
    #     x = EDR.PS.2[[xx]]; print(xx)
    lapply(x, function(y){ #for each repeat
      #     lapply(1:length(x), function(yy){ #for each repeat
      #derive parameter of inverse power law curve
      #       y = x[[yy]]; print(yy)
      EDR.curve.para = curve.fitting(n = N, r = R, EDR = y)

      NR = Optim.Function.2(a = Cost.N, b = Cost.R, Power = Power, EDR.curve.para)
      # N.optim = NR[1]; R.optim = NR[2]; EDR.optim = NR[3]; cost.optim = NR[4]

      EDR.curve.predict = matrix(apply(cbind(N, R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read))
      rownames(EDR.curve.predict) = read
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      mse = sqrt(mean((y-EDR.curve.predict)[, -1]^2))

      EDR.curve.predict = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read.predict))
      rownames(EDR.curve.predict) = read.predict
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      return(list(optim = NR, predict = EDR.curve.predict, mse = mse, EDR.curve.para = EDR.curve.para))
    })
  })

  #true optim
  NR.optim.true = Optim.Function.2(a = Cost.N, b = Cost.R, Power = Power, EDR.curve.para = true.para)

  EDR.curve.true = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], true.para)), nrow = length(read.predict))
  rownames(EDR.curve.true) = read.predict
  colnames(EDR.curve.true) = paste("Predict", samplesize.predict, sep = "")

  #   library(rgl); open3d()
  #   plot3d(x = predict.N, y = predict.R, z = EDR.curve.true, add = F, col = 2, type = "p", xlab = "N", ylab = "R", zlab = "EDR", ylim = c(0, 1))
  #   plot3d(x = predict.N, y = predict.R, z = EDR.T, add = T, col = 1, type = "p", xlab = "N", ylab = "R", zlab = "EDR")
  mse.pilot = sapply(NR.optim.pilot, function(x) sapply(x, function(y) y$mse))
  mse.T = sqrt(mean((EDR.T-EDR.curve.true)[, -1]^2)); mse.T #drop predict 10, the worst case of curve fitting

  mse.pilot.T.fitted = sapply(NR.optim.pilot, function(x) sapply(x, function(y) sqrt(mean((y$predict-EDR.curve.true)[, -1]^2))))


  return(list(mse.pilot = mse.pilot, mse.T = mse.T, true.para = true.para, mse.pilot.T.fitted = mse.pilot.T.fitted, NR.optim.true = NR.optim.true, NR.optim.pilot = NR.optim.pilot))
})

read.list = list(design1 = seq(10^6, 8*10^6, length = 6), design2 = seq(10^6, 4*10^6, length = 6), design3 = seq(10^6, 2*10^6, length = 6), design4 = seq(10^5, 10^6, length = 6))

#=========================================================#
#               different pilot design
#=========================================================#
EDR.result.pilot.design = lapply(1:4, function(design){
  read = read.list[[design]]
  print(design)
  # EDR.result.Storey = lapply(1:4, function(fold.change){
  # EDR.result.down = lapply(1:4, function(fold.change){
  EDR.PS.2 = lapply(1, function(j){ #for different pilot sample size
    EDR = lapply(1:n.repeat, function(i){ #for each repeat
      load(paste("./PilotDesign/(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", c(2, 4, 8, 16)[design], ".predict.5-100.read.1M-6M.M10.2way.Fold.2.dispersion.", dispersion, ".CDD.", c(8, 4, 2, 1)[design], "M.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.new.RData",sep=""))
      EDR = t(matrix(apply(sapply(Pilot.power$Result, function(x) sapply(x, function(y) y[3,])), 1, mean), nrow = length(samplesize.predict)))
      rownames(EDR) = read; colnames(EDR) = paste("Predict", samplesize.predict, sep = "")
      return(EDR)
    })
    names(EDR) = paste("Repeat", 1:n.repeat, sep = "")
    return(EDR)
  })
  names(EDR.PS.2) = paste("Pilot", c(2, 4, 8, 16)[design], sep = "")

  #==================================================#
  #True
  #==================================================#
  EDR.T = t(sapply(1:length(read.predict), function(r){
    load(paste("True.Dispersion.50.Depth.", r, "M.rdata", sep = ""))
    EDR.T = sapply(True.Dispersion.50[[2]], function(x) mean(x[,1]))
    EDR.T = EDR.T[-1]
    names(EDR.T) = samplesize.predict
    return(EDR.T)
  }))

  rownames(EDR.T) = read.predict
  #==================================================#
  #Surface fitting
  #==================================================#
  #pilot
  N = rep(samplesize.predict, each = length(read))
  R = rep(read, length(samplesize.predict))

  #true
  predict.N = rep(samplesize.predict, each = length(read.predict))
  predict.R = rep(read.predict, length(samplesize.predict))

  true.para = curve.fitting(n = predict.N, r = predict.R, EDR = EDR.T)

  #pilot optim
  NR.optim.pilot = lapply(EDR.PS.2, function(x){ #for each pilot data sample size
    lapply(x, function(y){ #for each repeat
      #derive parameter of inverse power law curve
      EDR.curve.para = curve.fitting(n = N, r = R, EDR = y)

      NR = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para)
      # N.optim = NR[1]; R.optim = NR[2]; EDR.optim = NR[3]; cost.optim = NR[4]
      EDR.pilot.NR.on.true = EDR.curve(NR[1], NR[2], true.para)
      EDR.curve.predict = matrix(apply(cbind(N, R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read))
      rownames(EDR.curve.predict) = read
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      mse = sqrt(mean((y-EDR.curve.predict)[, -1]^2))

      EDR.curve.predict = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read.predict))
      rownames(EDR.curve.predict) = read.predict
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      return(list(optim = NR, predict = EDR.curve.predict, EDR.pilot.NR.on.true = EDR.pilot.NR.on.true, mse = mse, EDR.curve.para = EDR.curve.para))
    })
  })

  #true optim
  NR.optim.true = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para = true.para)

  EDR.curve.true = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], true.para)), nrow = length(read.predict))
  rownames(EDR.curve.true) = read.predict
  colnames(EDR.curve.true) = paste("Predict", samplesize.predict, sep = "")

  #   library(rgl); open3d()
  #   plot3d(x = predict.N, y = predict.R, z = EDR.curve.true, add = F, col = 2, type = "p", xlab = "N", ylab = "R", zlab = "EDR", ylim = c(0, 1))
  #   plot3d(x = predict.N, y = predict.R, z = EDR.T, add = T, col = 1, type = "p", xlab = "N", ylab = "R", zlab = "EDR")
  mse.pilot = sapply(NR.optim.pilot, function(x) sapply(x, function(y) y$mse))
  mse.T = sqrt(mean((EDR.T-EDR.curve.true)[, -1]^2)); mse.T #drop predict 10, the worst case of curve fitting

  mse.pilot.T.fitted = sapply(NR.optim.pilot, function(x) sapply(x, function(y) sqrt(mean((y$predict-EDR.curve.true)[, -1]^2))))


  return(list(mse.pilot = mse.pilot, mse.T = mse.T, mse.pilot.T.fitted = mse.pilot.T.fitted, NR.optim.true = NR.optim.true, NR.optim.pilot = NR.optim.pilot))
})

#=========================================================
#compare RNAseqPower, scotty, fix cost, find optimal
#=========================================================
Cost = 80000
EDR.result.RNAseqP = lapply(1, function(fold.change){
  # EDR.result.Storey = lapply(1:4, function(fold.change){
  # EDR.result.down = lapply(1:4, function(fold.change){
  EDR.PS.2 = lapply(1:3, function(j){ #for different pilot sample size
    load(paste("RNASeqPower.N.", c(2, 4, 8)[j], ".rdata", sep = ""))
    EDR = get(paste("RNASeqPower.N.", c(2, 4, 8)[j], sep = ""))
    EDR = lapply(EDR, function(x){
      x = x[1:6, -1] #drop predict 5
      read = as.numeric(rownames(x))*10000
      rownames(x) = read
      colnames(x) = paste("Predict", samplesize.predict, sep = "")
      return(x)
    })
    names(EDR) = paste("Repeat", 1:length(EDR), sep = "")
    return(EDR)
  })
  read = as.numeric(rownames(EDR.PS.2[[1]][[1]]))
  #==================================================#
  #True
  #==================================================#
  EDR.T = t(sapply(1:length(read.predict), function(r){
    load(paste("True.Dispersion.50.Depth.", r, "M.rdata", sep = ""))
    EDR.T = sapply(True.Dispersion.50[[2]], function(x) mean(x[,1]))
    EDR.T = EDR.T[-1]
    names(EDR.T) = samplesize.predict
    return(EDR.T)
  }))

  rownames(EDR.T) = read.predict
  #==================================================#
  #Surface fitting
  #==================================================#
  #pilot
  N = rep(samplesize.predict, each = length(read))
  R = rep(read, length(samplesize.predict))

  #true
  predict.N = rep(samplesize.predict, each = length(read.predict))
  predict.R = rep(read.predict, length(samplesize.predict))

  true.para = curve.fitting(n = predict.N, r = predict.R, EDR = EDR.T)

  #pilot optim
  NR.optim.pilot = lapply(EDR.PS.2, function(x){ #for each pilot data sample size
    lapply(x, function(y){ #for each repeat
      #derive parameter of inverse power law curve
      EDR.curve.para = curve.fitting(n = N, r = R, EDR = y)

      NR = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para)
      # N.optim = NR[1]; R.optim = NR[2]; EDR.optim = NR[3]; cost.optim = NR[4]
      EDR.pilot.NR.on.true = EDR.curve(NR[1], NR[2], true.para)
      EDR.curve.predict = matrix(apply(cbind(N, R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read))
      rownames(EDR.curve.predict) = read
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      mse = sqrt(mean((y-EDR.curve.predict)[, -1]^2))

      EDR.curve.predict = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read.predict))
      rownames(EDR.curve.predict) = read.predict
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      return(list(optim = NR, predict = EDR.curve.predict, EDR.pilot.NR.on.true = EDR.pilot.NR.on.true, mse = mse, EDR.curve.para = EDR.curve.para))
    })
  })

  #true optim
  NR.optim.true = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para = true.para)

  EDR.curve.true = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], true.para)), nrow = length(read.predict))
  rownames(EDR.curve.true) = read.predict
  colnames(EDR.curve.true) = paste("Predict", samplesize.predict, sep = "")

  #   library(rgl); open3d()
  #   plot3d(x = predict.N, y = predict.R, z = EDR.curve.true, add = F, col = 2, type = "p", xlab = "N", ylab = "R", zlab = "EDR", ylim = c(0, 1))
  #   plot3d(x = predict.N, y = predict.R, z = EDR.T, add = T, col = 1, type = "p", xlab = "N", ylab = "R", zlab = "EDR")
  mse.pilot = sapply(NR.optim.pilot, function(x) sapply(x, function(y) y$mse))
  mse.T = sqrt(mean((EDR.T-EDR.curve.true)[, -1]^2)); mse.T #drop predict 10, the worst case of curve fitting

  mse.pilot.T.fitted = sapply(NR.optim.pilot, function(x) sapply(x, function(y) sqrt(mean((y$predict-EDR.curve.true)[, -1]^2))))


  return(list(mse.pilot = mse.pilot, mse.T = mse.T, true.para = true.para, mse.pilot.T.fitted = mse.pilot.T.fitted, NR.optim.true = NR.optim.true, NR.optim.pilot = NR.optim.pilot))
})

EDR.result.Scotty = lapply(1, function(fold.change){
  # EDR.result.Storey = lapply(1:4, function(fold.change){
  # EDR.result.down = lapply(1:4, function(fold.change){
  EDR.PS.2 = lapply(1:3, function(j){ #for different pilot sample size
    load(paste("Scotty.N.", c(2, 4, 8)[j], ".rdata", sep = ""))
    EDR = EDR = get(paste("Scotty.N.", c(2, 4, 8)[j], sep = ""))
    EDR = lapply(EDR, function(x){
      x = t(x)
      x = x[1:6, -1] #drop predict 5
      read = as.numeric(rownames(x))*10000
      rownames(x) = read
      colnames(x) = paste("Predict", samplesize.predict, sep = "")
      return(x)
    })
    names(EDR) = paste("Repeat", 1:length(EDR), sep = "")
    return(EDR)
  })
  read = as.numeric(rownames(EDR.PS.2[[1]][[1]]))
  #==================================================#
  #True
  #==================================================#
  EDR.T = t(sapply(1:length(read.predict), function(r){
    load(paste("True.Dispersion.50.Depth.", r, "M.rdata", sep = ""))
    EDR.T = sapply(True.Dispersion.50[[2]], function(x) mean(x[,1]))
    EDR.T = EDR.T[-1]
    names(EDR.T) = samplesize.predict
    return(EDR.T)
  }))

  rownames(EDR.T) = read.predict
  #==================================================#
  #Surface fitting
  #==================================================#
  #pilot
  N = rep(samplesize.predict, each = length(read))
  R = rep(read, length(samplesize.predict))

  #true
  predict.N = rep(samplesize.predict, each = length(read.predict))
  predict.R = rep(read.predict, length(samplesize.predict))

  true.para = curve.fitting(n = predict.N, r = predict.R, EDR = EDR.T)

  #pilot optim
  NR.optim.pilot = lapply(EDR.PS.2, function(x){ #for each pilot data sample size
    lapply(x, function(y){ #for each repeat
      #derive parameter of inverse power law curve
      EDR.curve.para = curve.fitting(n = N, r = R, EDR = y)

      NR = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para)
      # N.optim = NR[1]; R.optim = NR[2]; EDR.optim = NR[3]; cost.optim = NR[4]
      EDR.pilot.NR.on.true = EDR.curve(NR[1], NR[2], true.para)
      EDR.curve.predict = matrix(apply(cbind(N, R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read))
      rownames(EDR.curve.predict) = read
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      mse = sqrt(mean((y-EDR.curve.predict)[, -1]^2))

      EDR.curve.predict = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], EDR.curve.para)), nrow = length(read.predict))
      rownames(EDR.curve.predict) = read.predict
      colnames(EDR.curve.predict) = paste("Predict", samplesize.predict, sep = "")

      return(list(optim = NR, predict = EDR.curve.predict, EDR.pilot.NR.on.true = EDR.pilot.NR.on.true, mse = mse, EDR.curve.para = EDR.curve.para))
    })
  })

  #true optim
  NR.optim.true = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para = true.para)

  EDR.curve.true = matrix(apply(cbind(predict.N, predict.R), 1, function(x) EDR.curve(x[1], x[2], true.para)), nrow = length(read.predict))
  rownames(EDR.curve.true) = read.predict
  colnames(EDR.curve.true) = paste("Predict", samplesize.predict, sep = "")

  #   library(rgl); open3d()
  #   plot3d(x = predict.N, y = predict.R, z = EDR.curve.true, add = F, col = 2, type = "p", xlab = "N", ylab = "R", zlab = "EDR", ylim = c(0, 1))
  #   plot3d(x = predict.N, y = predict.R, z = EDR.T, add = T, col = 1, type = "p", xlab = "N", ylab = "R", zlab = "EDR")
  mse.pilot = sapply(NR.optim.pilot, function(x) sapply(x, function(y) y$mse))
  mse.T = sqrt(mean((EDR.T-EDR.curve.true)[, -1]^2)); mse.T #drop predict 10, the worst case of curve fitting

  mse.pilot.T.fitted = sapply(NR.optim.pilot, function(x) sapply(x, function(y) sqrt(mean((y$predict-EDR.curve.true)[, -1]^2))))


  return(list(mse.pilot = mse.pilot, mse.T = mse.T, true.para = true.para, mse.pilot.T.fitted = mse.pilot.T.fitted, NR.optim.true = NR.optim.true, NR.optim.pilot = NR.optim.pilot))
})

#===================#
#MSE plot
#===================#
pdf(paste("(Dispersion", dispersion, ")MSE.pilot.true.pilot.vs.true.CDD.pilot.design.N2-16.R8-1.pdf", sep = ""), width = 7, height = 7/3*4)

EDR.result = EDR.result.pilot.design
par(mfrow = c(1, 3), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
mse.pilot = sapply(1:4, function(design) EDR.result[[design]]$mse.pilot); colnames(mse.pilot) = paste("Pilot", c(2, 4, 8, 16), sep = "")
mse.T = EDR.result[[1]]$mse.T
mse.pilot.T.fitted = sapply(1:4, function(design) EDR.result[[design]]$mse.pilot.T.fitted); colnames(mse.pilot.T.fitted) = paste("Pilot", c(2, 4, 8, 16), sep = "")

boxplot(mse.pilot, ylab = "MSE", main = "(Pilot) Goodness of fit", ylim = c(0, .15))
boxplot(mse.T, ylab = "MSE", main = "(True) Goodness of fit", ylim = c(0, .15))
boxplot(mse.pilot.T.fitted, ylab = "MSE", main = "(Pilot vs. True) Goodness of fit", ylim = c(0, .15))

dev.off()

#===================================#
#optim plot - given power, see cost
#===================================#
pdf(paste("(Dispersion", dispersion, ")Cost.optim.N.R.CDD.6.5M.pdf", sep = ""), width = 7/2*3, height = 7*2)
par(mfrow = c(4, 3), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
EDR.result = EDR.result.CDD.cost
# EDR.result = EDR.result.Storey
read.pilot = 6.5*10^6
for(j in 1:4){
  NR.optim.true = EDR.result[[j]]$NR.optim.true
  NR.optim.pilot = EDR.result[[j]]$NR.optim.pilot
  true.para = EDR.result[[j]]$true.para
  #   Cost = NR.optim.true[4]
  for(i in 1:(length(NR.optim.pilot)-1)){ #for different pilot sample size
    #     plot(0, type = "n", ylim = range(0, 12*10^6), xlim = c(0, max(samplesize.predict)), xlab = "N", ylab = "R", main = names(NR.optim.pilot)[i])
    plot(0, type = "n", ylim = range(0, 12*10^6), xlim = c(0, 60), xlab = "N", ylab = "R", main = names(NR.optim.pilot)[i])
    points(NR.optim.true[1], NR.optim.true[2], col = "purple", pch = 4, cex = 1, lwd = 5)
    for(Power in seq(0.71, 0.89, .01)) curve(power.R.function, col = "lightgray", lty = 3, add = T, xlim = c(.1, 100))

    NR = sapply(NR.optim.pilot[[i]], function(x) x$optim)
    points(NR[1,]+runif(10, -1, 0), NR[2,], col = i, lwd = 2)

    Power = power1 = 0.8; curve(power.R.function, col = 2, add = T, xlim = c(.1, 100)) #power curve
    Power = power2 = 0.7; curve(power.R.function, col = 2, lty = 2, add = T, xlim = c(.1, 100)) #power curve
    Power = power3 = 0.9; curve(power.R.function, col = 2, lty = 3, add = T, xlim = c(.1, 100)) #power curve
    #     Power = power3 = EDR.curve(10, 8*10^6, true.para); curve(power.R.function, col = "green", lty = 3, add = T, xlim = c(.1, 100)) #power curve

    Cost = cost1 = NR.optim.true[4]; curve(cost.R.function, col = 4, add = T, xlim = c(.1, 100)) #cost curve
    Cost = cost2 = NR.optim.true[4]*.8; curve(cost.R.function, col = 4, lty = 2, add = T, xlim = c(.1, 100)) #cost curve
    Cost = cost3 = NR.optim.true[4]*1.2; curve(cost.R.function, col = 4, lty = 3, add = T, xlim = c(.1, 100)) #cost curve
    points(samplesize.pilot[i], read.pilot, col = "Orange", pch = 8, cex = 1, lwd = 2)
    #   readline()
    #   segments(x0=NR[1,], y0=NR[2,], x1=NR.T[1], y1=NR.T[2])
    #     legend("topleft", c("True optimal N & R", "Pilot N & R", "Cost function"), bg = "white",
    legend("topright", c("True optimal N & R", "Pilot N & R",
                         paste("Power = ", round(power1, 2), sep = ""),
                         paste("Power = ", round(power2, 2), sep = ""),
                         paste("Power = ", round(power3, 2), sep = ""),
                         paste("Cost = ", round(cost1), sep = ""),
                         paste("Cost = ", round(cost2), sep = ""),
                         paste("Cost = ", round(cost3), sep = "")), bg = "white",
           cex = .9, pch = c(4, 8, rep(NA, 6)), col = c("purple", "orange", rep("red", 3), rep("blue", 3)), lty = c(0, 0, 1:3, 1:3), merge = T, pt.lwd = c(5, 2))
  }
}
dev.off()

#================================================#
#optim plot - given cost see power, pilot design
#================================================#
pdf(paste("(Dispersion", dispersion, ")EDR.optim.N.R.CDD.pilot.design.N2-16.R8-1.pdf", sep = ""), width = 7, height = 7)
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
EDR.result = EDR.result.pilot.design

for(j in 1:4){ #for different pilot design
  read = read.list[[j]]
  NR.optim.true = EDR.result[[j]]$NR.optim.true
  NR.optim.pilot = EDR.result[[j]]$NR.optim.pilot
  i = 1
  plot(0, type = "n", ylim = range(0, 12*10^6), xlim = c(0, 60), xlab = "N", ylab = "R", main = names(NR.optim.pilot)[i])
  points(NR.optim.true[1], NR.optim.true[2], col = "purple", pch = 4, cex = 1, lwd = 5)

  NR = sapply(NR.optim.pilot[[i]], function(x) x$optim)
  points(NR[1,]+runif(10, -1, 0), NR[2,], col = j, lwd = 2)
  curve(cost.R.function, col = 2, add = T, xlim = c(5, 100))
  points(samplesize.pilot[i], max(read), col = "Orange", pch = 8, cex = 1, lwd = 2)
  #   readline()
  #   segments(x0=NR[1,], y0=NR[2,], x1=NR.T[1], y1=NR.T[2])
  legend("topleft", c("True optimal N & R", "Pilot N & R", "Cost function"), bg = "white",
         cex = .9, pch = c(4, 8, NA), col = c("purple", "orange", "red"), lty = c(0, 0, 1), merge = T, pt.lwd = c(5, 2))
}
dev.off()

#========================================================#
#optim plot - given cost see power, scotty & RNAseqPower
#========================================================#
pdf(paste("(Dispersion", dispersion, ")EDR.optim.N.R.CDD.RNAseqP.Scotty.pdf", sep = ""), width = 7/2*3, height = 7/2*3)
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
EDR.result.list = list(RNAseqP = EDR.result.RNAseqP, Scotty = EDR.result.Scotty, SeqDesign = EDR.result.CDD[2]) #fold change 2
read.max = c(.6*10^7, .6*10^7, 6.5*10^6)
method = c("RNAseqPower", "Scotty", "SeqDesign")
for(j in 1:3){ #different method
  EDR.result = EDR.result.list[[j]]
  NR.optim.true = EDR.result[[1]]$NR.optim.true
  NR.optim.pilot = EDR.result[[1]]$NR.optim.pilot
  true.para = EDR.result[[1]]$true.para
  for(i in 1:3){ #different samplesize
    plot(0, type = "n", ylim = range(0, 12*10^6), xlim = c(0, 60), xlab = "N", ylab = "R", main = paste(method[j], "Pilot", c(2, 4, 8)[i]))
    points(NR.optim.true[1], NR.optim.true[2], col = "purple", pch = 4, cex = 1, lwd = 5)

    NR = sapply(NR.optim.pilot[[i]], function(x) x$optim)
    points(NR[1,]+runif(10, -1, 0), NR[2,], col = i, lwd = 2)
    curve(cost.R.function, col = "blue", add = T, xlim = c(5, 100))

    Power = power1 = 0.8; curve(power.R.function, col = 2, add = T, xlim = c(.1, 100)) #power curve
    Power = power2 = 0.7; curve(power.R.function, col = 2, lty = 2, add = T, xlim = c(.1, 100)) #power curve
    Power = power3 = 0.9; curve(power.R.function, col = 2, lty = 3, add = T, xlim = c(.1, 100)) #power curve
    #     plot(0, type = "n", xlim = c(0, 60), ylim = c(0, 1.2*10^7))
    points(samplesize.pilot[i], read.max[j], col = "Orange", pch = 8, cex = 1, lwd = 2)
    legend("topright", c("True optimal N & R", "Pilot N & R",
                         paste("Power = ", round(power1, 2), sep = ""),
                         paste("Power = ", round(power2, 2), sep = ""),
                         paste("Power = ", round(power3, 2), sep = ""),
                         paste("Cost = ", round(Cost), sep = "")),
           cex = .9, pch = c(4, 8, rep(NA, 4)), col = c("purple", "orange", rep("red", 3), rep("blue", 1)), lty = c(0, 0, 1:3, 1), merge = T, pt.lwd = c(5, 2))
  }
}
dev.off()


#=========================================================#
#                         Boxplot
#=========================================================#
pdf(paste("(Dispersion", dispersion, ")Cost.optim.N.R.CDD.6.5M.pilot.design.N2-16.R8-1.boxplot.pdf", sep = ""), width = 7, height = 7)
par(mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

EDR.pilot = lapply(1:4, function(design){
  NR.optim.true = EDR.result[[design]]$NR.optim.true
  NR.optim.pilot = EDR.result[[design]]$NR.optim.pilot
  EDR.pilot = sapply(NR.optim.pilot, function(x) sapply(x, function(y) y$EDR.pilot.NR.on.true))

})

names(EDR.pilot) = paste("Pilot", c(2, 4, 8, 16), sep = "")
EDR.true = EDR.result[[1]]$NR.optim.true[3]
boxplot(EDR.pilot, ylab = "EDR", main = "Compare pilot design", ylim = c(0.5, 1)); abline(h = EDR.true, col = 2, lwd = 2)

dev.off()

#==========================================================#
for(fold.change in 1:4){
  fold.change = 1
  print(fold.change)
  EDR.PS.2 = lapply(samplesize.pilot, function(j){
    EDR = lapply(1:n.repeat, function(i){ #for each repeat
      load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.RData",sep=""))
      #       load(paste("(EDR)New.setting.Pilot.repeat.", i, ".samplesize.", j, ".predict.5-100.read.1M-6M.M10.2way.Fold.", fold.change, ".dispersion.", dispersion, ".CDD.new.RData",sep=""))
      EDR = t(matrix(apply(sapply(Pilot.power$Result, function(x) sapply(x, function(y) y[3,])), 1, mean), nrow = length(samplesize.predict)))
      rownames(EDR) = read; colnames(EDR) = paste("Predict", samplesize.predict, sep = "")
      return(EDR)
    })
    names(EDR) = paste("Repeat", 1:n.repeat, sep = "")
    return(EDR)
  })
  names(EDR.PS.2) = paste("Pilot", samplesize.pilot, sep = "")

  if(0){
    EDR.T = t(sapply(length(read.predict):1, function(r){
      load(paste("True.result.", r, "SamplePerLane.rdata", sep = ""))
      a = get(paste("True.result.", r, "SamplePerLane", sep = ""))
      EDR.T = sapply(a, function(x) mean(x[,9]))
      #     EDR.T = EDR.T[-1]
      names(EDR.T) = samplesize.predict
      return(EDR.T)
    }))
  }

  # New setting -------------------------------------------------------------
  load("EDR.True.Dispersion.50.N.4.sd.0.10.2sample.perLane.Result.rdata")
  EDR.T = sapply(EDR.True.Dispersion.50.N.4.sd.0.10.2sample.perLane.Result, function(x){ #for each repeat
    sapply(x, function(y){ #for each sample size
      y[1,]
    })
  })
  EDR.T = matrix(apply(EDR.T, 1, mean, na.rm = T), nrow = 10)
  rownames(EDR.T) = read.predict; colnames(EDR.T) = names(EDR.True.Dispersion.50.N.4.sd.0.10.2sample.perLane.Result[[1]])

  #==================================================#
  #Surface fitting
  #==================================================#
  #pilot
  N = rep(samplesize.predict, each = length(read))
  R = rep(read, length(samplesize.predict))

  #true
  predict.N = rep(samplesize.predict, each = length(read.predict))
  predict.R = rep(read.predict, length(samplesize.predict))

  true.para = curve.fitting(n = predict.N, r = predict.R, EDR = EDR.T)

  # Cost.seq = seq(10^4, 10^6, 10^4)
  Cost.seq = c(seq(10^4, 5*10^4, length = 20), seq(5*10^4, 2*10^5, 10^4))

  library(snowfall)
  sfInit(cpus = 3, parallel = T)
  sfExportAll(); sfLibrary(DEoptim)
  optim.pilot = sfLapply(Cost.seq, function(Cost){
    #==================================================#
    #pilot optim
    if(0){
      NR.optim.pilot = lapply(EDR.PS.2, function(x){ #for each pilot data sample size
        sapply(x, function(y){ #for each repeat
          #derive parameter of inverse power law curve
          EDR.curve.para = curve.fitting(n = N, r = R, EDR = y)

          NR = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para)

          return(c(optim = NR))
        })
      })
    }
    #true optim
    NR.optim.true = Optim.Function(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para = true.para)

    #     return(list(NR.optim.true = NR.optim.true, NR.optim.pilot = NR.optim.pilot))
    return(list(NR.optim.true = NR.optim.true, NR.optim.pilot = NULL))
  })

  #   N.restrict.vec = c(20, 30, 40, 60, 80)

  #Under restriction
  N.restrict.vec = c(40)
  optim.true = lapply(N.restrict.vec, function(N.restrict){
    print(N.restrict)
    #   sfExportAll(); sfLibrary(DEoptim)
    optim.pilot.restrict = sfLapply(Cost.seq, function(Cost){
      #==================================================#
      #pilot optim
      if(0){
        NR.optim.pilot = lapply(EDR.PS.2, function(x){ #for each pilot data sample size
          sapply(x, function(y){ #for each repeat
            #derive parameter of inverse power law curve
            EDR.curve.para = curve.fitting(n = N, r = R, EDR = y)

            NR = Optim.Function.restrict.N(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para, N.restrict = N.restrict)

            return(c(optim = NR))
          })
        })
      }
      #true optim
      NR.optim.true = Optim.Function.restrict.N(a = Cost.N, b = Cost.R, Cost = Cost, EDR.curve.para = true.para, N.restrict = N.restrict)

      #   return(list(NR.optim.true = NR.optim.true, NR.optim.pilot = NR.optim.pilot))
      return(list(NR.optim.true = NR.optim.true, N.restrict = N.restrict))
    })
  })
  # NR.pilot = lapply(optim.pilot, function(x) x$NR.optim.pilot); NR.pilot.restrict = lapply(optim.pilot.restrict, function(x) x$NR.optim.pilot)
  # NR.true = lapply(optim.pilot, function(x) x$NR.optim.true); NR.true.restrict = lapply(optim.pilot.restrict, function(x) x$NR.optim.true);
  NR.true = sapply(optim.pilot, function(x) x$NR.optim.true); NR.true.restrict = lapply(optim.true, function(y) sapply(y, function(x) x$NR.optim.true))
  colnames(NR.true) = paste("Cost = ", Cost.seq, sep = ""); names(NR.true.restrict) = N.restrict.vec

  if(0){
    Cost.EDR.pilot = lapply(1:4, function(i){ #for each pilot sample size
      Cost.EDR = sapply(NR.pilot, function(x){ #for each cost
        x[[i]][3,] #for each pilot sample size
      })
      colnames(Cost.EDR) = Cost.seq
      return(Cost.EDR)
    })
    Cost.EDR.pilot.restrict = lapply(1:4, function(i){ #for each pilot sample size
      Cost.EDR = sapply(NR.pilot.restrict, function(x){ #for each cost
        x[[i]][3,] #for each pilot sample size
      })
      colnames(Cost.EDR) = Cost.seq
      return(Cost.EDR)
    })
    names(Cost.EDR.pilot) = names(Cost.EDR.pilot.restrict) = names(NR.pilot[[1]])

    Cost.N.pilot = lapply(1:4, function(i){ #for each pilot sample size
      Cost.EDR = sapply(NR.pilot, function(x){ #for each cost
        x[[i]][1,] #for each pilot sample size
      })
      colnames(Cost.EDR) = Cost.seq
      return(Cost.EDR)
    })
    Cost.N.pilot.restrict = lapply(1:4, function(i){ #for each pilot sample size
      Cost.EDR = sapply(NR.pilot.restrict, function(x){ #for each cost
        x[[i]][1,] #for each pilot sample size
      })
      colnames(Cost.EDR) = Cost.seq
      return(Cost.EDR)
    })

    names(Cost.EDR.pilot) = names(Cost.N.pilot.restrict) = names(NR.pilot[[1]])

    Cost.R.pilot = lapply(1:4, function(i){ #for each pilot sample size
      Cost.EDR = sapply(NR.pilot, function(x){ #for each cost
        x[[i]][2,] #for each pilot sample size
      })
      colnames(Cost.EDR) = Cost.seq
      return(Cost.EDR)
    }); Cost.R.pilot.restrict = lapply(1:4, function(i){ #for each pilot sample size
      Cost.EDR = sapply(NR.pilot.restrict, function(x){ #for each cost
        x[[i]][2,] #for each pilot sample size
      })
      colnames(Cost.EDR) = Cost.seq
      return(Cost.EDR)
    })

    names(Cost.R.pilot) = names(Cost.R.pilot.restrict) = names(NR.pilot[[1]])
  }

  Cost.EDR.true = NR.true[3,]
  Cost.EDR.true.restrict = sapply(NR.true.restrict, function(x) x[3,])
  rownames(Cost.EDR.true.restrict) = Cost.seq

  #   index = c(10, 20, 30)
  # abline(v = Cost.seq[index], col = 2)

  #N&R
  # N.true = sapply(NR.true, function(x) x[1]); N.true.restrict = sapply(NR.true.restrict, function(x) x[1])
  N.true = NR.true[1,]; N.true.restrict = sapply(NR.true.restrict, function(x) x[1,]); rownames(N.true.restrict) = Cost.seq
  # R.true = sapply(NR.true, function(x) x[2]); R.true.restrict = sapply(NR.true.restrict, function(x) x[2])
  R.true = NR.true[2,]; R.true.restrict = sapply(NR.true.restrict, function(x) x[2,]); rownames(R.true.restrict) = Cost.seq


  if(0){
    # pdf("try.pdf", width = 7, height = 7/2)
    # pdf(paste("try", fold.change, "_2.pdf", sep = ""), width = 7, height = 7/2*4)
    pdf(paste("try", fold.change, "_2.pdf", sep = ""), width = 7, height = 7/2*1)
    # par(mfrow = c(4, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
    par(mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 1), mgp = c(1.5, .5, 0))

    # for(j in 1:4){ #for different pilot samplesize
    plot(Cost.seq, Cost.EDR.true, col = 1, type = "l", ylim = c(0.75, 1.05), xlab = "Given Cost", ylab = "optimial EDR",
         #        main = paste("Pilot N = ", c(2, 4, 8, 16)[j], sep = ""), lwd = 2)
         main = "", lwd = 2, cex.axis = .75)
    for(i in 1:ncol(Cost.EDR.true.restrict)) lines(Cost.seq, Cost.EDR.true.restrict[,i], lty = 2, col = i+1, lwd = 2)
    #   for(i in 1:nrow(Cost.EDR.pilot[[j]])) lines(Cost.seq, Cost.EDR.pilot[[j]][i,], lty = 2)
    #   for(i in 1:nrow(Cost.EDR.pilot.restrict[[j]])) lines(Cost.seq, Cost.EDR.pilot.restrict[[j]][i,], lty = 2, col = "gray")
    points(Cost.seq[index], Cost.EDR.true[index], pch = c(8, 9, 11), col = 1, lwd = 2, cex = 1)
    for(i in 1:ncol(Cost.EDR.true.restrict)) points(Cost.seq[index], Cost.EDR.true.restrict[index, i], pch = c(8, 9, 11), col = 1+i, lwd = 2, cex = 1)
    legend("topleft", c("No constraint", paste("N <= ", N.restrict.vec, sep = "")), cex = .5, col = 1:6, lty = c(1, rep(2, 5)))

    # plot(N.true[index], R.true[index], pch = c(8, 9, 11), col = c("blue", "green", "purple"), lwd = 3, xlab = "optimal N", ylab = "optimal R", ylim = c(0, 5*10^6))
    plot(N.true[index], R.true[index], pch = c(8, 9, 11), col = 1, lwd = 2, cex = 1,
         xlab = "optimal N", ylab = "optimal R", cex.axis = .75,
         #        main = paste("Pilot N = ", c(2, 4, 8, 16)[j], sep = ""),
         main = "",
         #        xlim = c(0, 200), ylim = c(0, 5*10^7))
         xlim = c(0, 150), ylim = c(0, 2*10^8))
    #   points(Cost.N.pilot[[j]][, index] + runif(30, -1, 1), Cost.R.pilot[[j]][, index], pch = rep(c(8, 9, 11), each = 10))

    for(i in 1:ncol(Cost.EDR.true.restrict)) points(N.true.restrict[index, i], R.true.restrict[index, i], pch = c(8, 9, 11), col = i+1, lwd = 2, cex = 1)
    #   points(Cost.N.pilot.restrict[[j]][, index] + runif(30, -1, 1), Cost.R.pilot.restrict[[j]][, index], pch = rep(c(15:17), each = 10))
    # }
    dev.off()
  }
} # end of fold.change

#==================================================#
#draw each biological question one by one
#==================================================#
#================================#
#               Q1
#================================#
pdf(paste("Q1.pdf", sep = ""), width = 7, height = 7/2*1)
# par(mfrow = c(4, 2), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 1), mgp = c(1.5, .5, 0))

index = which(names(N.true) == "Cost = 1e+05")
# for(j in 1:4){ #for different pilot samplesize
plot(Cost.seq, Cost.EDR.true, col = 1, type = "l", ylim = c(.3, 1), xlab = "Optimal Cost", ylab = "Given EDR",
     main = "Q1", lwd = 2, cex.axis = .75)
# abline(h = 0.95, lty = 2, lwd = 2, col = 2);
abline(h = Cost.EDR.true[index], lty = 2, lwd = 2, col = 2);
# points(45000, 0.95, pch = 8, col = 2, lwd = 2, cex = 1)
points(100000, Cost.EDR.true[index], pch = 8, col = 2, lwd = 2, cex = 1)


# plot(N.true[index], R.true[index], pch = c(8, 9, 11), col = c("blue", "green", "purple"), lwd = 3, xlab = "optimal N", ylab = "optimal R", ylim = c(0, 5*10^6))
plot(N.true[index], R.true[index], pch = 8, col = 1, lwd = 2, cex = 1,
     xlab = "optimal N", ylab = "optimal R", cex.axis = .75,
     main = "",
     xlim = c(0, 200), ylim = c(0, 10^8)
)
dev.off()

#================================#
#               Q2
#================================#
pdf(paste("Q2.pdf", sep = ""), width = 7, height = 7/2*1)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 1), mgp = c(1.5, .5, 0))
index = 31
# for(j in 1:4){ #for different pilot samplesize
plot(Cost.seq, Cost.EDR.true, col = 1, type = "l", xlab = "Given Cost", ylab = "Optimal EDR",
     main = "Q2", lwd = 2, cex.axis = .75, ylim = c(.3, 1))
# abline(v = 45000, lty = 2, lwd = 2, col = 2);
# points(45000, 0.95, pch = 8, col = 2, lwd = 2, cex = 1)
abline(v = Cost.seq[index], lty = 2, lwd = 2, col = 2);
points(Cost.seq[index], Cost.EDR.true[index], pch = 8, col = 2, lwd = 2, cex = 1)

# plot(N.true[index], R.true[index], pch = c(8, 9, 11), col = c("blue", "green", "purple"), lwd = 3, xlab = "optimal N", ylab = "optimal R", ylim = c(0, 5*10^6))
plot(N.true[index], R.true[index], pch = 8, col = 1, lwd = 2, cex = 1,
     xlab = "optimal N", ylab = "optimal R", cex.axis = .75,
     main = "",
     xlim = c(0, 200), ylim = c(0, 10^8)
)
dev.off()

#================================#
#               Q3
#================================#
pdf(paste("Q3.pdf", sep = ""), width = 7, height = 7/2*1)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 1), mgp = c(1.5, .5, 0))
index = 26
# for(j in 1:4){ #for different pilot samplesize
plot(Cost.seq, Cost.EDR.true.restrict[,1], col = 1, type = "l", ylim = c(.3, 1), xlab = "Optimal Cost", ylab = "Given EDR",
     main = "Q3", lwd = 2, cex.axis = .75)
abline(h = Cost.EDR.true.restrict[index,1], lty = 2, lwd = 2, col = 2);
points(Cost.seq[index], Cost.EDR.true.restrict[index,1], pch = 8, col = 2, lwd = 2, cex = 1)

# plot(N.true[index], R.true[index], pch = c(8, 9, 11), col = c("blue", "green", "purple"), lwd = 3, xlab = "optimal N", ylab = "optimal R", ylim = c(0, 5*10^6))
plot(N.true.restrict[index], R.true.restrict[index], pch = 8, col = 1, lwd = 2, cex = 1,
     xlab = "optimal N", ylab = "optimal R", cex.axis = .75,
     main = "", type = "n",
     xlim = c(0, 200), ylim = c(0, 10^8))
abline(v = 40, lty = 2, col = 2, lwd = 2); polygon(x = c(0, 0, 40, 40), y = c(3*10^8, 0, 0, 3*10^8), density = 10, border = NA, col = 2)
points(N.true.restrict[index], R.true.restrict[index], pch = 8, col = 1, lwd = 2, cex = 1)
dev.off()
# plot(1:10, 1:10)
# polygon(x = c(1, 1, 6, 6), y = c(8, 2, 2, 10), density = 20, border = NA)
#================================#
#               Q4
#================================#
pdf(paste("Q4.pdf", sep = ""), width = 7, height = 7/2*1)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 1, 1), mgp = c(1.5, .5, 0))
index = 23
# for(j in 1:4){ #for different pilot samplesize
#restrict
plot(Cost.seq, Cost.EDR.true.restrict[,1], col = 4, type = "l", ylim = c(.3, 1), xlab = "Optimal Cost", ylab = "Given EDR",
     main = "Q4", lwd = 2, cex.axis = .75)
#true
lines(Cost.seq, Cost.EDR.true, col = 1, lty = 2, lwd = 2)

# abline(h = Cost.EDR.true.restrict[index], lty = 2, lwd = 2, col = 2);
# points(Cost.seq[index], Cost.EDR.true[index], pch = 1, col = 2, lwd = 2, cex = 1.5) #true
# points(75000, Cost.EDR.true[index], pch = 2, col = 2, lwd = 2, cex = 1.5) #given power, restrict
# points(75000, 0.939, pch = 8, col = 2, lwd = 2, cex = 1.5) #given cost, restrict

# plot(N.true[index], R.true[index], pch = c(8, 9, 11), col = c("blue", "green", "purple"), lwd = 3, xlab = "optimal N", ylab = "optimal R", ylim = c(0, 5*10^6))
if(0){
  plot(40, 22500000, pch = 1, col = 2, lwd = 2, cex = 1,
       xlab = "optimal N", ylab = "optimal R", cex.axis = .75,
       main = "",
       xlim = c(0, 150), ylim = c(0, 2*10^8))
  points(50, 5000000, pch = 8, col = 2)
}
dev.off()


