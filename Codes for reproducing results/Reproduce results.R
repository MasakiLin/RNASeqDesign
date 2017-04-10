source("functions.R")
setwd("Data/")
#---------------------------------------
# 1. Illustration of T1-T5
#---------------------------------------
load("Two dimentsion hyperplane parameter.RData")

Cost.N = 500*2 #cost per sample
Cost.R = (25/10^6)*2 #cost per read

N.discrete = seq(50, 200, 10)
R.discrete = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid = rep(N.discrete, each = length(R.discrete))
R.grid = rep(R.discrete, times = length(N.discrete))

Cost.discrete = apply(cbind(N.grid, R.grid), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete = apply(cbind(N.grid, R.grid), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete[EDR.discrete < 0] = 0

index.original = find.admissible(Cost.discrete, EDR.discrete)

pdf("Figure 2.pdf", width = 15, height = 8)
par(mfcol = c(2, 5), mgp = c(2.8, 1.2, 0), mar = c(4, 5, 3, 1) + .1, cex.axis = 1.8, cex.lab = 2, cex.main = 3)
# --------------------
# T1
# --------------------
Cost = 200000
T1(Cost.discrete, EDR.discrete, index.original, N.grid, R.grid, R.discrete, Cost)
# --------------------
# T2
# --------------------
EDR = 0.85

T2(Cost.discrete, EDR.discrete, index.original, N.grid, R.grid, R.discrete, EDR)
# --------------------
# T3
# --------------------
EDR = 0.85

N.discrete.res = seq(50, 80, 10)
R.discrete.res = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid.res = rep(N.discrete.res, each = length(R.discrete.res))
R.grid.res = rep(R.discrete.res, times = length(N.discrete.res))

Cost.discrete.res = apply(cbind(N.grid.res, R.grid.res), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res = apply(cbind(N.grid.res, R.grid.res), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res[EDR.discrete.res < 0] = 0

index.original.res = find.admissible(Cost.discrete.res, EDR.discrete.res)

T3(Cost.discrete.res, EDR.discrete.res, index.original.res, N.grid.res, R.grid.res, R.discrete, EDR)

# --------------------
# T4
# --------------------
N.discrete.res.130 = seq(50, 130, 10)
R.discrete.res.130 = 6*10^7/c(4, 2, 1, 2/3, 1/2)
N.grid.res.130 = rep(N.discrete.res.130, each = length(R.discrete.res.130))
R.grid.res.130 = rep(R.discrete.res.130, times = length(N.discrete.res.130))

Cost.discrete.res.130 = apply(cbind(N.grid.res.130, R.grid.res.130), 1, function(x) Cost.N.R(x[1], x[2]))
EDR.discrete.res.130 = apply(cbind(N.grid.res.130, R.grid.res.130), 1, function(x) EDR.curve(x[1], x[2], true.para))
EDR.discrete.res.130[EDR.discrete.res.130 < 0] = 0

index.original.res.130 = find.admissible(Cost.discrete.res.130, EDR.discrete.res.130)

T4(Cost.discrete.res, EDR.discrete.res, index.original.res, Cost.discrete.res.130, EDR.discrete.res.130, index.original.res.130,
   N.grid.res, R.grid.res, N.grid.res.130, R.grid.res.130, R.discrete, EDR)

# --------------------
# T5
# --------------------
T5(Cost.discrete, EDR.discrete, R.grid, R.discrete)
dev.off()

# ---------------------------
# 2. Simulation study
# ---------------------------
# ----------------------------
# 1.1. simulate true data
# ----------------------------
load("Parameter setting estimated from HIP data 20M.RData")
load("Normalized.HIP.data.rdata")

par.setting$lfc.mean.sd = c(0, 0.2)

par.setting$Prop.DE = 0.15
par.setting$lfc = c(-log2(1.4), log2(1.4))

num.repeat = 20
sample.size = c(2, 4, 8, 16)
target.N = c(5, 10, 20, 30, 40, 50, 100)
library(snowfall)

for(dispersion in c(2, 5, 10, 20)){

  par.setting$dispersion = dispersion

  sfInit(parallel = TRUE, cpus = 20)

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
                               FDR = 0.05, filter = T, filter.level = 5)
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

  library(snowfall)
  sfInit(parallel = TRUE, cpus = 20)

  num.repeat = 20
  sample.size = c(2, 4, 8, 16)
  target.N = c(5, 10, 20, 30, 40, 50, 100)

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

}
sfStop()

#predict EDR
n.repeat = 20

for(dispersion in c(2, 5, 10, 20)){
  print(paste("Dispersion = ", dispersion, sep = ""))

  load(paste("Simulated.data.20M.20.repeat.dispersion.", dispersion, ".lfc.mean.0.sd.0.2.N.2.rdata")
  load(paste("Simulated.data.20M.20.repeat.dispersion.", dispersion, ".lfc.mean.0.sd.0.2.N.4.rdata")
  load(paste("Simulated.data.20M.20.repeat.dispersion.", dispersion, ".lfc.mean.0.sd.0.2.N.8.rdata")
  load(paste("Simulated.data.20M.20.repeat.dispersion.", dispersion, ".lfc.mean.0.sd.0.2.N.16.rdata")

  library(snowfall)
  sfInit(parallel = TRUE, cpus = 20)
  set.seed(12345)

  sfExportAll()

  for(i in 1:length(sample.size)){ #16
    print(i)
    t1=Sys.time()
    sfExport("i")
    assign(paste("RNASeqDesign.", sample.size[i], sep = ""),
           sfLapply(1:n.repeat,function(x){
             Data = get(paste("simulation.N.", sample.size[i], sep = ""))[[x]]
             n = ncol(Data)/2
             status = c(rep("Case", n), rep("Control", n))
             Estimate.EDR.from.pilot(Data, status, group.name = c("Control", "Case"),
                                     FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                     target.N = target.N)
           }))
    t2=Sys.time()
    print(t2-t1)
    do.call(save, list(paste("RNASeqDesign.", sample.size[i], sep = ""),
                       file = paste("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.", dispersion, ".lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.",
                                    sample.size[i], ".predict.EDR.not.considering.lfc.rdata", sep = "")))
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
sfExportAll()

for(i in 1:4){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Poisson.", i, sep = ""),
         sfSapply(1:length(Data.2),function(x){
           Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           Poisson(Data, status, target.N = target.N, rho = 1.4, DE.prop = .15)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Poisson.", i, sep = ""),
                     file = paste("(Simulation)Poisson.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:4){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("RNASeqPower.", i, sep = ""),
         sfSapply(1:length(Data.2),function(x){
           Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           RNASeqPower(Data, status, target.N = target.N, n.prop = 1,
                       effect = 1.4, alpha = 0.0001)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("RNASeqPower.", i, sep = ""),
                     file = paste("(Simulation)RNASeqPower.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
}

#NB
for(i in 1:4){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("NB.", i, sep = ""),
         sfSapply(1:length(Data.2),function(x){
           Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           NB.exact(Data, status, target.N = target.N, rho = 1.4, DE.prop = .15)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("NB.", i, sep = ""),
                     file = paste("(Simulation)NB.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".rdata", sep = "")))
}

#Scotty, prepare files for running on Matlab
for(i in 1:4){
  for(j in 1:length(Data.2)){ #repeatment
    print(paste("i", i, "j", j))

    Data = get(paste("Data.", sample.size[i], sep = ""))[[j]]
    n = ncol(Data)/2
    status = c(rep("Case", n), rep("Control", n))

    colnames(Data) = c(paste("Case_", 1:n, sep = ""), paste("Control_", 1:n, sep = ""))

    Data = Data[apply(Data, 1, mean) > 5,]

    Data = rbind(c("Gene", colnames(Data)), cbind(row.names(Data), Data))

    write.table(Data, file = paste("Scotty.dispersion.", round(dispersion, 2), ".N.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#PROPER
for(i in 1:4){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("PROPER.", i, sep = ""),
         sfSapply(1:length(Data.2),function(x){
           Data = get(paste("Data.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           PROPER(Data, status, n.prop = 1, target.N = target.N,
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
for(i in 1:length(sample.size)) load(paste("(Simulation)True.EDR.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.", sample.size[i], ".not.considering.lfc.rdata", sep = ""))
True = sapply(1:length(sample.size), function(i){
  x = get(paste("True.", sample.size[i], sep = ""))
  x = x[sapply(x, length) != 1]
  mean(sapply(x, function(y) y$Result[3]))
}); True

# RNASeqDesign-----------
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.2.predict.EDR.not.considering.lfc.rdata")
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.4.predict.EDR.not.considering.lfc.rdata")
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.8.predict.EDR.not.considering.lfc.rdata")

for(i in c(2, 4, 8)) assign(paste("p.", i, sep = ""),
                            Each.EDR(get(paste("RNASeqDesign.", i, sep = "")), method = "RNASeqDesign",
                                     True = True, target.N = sample.size,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))
# Poisson------------
load("(Simulation)Poisson.dispersion.5.N.2.rdata")
load("(Simulation)Poisson.dispersion.5.N.4.rdata")
load("(Simulation)Poisson.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.Poisson.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("Poisson.", i, sep = "")), target.N = sample.size,
                                     method = "Poisson", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# RNASeqPower--------
load("(Simulation)RNASeqPower.dispersion.5.N.2.rdata")
load("(Simulation)RNASeqPower.dispersion.5.N.4.rdata")
load("(Simulation)RNASeqPower.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.RNASeqPower.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("RNASeqPower.", i, sep = "")), target.N = sample.size,
                                     method = "RNASeqPower", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# NB------------
load("(Simulation)NB.dispersion.5.N.2.rdata")
load("(Simulation)NB.dispersion.5.N.4.rdata")
load("(Simulation)NB.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.NB.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("NB.", i, sep = "")), target.N = sample.size,
                                     method = "NB", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# PROPER------------
load("(Simulation)PROPER.default.dispersion.5.N.2.rdata")
load("(Simulation)PROPER.default.dispersion.5.N.4.rdata")
load("(Simulation)PROPER.default.dispersion.5.N.8.rdata")
for(i in c(1, 2, 3)) assign(paste("p.PROPER.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("PROPER.", i, sep = "")), target.N = sample.size,
                                     method = "PROPER", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))
# Scotty----------
for(i in c(1, 2, 3)){
  tmp=read.table(paste("Simulation_Scotty_dispersion_5_N_", pilot.size[i], ".txt", sep = ""), sep = ",")
  assign(paste("Scotty.", i, sep = ""), {
    x = sapply(1:10, function(j) tmp[(1+149*(j-1)):(149*(j-1)+149),][c(sample.size)-1,10]/100)
  })
}

for(i in c(1, 2, 3)) assign(paste("p.Scotty.", pilot.size[i], sep = ""),
                            Each.EDR(get(paste("Scotty.", i, sep = "")), target.N = sample.size,
                                     method = "Scotty", True = True,
                                     variation.True = F, sd.True = NULL, output.MSE = T, Power.max = 0.8))

# plot-----------------------------------------------------------------------------------------------
pdf("Figure 3.pdf", width = 18, height = 9)
multiplot(p.Poisson.2[[1]], p.RNASeqPower.2[[1]], p.NB.2[[1]], p.Scotty.2[[1]], p.PROPER.2[[1]], p.2[[1]],
          p.Poisson.4[[1]], p.RNASeqPower.4[[1]], p.NB.4[[1]], p.Scotty.4[[1]], p.PROPER.4[[1]], p.4[[1]],
          p.Poisson.8[[1]], p.RNASeqPower.8[[1]], p.NB.8[[1]], p.Scotty.8[[1]], p.PROPER.8[[1]], p.8[[1]],

          cols=6,layout=matrix(1:18,ncol=6,byrow=T))

dev.off()

# Benchmark 1
MSE.EDR=matrix(c(p.Poisson.2[[2]],p.Poisson.4[[2]],p.Poisson.8[[2]],
                 p.RNASeqPower.2[[2]],p.RNASeqPower.4[[2]],p.RNASeqPower.8[[2]],
                 p.NB.2[[2]],p.NB.4[[2]],p.NB.8[[2]],
                 p.Scotty.2[[2]],p.Scotty.4[[2]],p.Scotty.8[[2]],
                 p.PROPER.2[[2]],p.PROPER.4[[2]],p.PROPER.8[[2]],
                 p.2[[2]], p.4[[2]], p.8[[2]]),nrow=6,byrow=T)

# Benchmark 2
MSE.N=matrix(c(p.Poisson.2[[3]],p.Poisson.4[[3]],p.Poisson.8[[3]],
               p.RNASeqPower.2[[3]],p.RNASeqPower.4[[3]],p.RNASeqPower.8[[3]],
               p.NB.2[[3]],p.NB.4[[3]],p.NB.8[[3]],
               p.Scotty.2[[3]],p.Scotty.4[[3]],p.Scotty.8[[3]],
               p.PROPER.2[[3]],p.PROPER.4[[3]],p.PROPER.8[[3]],
               p.2[[3]], p.4[[3]], p.8[[3]]),nrow=6,byrow=T)

#----------------------------
# 2. Real data applications
#----------------------------
#---------------------------------------
# 2.1. HIP mice data
#---------------------------------------
load("HIP.n.2.data.rdata")
load("HIP.n.4.data.rdata")
load("HIP.n.6.data.rdata")
load("HIP.n.8.data.rdata")
load("HIP.n.10.data.rdata")
load("HIP.n.12.data.rdata")

sample.size = c(2, 4, 6, 8, 10, 12)

num.resample = 10

#RNARNASeqDesign
library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
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
                                   target.N = c(n:12,20,30,40,50,100))
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("HIP.n.", sample.size[i], sep = ""),
                     file = paste("(HIP)RNASeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.", sample.size[i], ".predict.EDR.not.considering.lfc.rdata", sep = "")))
}

i=6
t1=Sys.time()
Data = get(paste("HIP.n.", sample.size[i], sep = ""))
n = ncol(Data)/2
status = c(rep("Case", n), rep("Control", n))
HIP.n.12 = Estimate.EDR.from.pilot(Data, status, group.name = c("Control", "Case"),
                                   FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                   target.N = c(n:12,20,30,40,50,100))
t2=Sys.time()
print(t2-t1)
do.call(save, list(paste("HIP.n.", sample.size[i], sep = ""),
                   file = paste("(HIP)RNASeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.", sample.size[i], ".predict.EDR.not.considering.lfc.rdata", sep = "")))

#other methods
load("HIP.n.2.data.rdata")
load("HIP.n.4.data.rdata")
load("HIP.n.6.data.rdata")
load("HIP.n.8.data.rdata")
load("HIP.n.10.data.rdata")
load("HIP.n.12.data.rdata")
HIP.n.12 = list(HIP.n.12)
#----------------------------------------------------
sample.size = c(2, 4, 6, 8, 10, 12)

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfExportAll()

#Poisson
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Poisson.", i, sep = ""),
         sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           Poisson(Data, status, target.N = c(n:12,20,30,40,50,100), rho = 1.4, DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Poisson.", i, sep = ""),
                     file = paste("(HIP)Poisson.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("RNASeqPower.", i, sep = ""),
         sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           RNASeqPower(Data, status, target.N = c(n:12,20,30,40,50,100), n.prop = 1,
                       effect = 1.4, alpha = 0.0001)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("RNASeqPower.", i, sep = ""),
                     file = paste("(HIP)RNASeqPower.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
}

#NB
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("NB.", i, sep = ""),
         sfSapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           NB.exact(Data, status, target.N = c(n:12,20,30,40,50,100), rho = 1.4, DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("NB.", i, sep = ""),
                     file = paste("(HIP)NB.fold.change.1.4.N.", sample.size[i], ".rdata", sep = "")))
}

#Scotty, on Matlab
for(i in 1:6){
  for(j in 1:length(get(paste("HIP.n.", sample.size[i], sep = "")))){ #repeatment
    print(paste("i", i, "j", j))

    Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[j]]
    n = ncol(Data)/2
    status = c(rep("Case", n), rep("Control", n))

    colnames(Data) = c(paste("Case_", 1:n, sep = ""), paste("Control_", 1:n, sep = ""))

    Data = Data[apply(Data, 1, mean) > 5,]

    Data = rbind(c("Gene", colnames(Data)), cbind(row.names(Data), Data))

    write.table(Data, file = paste("Scotty.fold.change.1.4.N.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#PROPER
for(i in 1:5){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("PROPER.", i, sep = ""),
         sfSapply(1:num.repeat, function(x){
           Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))

           PROPER(Data, status, n.prop = 1, target.N = c(n:12,20,30,40,50,100),
                  effect = log(1.4), DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("PROPER.", i, sep = ""),
                     file = paste("(HIP)PROPER.fold.change.1.4.N.", sample.size[i], ".default.rdata", sep = "")))
}

i=6
print(i)
t1=Sys.time()

assign(paste("PROPER.", i, sep = ""),
       sapply(1:length(get(paste("HIP.n.", sample.size[i], sep = ""))),function(x){
         Data = get(paste("HIP.n.", sample.size[i], sep = ""))[[x]]
         n = ncol(Data)/2
         status = c(rep("Case", n), rep("Control", n))

         PROPER(Data, status, n.prop = 1, target.N = c(n:12,20,30,40,50,100),
                effect = log(1.4), DE.prop = 0.1)
       }))
t2=Sys.time()
print(t2-t1)
do.call(save, list(paste("PROPER.", i, sep = ""),
                   file = paste("(HIP)PROPER.fold.change.1.4.N.", sample.size[i], ".default.rdata", sep = "")))
#-----------------------------------------------
# plot
#-----------------------------------------------
sample.size = c(2, 4, 6, 8, 10, 20) #for predicted EDR
target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
# RNASeqDesign-----------
load("(HIP)RNASeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.2.predict.EDR.not.considering.lfc.rdata")
load("(HIP)RNASeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.4.predict.EDR.not.considering.lfc.rdata")
load("(HIP)RNASeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.10.predict.EDR.not.considering.lfc.rdata")
load("(HIP)RNASeqDesign.filter.by.5.fold.change.1.4.p.calibrated.2B.automate.CDD.pilot.12.predict.EDR.not.considering.lfc.rdata")
True.RNASeqDesign=apply(sapply(HIP.n.12[[1]], function(x) x[3,]), 1, median)

for(i in c(2, 4, 10)) assign(paste("p.", i, sep = ""),
                             Each.EDR(get(paste("HIP.n.", i, sep = "")), method = "RNASeqDesign",
                                      target.N = c(i:12,20,30,40,50,100),
                                      target.N.true = c(12, 20, 30, 40, 50, 100),
                                      True = True.RNASeqDesign,
                                      variation.True = F, sd.True = NULL))
# Poisson------------
load("(HIP)Poisson.fold.change.1.4.N.2.rdata")
load("(HIP)Poisson.fold.change.1.4.N.4.rdata")
load("(HIP)Poisson.fold.change.1.4.N.10.rdata")
load("(HIP)Poisson.fold.change.1.4.N.12.rdata")
for(i in c(1, 2, 5)) assign(paste("p.Poisson.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("Poisson.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     method = "Poisson", True = Poisson.6, target.N.true = c(12, 20, 30, 40, 50, 100),
                                     variation.True = F, sd.True = NULL))

# RNASeqPower--------
load("(HIP)RNASeqPower.fold.change.1.4.N.2.rdata")
load("(HIP)RNASeqPower.fold.change.1.4.N.4.rdata")
load("(HIP)RNASeqPower.fold.change.1.4.N.10.rdata")
load("(HIP)RNASeqPower.fold.change.1.4.N.12.rdata")

for(i in c(1, 2, 5)) assign(paste("p.RNASeqPower.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("RNASeqPower.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "RNASeqPower", True = RNASeqPower.6,
                                     variation.True = F, sd.True = NULL))

# NB------------
load("(HIP)NB.fold.change.1.4.N.2.rdata")
load("(HIP)NB.fold.change.1.4.N.4.rdata")
load("(HIP)NB.fold.change.1.4.N.10.rdata")
load("(HIP)NB.fold.change.1.4.N.12.rdata")
for(i in c(1, 2, 5)) assign(paste("p.NB.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("NB.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "NB", True = NB.6,
                                     variation.True = F, sd.True = NULL))

# PROPER------------
load("(HIP)PROPER.fold.change.1.4.N.2.default.rdata")
load("(HIP)PROPER.fold.change.1.4.N.4.default.rdata")
load("(HIP)PROPER.fold.change.1.4.N.10.default.rdata")
load("(HIP)PROPER.fold.change.1.4.N.12.default.rdata")
for(i in c(1, 2, 5)) assign(paste("p.PROPER.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("PROPER.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "PROPER", True = PROPER.6,
                                     variation.True = F, sd.True = NULL))
# Scotty----------
for(i in c(1, 2, 5)){
  tmp=read.table(paste("scotty_fold_change_1.4_result_n_", sample.size[i], ".txt", sep = ""), sep = ",")
  assign(paste("Scotty.", i, sep = ""), {
    x = sapply(1:10, function(j) tmp[(1+149*(j-1)):(149*(j-1)+149),][c(sample.size[i]:12,20,30,40,50,100)-1,10]/100)
  })
}
tmp=read.table(paste("scotty_fold_change_1.4_result_n_12.txt", sep = ""), sep = ",")
Scotty.6 = tmp[c(12,20,30,40,50,100)-1,10]/100

for(i in c(1, 2, 5)) assign(paste("p.Scotty.", sample.size[i], sep = ""),
                            Each.EDR(get(paste("Scotty.", i, sep = "")), target.N = c(sample.size[i]:12,20,30,40,50,100),
                                     target.N.true = c(12, 20, 30, 40, 50, 100),
                                     method = "Scotty", True = Scotty.6,
                                     variation.True = F, sd.True = NULL))
#---------------------------------------
# 2.2. TCGA ER data
#---------------------------------------
# subsampling
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))

repeat.num = 30
sample.size = c(2:20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130:150)
set.seed(12345)
Indicator.ER.pred = lapply(1:repeat.num,function(x) Sampling.function(sample.size, Data = TCGA.ER$Data, Status = TCGA.ER$status, 10, sample.ratio = 3))

for(i in 1:length(sample.size)) assign(paste("Indicator.ER.pred.", i, sep = ""),
                                       lapply(Indicator.ER.pred, function(x) list(x[[1]][[i]],x[[2]])))

a = paste("Indicator.ER.pred.", 1:length(sample.size), sep = "")
do.call(save, list(list = a, file = "TCGA.ER.indicator.N2-20.30.40.50.60.70.80.90.100.120.130-150.RData"))

# derive number of DE genes based on subsampling data
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))
load("TCGA.ER.indicator.N2-20.30.40.50.60.70.80.90.100.120.130-150.RData")
repeat.num = 30
sample.size = c(2:20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130:150)

library(snowfall)
sfInit(parallel = TRUE, cpus = 30)

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
  do.call(save, list(paste("ER.DE.", i, sep = ""), file = paste("TCGA.ER.N", sample.size[i], ".DE.rdata", sep = "")))
}


sample.size = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) #for predicted EDR
N.DE = c(2:20, 30, 40, 50, 60, 70, 80, 90, 100)

for(i in sample.size[1:6]) load(paste("TCGA.ER.N", i, ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.rdata", sep = ""))
for(i in N.DE) load(paste("TCGA.ER.N", i, ".DE.rdata", sep = ""))

FDR = 0.05
LFC = round(log2(1.40), 3); LFC

pdf("Figure S7.pdf")
result = DE.plot(N.DE, FDR = FDR, LFC = LFC, y.min = 0, y.max = 9000, q.full = NULL, lfc.full = NULL, DE.num.FDR.adjust = T)
abline(h = 6500, col = 3)
dev.off()

num.sig = result[[1]]
median.DE = result[[2]]; names(median.DE) = N.DE

drop.index.2 = which(!N.DE %in% c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

True.EDR = t(num.sig)[-drop.index.2,]/6500; True.EDR[True.EDR > 1] = 1; rownames(True.EDR) = N.DE[-drop.index.2]

True.CI = cbind(apply(True.EDR, 1, function(x) quantile(x, .25)), apply(True.EDR, 1, function(x) quantile(x, .75)))

#---------------------------
# predict EDR
#---------------------------
TCGA.ER = get(load("TCGA.ER.positive.negative.rdata"))

sample.size = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

set.seed(12345)
Indicator.ER.pred = lapply(1:10,function(x) Sampling.function(sample.size, Data = TCGA.ER$Data, Status = TCGA.ER$status, filter.cutoff = 10, sample.ratio = 3))

for(i in 1:length(sample.size)) assign(paste("Indicator.ER.pred.", i, sep = ""),
                                       lapply(Indicator.ER.pred, function(x) list(x[[1]][[i]],x[[2]])))

save(Indicator.ER.pred.1, Indicator.ER.pred.2, Indicator.ER.pred.3,
     Indicator.ER.pred.4, Indicator.ER.pred.5, Indicator.ER.pred.6,
     Indicator.ER.pred.7, Indicator.ER.pred.8, Indicator.ER.pred.9,
     Indicator.ER.pred.10, Indicator.ER.pred.11, Indicator.ER.pred.12,
     Indicator.ER.pred.13, Indicator.ER.pred.14,
     file = "TCGA.ER.indicator.UnBalance.3to1.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

load("TCGA.ER.indicator.UnBalance.3to1.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

sample.size = target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

#RNARNASeqDesign
library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfExportAll()

for(i in 1:6){
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
                                   target.N = target.N)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("ER.Predict.", i, ".2B", sep = ""), file = paste("TCGA.ER.N", target.N[i], ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.not.considering.lfc.rdata", sep = "")))
}

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfExportAll()

#Poisson
for(i in 1:6){
print(i)
t1=Sys.time()
sfExport("i")
assign(paste("ER.Poisson.", i, sep = ""),
       sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
         Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
         Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
         status = TCGA.ER$status[Ind[[1]]]

         Poisson(Data, status, target.N = target.N, rho = 1.4, DE.prop = .3)
       }))
t2=Sys.time()
print(t2-t1)
do.call(save, list(paste("ER.Poisson.", i, sep = ""), file = paste("TCGA.ER.Poisson.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:6){
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
do.call(save, list(paste("ER.RNASeqPower.", i, sep = ""), file = paste("TCGA.ER.RNASeqPower.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#NB
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("ER.NB.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
           Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.ER$status[Ind[[1]]]

           NB.exact(Data, status, target.N = target.N, rho = 1.4, DE.prop = 0.3)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("ER.NB.", i, sep = ""), file = paste("TCGA.ER.NB.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#PROPER
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("ER.PROPER.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.ER.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.ER.pred.", i, sep = ""))[[x]]
           Data = TCGA.ER$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.ER$status[Ind[[1]]]

           PROPER(Data, status, n.prop = 3, target.N = target.N, effect = log(1.4), DE.prop = 0.3)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("ER.PROPER.", i, sep = ""), file = paste("TCGA.ER.PROPER.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#Scotty, on Matlab
for(i in 1:6){
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

    write.table(Data, file = paste("ER.Data.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#---------------------------
# plot
#---------------------------
target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
target.N.true = c(20, 30, 40, 50, 60, 70, 80, 90, 100)
# RNASeqDesign-----
for(i in sample.size[1:6]) load(paste("TCGA.ER.N", i, ".predict.EDR.2B.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.Predict.", i, ".2B", sep = "")),
                                     method = "RNASeqDesign", target.N = target.N,
                                     target.N.true = target.N.true,
                                     pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Poisson-------
for(i in sample.size[1:6]) load(paste("TCGA.ER.Poisson.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.Poisson.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.Poisson.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "Poisson", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# RNASeqPower-------
for(i in sample.size[1:6]) load(paste("TCGA.ER.RNASeqPower.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.RNASeqPower.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.RNASeqPower.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "RNASeqPower", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# NB-------
for(i in sample.size[1:6]) load(paste("TCGA.ER.NB.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.NB.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.NB.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "NB", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))
# PROPER-------
for(i in sample.size[1:6]) load(paste("TCGA.ER.PROPER.fold.change.1.4.N", i, ".default.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.PROPER.", sample.size[i], ".ER", sep = ""),
                            Each.EDR(get(paste("ER.PROPER.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "PROPER", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Scotty------
for(i in 1:6){
  tmp=read.table(paste("er_result_fc_1.4_n_", sample.size[i], ".txt", sep = ""),sep = ",")
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

#---------------------------------------
# 2.3. TCGA Stage data
#---------------------------------------
# subsampling
TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))

repeat.num = 30
sample.size = c(2:20, 30, 40, 50, 60, 70, 80:120)
set.seed(12345)
Indicator.Stage.pred = lapply(1:repeat.num,function(x) Sampling.function(sample.size, Data = TCGA.Stage$Data, Status = TCGA.Stage$status, 10))

for(i in 1:length(sample.size)) assign(paste("Indicator.Stage.pred.", i, sep = ""),
                                       lapply(Indicator.Stage.pred, function(x) list(x[[1]][[i]],x[[2]])))

a = paste("Indicator.Stage.pred.", 1:length(sample.size), sep = "")
do.call(save, list(list = a, file = "TCGA.Stage.indicator.N2-20.30.40.50.60.70.80-120.RData"))

# derive number of DE genes based on subsampling data

TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))
load("TCGA.Stage.indicator.N2-20.30.40.50.60.70.80-120.RData")
repeat.num = 30
sample.size = c(2:20, 30, 40, 50, 60, 70, 80:120)

library(snowfall)
sfInit(parallel = TRUE, cpus = 30)

sfExportAll()

for(i in 1:length(sample.size)){
  print(i); sfExport("i")
  t1=Sys.time()
  assign(paste("Stage.DE.", i, sep = ""),
         sfLapply(1:repeat.num,function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]
           DE(Data, status, group.name = c("Early", "Late"), method = "GLM", FDR = 0.05)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.DE.", i, sep = ""), file = paste("TCGA.Stage.N", sample.size[i], ".DE.rdata", sep = "")))
}

sample.size = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
for(i in sample.size[1:6]) load(paste("TCGA.Stage.N", i, ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.rdata", sep = ""))

for(i in N.DE) load(paste("TCGA.Stage.N", i, ".DE.rdata", sep = ""))

FDR = 0.05
LFC = round(log2(1.40), 3)

pdf("Figure S8.pdf")
result = DE.plot(N.DE, FDR = FDR, LFC = LFC, dataset = "Stage", y.min = 300, y.max = 2100, q.full = NULL, lfc.full = NULL, DE.num.FDR.adjust = T)
abline(h = 1500, col = 3)
dev.off()

num.sig = result[[1]]
median.DE = result[[2]]

drop.index.2 = which(!N.DE %in% c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

True.EDR = t(num.sig)[-drop.index.2,]/1500; True.EDR[True.EDR > 1] = 1; rownames(True.EDR) = N.DE[-drop.index.2]

True.CI = cbind(apply(True.EDR, 1, function(x) quantile(x, .25)), apply(True.EDR, 1, function(x) quantile(x, .75)))

#---------------------------
# predict EDR
#---------------------------
TCGA.Stage = get(load("TCGA.Stage.early.late.rdata"))

set.seed(12345)
Indicator.Stage.pred = lapply(1:10,function(x) Sampling.function(c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), Data = TCGA.Stage$Data, Status = TCGA.Stage$status, 10))

for(i in 1:length(sample.size)) assign(paste("Indicator.Stage.pred.", i, sep = ""),
                                       lapply(Indicator.Stage.pred, function(x) list(x[[1]][[i]],x[[2]])))

save(Indicator.Stage.pred.1, Indicator.Stage.pred.2, Indicator.Stage.pred.3,
     Indicator.Stage.pred.4, Indicator.Stage.pred.5, Indicator.Stage.pred.6,
     Indicator.Stage.pred.7, Indicator.Stage.pred.8, Indicator.Stage.pred.9,
     Indicator.Stage.pred.10, Indicator.Stage.pred.11, Indicator.Stage.pred.12,
     Indicator.Stage.pred.13, Indicator.Stage.pred.14,
     file = "TCGA.Stage.indicator.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

load("TCGA.Stage.indicator.N2.4.6.8.10.20.30.40.50.60.70.80.90.100.RData")

sample.size = target.N = c(2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

library(snowfall)
sfInit(parallel = TRUE, cpus = 10)

sfExportAll()

for(i in 1:6){
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
                                   target.N = target.N)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.Predict.", i, ".2B", sep = ""), file = paste("TCGA.Stage.N", sample.size[i], ".filter.by.5.fold.change.1.4.p.calibrated.predict.EDR.2B.automate.CDD.not.considering.lfc.rdata", sep = "")))
}


library(snowfall)
sfInit(parallel = TRUE, cpus = 10)
sfExportAll()

#Poisson
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.Poisson.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]

           Poisson(Data, status, target.N = target.N, rho = 1.4, DE.prop = .1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.Poisson.", i, sep = ""), file = paste("TCGA.Stage.Poisson.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#RNASeqPower
for(i in 1:6){
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
  do.call(save, list(paste("Stage.RNASeqPower.", i, sep = ""), file = paste("TCGA.Stage.RNASeqPower.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#NB
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.NB.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]

           NB.exact(Data, status, target.N = target.N, rho = 1.4, DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.NB.", i, sep = ""), file = paste("TCGA.Stage.NB.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}
#PROPER
for(i in 1:6){
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("Stage.PROPER.", i, sep = ""),
         sfLapply(1:length(get(paste("Indicator.Stage.pred.", i, sep = ""))),function(x){
           Ind = get(paste("Indicator.Stage.pred.", i, sep = ""))[[x]]
           Data = TCGA.Stage$Data[-Ind[[2]], Ind[[1]]]
           status = TCGA.Stage$status[Ind[[1]]]

           PROPER(Data, status, n.prop = 1, target.N = target.N, effect = log(1.4), DE.prop = 0.1)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("Stage.PROPER.", i, sep = ""), file = paste("TCGA.Stage.PROPER.fold.change.1.4.N", target.N[i], ".rdata", sep = "")))
}

#Scotty
for(i in 1:6){
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

    write.table(Data, file = paste("Stage.Data.", sample.size[i], ".", j, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}

#---------------------------
# plot
#---------------------------
# RNASeqDesign-----
for(i in sample.size[1:6]) load(paste("TCGA.Stage.N", i, ".predict.EDR.2B.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.Predict.", i, ".2B", sep = "")),
                                     method = "RNASeqDesign",
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Poisson-------
for(i in sample.size[1:6]) load(paste("TCGA.Stage.Poisson.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.Poisson.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.Poisson.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "Poisson", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# RNASeqPower-------
for(i in sample.size[1:6]) load(paste("TCGA.Stage.RNASeqPower.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.RNASeqPower.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.RNASeqPower.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "RNASeqPower", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# NB-------
for(i in sample.size[1:6]) load(paste("TCGA.Stage.NB.fold.change.1.4.N", i, ".rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.NB.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.NB.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "NB", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))
# PROPER-------
for(i in sample.size[1:6]) load(paste("TCGA.Stage.PROPER.fold.change.1.4.N", i, ".default.rdata", sep = ""))
for(i in c(2, 5, 6)) assign(paste("p.PROPER.", sample.size[i], ".Stage", sep = ""),
                            Each.EDR(get(paste("Stage.PROPER.", i, sep = "")),
                                     target.N = target.N,
                                     target.N.true = target.N.true,
                                     method = "PROPER", pilot.n = sample.size[i], True = True.EDR,
                                     True.upper = True.CI[,1], True.lower = True.CI[,2],
                                     variation.True = T, sd.True = NULL))

# Scotty------
for(i in 1:6){
  tmp=read.table(paste("stage_result_fc_1.4_n_", sample.size[i], ".txt", sep = ""),sep = ",")
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
pdf("Figure 4.pdf", width = 18, height = 27)
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


#-----------------------------------------
# 3. Benchmark 3 - two way transformation
#-----------------------------------------
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.2.rdata")
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.4.rdata")
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.8.rdata")
load("Simulated.data.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.16.rdata")
#----------------------------------------------------
sample.size = c(2, 4, 8, 16)
target.N = c(5, 10, 20, 30, 40, 50, 100)
target.R = c(40, 80, 120, 160) * 10^6 #correspond to half lane, 1 lane, 1 and half lane, 2 lane

library(snowfall)
sfInit(parallel = TRUE, cpus = 20)

set.seed(12345)

n.repeat = length(simulation.N.2) #20

sfExportAll()

for(i in 1:length(sample.size)){ #16
  print(i)
  t1=Sys.time()
  sfExport("i")
  assign(paste("RNASeqDesign.", sample.size[i], sep = ""),
         sfLapply(1:n.repeat,function(x){
           Data = get(paste("simulation.N.", sample.size[i], sep = ""))[[x]]
           n = ncol(Data)/2
           status = c(rep("Case", n), rep("Control", n))
           Estimate.EDR.from.pilot(Data, status, group.name = c("Control", "Case"),
                                   FDR = 0.05, M = 20, filter = T, filter.level = 5,
                                   target.N = target.N, target.R = target.R)
         }))
  t2=Sys.time()
  print(t2-t1)
  do.call(save, list(paste("RNASeqDesign.", sample.size[i], sep = ""),
                     file = paste("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.",
                                  sample.size[i], ".predict.EDR.not.considering.lfc.2way.transformation.rdata", sep = "")))
}

#------------------------------
# summarize the result
#------------------------------
target.N = c(5, 10, 20, 30, 40, 50, 100)
target.R = paste(c(40, 80, 120, 160), "M", sep = "") #correspond to half lane, 1 lane, 1 and half lane, 2 lane

#load True
for(j in target.R){
  for(i in target.N){
    load(paste("(Simulation)True.EDR.", j, ".20.repeat.dispersion.5.lfc.mean.0.sd.0.2.N.", i,
               ".not.considering.lfc.rdata", sep = ""))
  }
}

#True EDR
True =
  sapply(target.R, function(j){
    sapply(target.N, function(i){
      x = get(paste("True.", i, ".", j, sep = ""))
      x = x[sapply(x, length) != 1]
      mean(sapply(x, function(y) y$Result[3]))
      # length(x)
    })
  }); True
rownames(True) = target.N
colnames(True) = target.R

#RNASeqDesign
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.2.predict.EDR.not.considering.lfc.2way.transformation.rdata")
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.4.predict.EDR.not.considering.lfc.2way.transformation.rdata")
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.8.predict.EDR.not.considering.lfc.2way.transformation.rdata")
load("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.5.lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.16.predict.EDR.not.considering.lfc.2way.transformation.rdata")

result.2 = sapply(1:20, function(i){ #repeatment
  ps = sapply(1:20, function(j){ #each posterior sampling
    sapply(RNASeqDesign.2[[i]]$Result[[j]], function(x) x[3,])
  })
  ps.median = matrix(apply(ps, 1, median), nrow = length(target.N))
})

result.4 = sapply(1:20, function(i){ #repeatment
  ps = sapply(1:20, function(j){ #each posterior sampling
    sapply(RNASeqDesign.4[[i]]$Result[[j]], function(x) x[3,])
  })
  ps.median = matrix(apply(ps, 1, median), nrow = length(target.N))
})

result.8 = sapply(1:20, function(i){ #repeatment
  ps = sapply(1:20, function(j){ #each posterior sampling
    sapply(RNASeqDesign.8[[i]]$Result[[j]], function(x) x[3,])
  })
  ps.median = matrix(apply(ps, 1, median), nrow = length(target.N))
})

result.16 = sapply(1:20, function(i){ #repeatment
  ps = sapply(1:20, function(j){ #each posterior sampling
    sapply(RNASeqDesign.16[[i]]$Result[[j]], function(x) x[3,])
  })
  ps.median = matrix(apply(ps, 1, median), nrow = length(target.N))
})

z = cbind(apply(result.2, 2, function(x) sqrt(mean((x - as.vector(True))^2))),
          apply(result.4, 2, function(x) sqrt(mean((x - as.vector(True))^2))),
          apply(result.8, 2, function(x) sqrt(mean((x - as.vector(True))^2))),
          apply(result.16, 2, function(x) sqrt(mean((x - as.vector(True))^2))))

colnames(z) = c("2", "4", "8", "16")
z = data.frame(z, check.names = F)

library(tidyr)
z1 = gather(z, Pilot, RMSE, factor_key = T)
library(ggplot2)

pdf("Figure S6.pdf")
ggplot(data = z1, aes(x = Pilot, y = RMSE)) + geom_boxplot() + xlab(bquote("Pilot n"[0]))
dev.off()

#---------------------------------------
# 4. Other dispersion (Figure S5)
#---------------------------------------
library(ggplot2)
library(drc)
library(DEoptim)

sample.size=c(5, 10, 20, 30, 40, 50, 100)
pilot.size = c(2, 4, 8, 16)
dispersion = c(2, 5, 10, 20)

#----------------------
#True EDR
#----------------------
True.EDR = lapply(dispersion, function(d){
  sapply(sample.size, function(n){
    print(n)
    x = get(load(paste("(Simulation)True.EDR.20M.20.repeat.dispersion.", d, ".lfc.mean.0.sd.0.2.N.", n, ".not.considering.lfc.rdata", sep = "")))
    x = x[sapply(x, length) != 1]
    mean(sapply(x, function(y) y$Result[3]))
  })
})
names(True.EDR) = dispersion

par(mfrow = c(1, length(dispersion)))
for(i in 1:length(dispersion)) plot(sample.size, True.EDR[[i]], type = "b", ylim = c(0, 1), main = dispersion[i])

#----------------
Power.max = 0.8
N.hat.true = 17

Predict.EDR = lapply(dispersion, function(d){
  print(d)
  lapply(pilot.size, function(n){
    x = get(load(paste("(Simulation)RNASeqDesign.20M.20.repeat.dispersion.", d, ".lfc.mean.0.sd.0.2.filter.by.5.2B.automate.CDD.pilot.", n, ".predict.EDR.not.considering.lfc.rdata", sep = "")))
    Each.EDR(x, method = "RNASeqDesign",
             True = True.EDR[[as.character(d)]], target.N = sample.size,
             variation.True = F, sd.True = NULL, output.MSE = T, Power.max = Power.max)[[1]]
  })
})

pdf("Figure S5.pdf", width = 12, height = 9)
multiplot(Predict.EDR[[1]][[1]], Predict.EDR[[1]][[2]], Predict.EDR[[1]][[3]], Predict.EDR[[1]][[4]],
          Predict.EDR[[2]][[1]], Predict.EDR[[2]][[2]], Predict.EDR[[2]][[3]], Predict.EDR[[2]][[4]],
          Predict.EDR[[3]][[1]], Predict.EDR[[3]][[2]], Predict.EDR[[3]][[3]], Predict.EDR[[3]][[4]],
          Predict.EDR[[4]][[1]], Predict.EDR[[4]][[2]], Predict.EDR[[4]][[3]], Predict.EDR[[4]][[4]],
          cols=5, layout = matrix(1:16, ncol = 4, byrow = F))
dev.off()


#-----------------------------------------------------
# 5. ROC curves of Wald test vs exact test (Figure S3)
#-----------------------------------------------------
library(snowfall)
sfInit(parallel=TRUE, cpus=20)
load("Normalized.HIP.data.rdata")

par.setting = Parameter.Estimate(Normalized.HIP.data, 650)

par.setting$lfc.mean.sd = c(0, 0.2)
par.setting$Prop.DE = 0.1
log.fold.change = rbind(c(-0.20,0.20),c(-0.26,0.26),c(-0.32,0.32),c(-0.38,0.38))

sfExportAll()
for(dispersion in c(40, 50, 60)){
  print(dispersion)

  par.setting$dispersion = dispersion

  sfExport("par.setting")

  assign(paste("HIP.Data.Dispersion.", dispersion, ".N.4", sep = ""), lapply(1:4,function(i){
    par.setting$lfc = log.fold.change[i,]
    Generate.data.Setting.adjust(4, 50, 10^4, par.setting)
  }))
  do.call(save, list(paste("HIP.Data.Dispersion.", dispersion, ".N.4", sep = ""),
                     file = paste("HIP.Data.Dispersion.", dispersion, ".N.4.new.rdata", sep = "")))
}

library(snowfall)
sfInit(parallel=TRUE, cpus=30)
sfExportAll()
sfLibrary(ROCR)
sfLibrary(edgeR)

pdf("Figure S3.pdf", height = 10.5, width = 14)
par(mfrow = c(3, 4))

auc.result = lapply(c(40, 50, 60), function(dispersion){
  lapply(1:4, function(i){ #for different effect size
    sfExport("dispersion", "i")
    print(paste("dispersion", dispersion, "effect size", i))
    p.value = sfLapply(1:50,function(j) Compare.TEST(get(paste("HIP.Data.Dispersion.", dispersion, ".N.4", sep = ""))[[i]][[j]]))

    exact = sapply(p.value, function(x) -x[[1]][,1])
    wald = sapply(p.value, function(x) -x[[1]][,3])

    labels = sapply(p.value, function(x) x[[2]])

    pred.1 = prediction(exact, labels)
    perf.1 <- performance(pred.1, measure = "tpr", x.measure = "fpr")

    pred.2 = prediction(wald, labels)
    perf.2 <- performance(pred.2, measure = "tpr", x.measure = "fpr")

    auc.1 <- performance(pred.1,"auc")
    auc.1 <- unlist(slot(auc.1, "y.values"))
    auc.2 <- performance(pred.2,"auc")
    auc.2 <- unlist(slot(auc.2, "y.values"))
    ## produce ROC plot
    plot(perf.1, col = "blue", avg = "vertical", spread.estimate = "boxplot", main = paste("Dispersion:", dispersion, "; LFC:", log.fold.change[i,2], sep = ""))
    plot(perf.2, col = "red", avg = "vertical", spread.estimate = "boxplot", add = T)

    return(list(auc.exact = auc.1, auc.wald = auc.2))
  })
})
dev.off()

auc.wald = lapply(auc.result, function(x) t(sapply(x, function(y) round(c(mean(y$auc.wald), sd(y$auc.wald)), 2))))
auc.exact = lapply(auc.result, function(x) t(sapply(x, function(y) round(c(mean(y$auc.exact), sd(y$auc.exact)), 2))))
#------------------------------------------------
# 6. Wald test vs exact test and LRT (Figure S4)
#------------------------------------------------
load("Data.SRA.PFC.rdata")
load("Data.SRA.HIP.rdata")
load("Data.SRA.STR.rdata")

pdf("Figure S4.pdf", height = 21, width = 14)
par(mfrow = c(3, 2))
Data.SRA.HIP=Data.SRA.HIP[-c((nrow(Data.SRA.HIP)-4):nrow(Data.SRA.HIP)),]
tmp.1=my.QQ(Data.SRA.HIP,"HIP")

Data.SRA.PFC=Data.SRA.PFC[-c((nrow(Data.SRA.PFC)-4):nrow(Data.SRA.PFC)),]
tmp.2=my.QQ(Data.SRA.PFC,"PFC")

Data.SRA.STR=Data.SRA.STR[-c((nrow(Data.SRA.STR)-4):nrow(Data.SRA.STR)),]
tmp.3=my.QQ(Data.SRA.STR,"STR")
dev.off()
