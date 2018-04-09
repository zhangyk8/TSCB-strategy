# Time efficiency experiments on synthetic datasets

### Calculate the running time in each step

time_step = function(data, stand = TRUE, method = "average", N = floor(sqrt(2)*ncol(data)/2), algo = hc, n_run = 100, ...){
  require(bnlearn)
  require(cluster)
  
  time_elapsed = matrix(0, ncol = 3, nrow = n_run)
  time_elapsed = apply(time_elapsed, MARGIN = 1, function(Time){
    
    Time[1] = system.time({
      data1 = as.data.frame(lapply(data, as.numeric))
      data1 = apply(data1, MARGIN = 2, function(column){
        column - mean(unique(column))
      })
      dissim = as.dist(matrix(1, nrow = ncol(data1), ncol = ncol(data1)) - cor(data1))
      dissim[is.na(dissim)] = 2
      cluster1 = agnes(dissim, diss = inherits(dissim, "dist"), stand = FALSE, method = method)
      Group = cutree(cluster1, k = N)
    })[3]
    
    Time[2] = system.time({
      whiteArcs = data.frame(from = c(), to = c())
      for(i in 1:N){
        index = (Group == i)
        
        if(sum(index) > 1){
          bn.cluster = algo(data[ ,index])
          whiteArcs = rbind(whiteArcs, bn.cluster$arcs)
        } 
      }
    })[3]
    
    Time[3] = system.time({
      if(nrow(whiteArcs) == 0) whiteArcs = NULL
      bn.cl1 = algo(data, whitelist = whiteArcs)
    })[3]
    return(Time)
  })
  time_elapsed = t(time_elapsed)
  colnames(time_elapsed) = c("cluster", "learn_in_cluster", "combine")
  return(time_elapsed)
}

# "alarm" data set
alarm_actual = readRDS("alarm.rds")
t_alarm1 = data.frame()
for(i in 1:50){
  data = rbn(alarm_actual, 2000)
  
  acc1 = numeric(ncol(data))
  for(i in 1:ncol(data)){
    acc1[i] = ACC(actual = arcs(alarm_actual), learning = cluster_BNL(data, N = i, algo = gs)$arcs, n = ncol(data))
  }
  
  t50 = time_step(data, algo = gs, n_run = 10, N = which.max(acc1))
  t50 = as.data.frame(t50)
  t50$tra = sapply(rep(0, 10), function(t){
    return(system.time(gs(data))[3])
  })
  
  t_alarm1 = rbind(t_alarm1, sapply(t50, mean))
}
colnames(t_alarm1) = c("cluster", "learn_in_cluster", "combine", "tra")

# "asia" data set
asia_actual = readRDS("asia.rds")
t_asia1 = data.frame()
for(i in 1:50){
  data = rbn(asia_actual, 2000)
  
  acc1 = numeric(ncol(data))
  for(i in 1:ncol(data)){
    acc1[i] = ACC(actual = arcs(asia_actual), learning = cluster_BNL(data, N = i, algo = gs)$arcs, n = ncol(data))
  }
  
  t50 = time_step(data, algo = gs, n_run = 10, N = which.max(acc1))
  t50 = as.data.frame(t50)
  t50$tra = sapply(rep(0, 10), function(t){
    return(system.time(gs(data))[3])
  })
  
  t_asia1 = rbind(t_asia1, sapply(t50, mean))
}
colnames(t_asia1) = c("cluster", "learn_in_cluster", "combine", "tra")

# "insurance" data set
insurance_actual = readRDS("insurance.rds")
t_insurance1 = data.frame()
for(i in 1:50){
  data = rbn(insurance_actual, 2000)
  
  acc1 = numeric(ncol(data))
  for(j in 1:ncol(data)){
    try((acc1[j] = ACC(actual = arcs(insurance_actual), learning = cluster_BNL(data, N = j, algo = gs)$arcs, n = ncol(data))), silent = TRUE)
    if(acc1[j]){
      j = j-1
      break
    }
  }
  
  t50 = time_step(data, algo = gs, n_run = 10, N = which.max(acc1))
  t50 = as.data.frame(t50)
  t50$tra = sapply(rep(0, 10), function(t){
    return(system.time(gs(data))[3])
  })
  
  t_insurance1 = rbind(t_insurance1, sapply(t50, mean))
}
colnames(t_insurance1) = c("cluster", "learn_in_cluster", "combine", "tra")

# "hepar2" data set
hepar2_actual = readRDS("hepar2.rds")
t_hepar21 = data.frame()
for(i in 1:50){
  data = rbn(hepar2_actual, 5000)
  
  acc1 = numeric(ncol(data)-9)
  for(i in 10:ncol(data)){
    acc1[i-9] = ACC(actual = arcs(hepar2_actual), learning = cluster_BNL(data, N = i, algo = gs)$arcs, n = ncol(data))
  }
  
  t50 = time_step(data, algo = gs, n_run = 10, N = which.max(acc1)+9)
  t50 = as.data.frame(t50)
  t50$tra = sapply(rep(0, 10), function(t){
    return(system.time(gs(data))[3])
  })
  
  t_hepar21 = rbind(t_hepar21, sapply(t50, mean))
}
colnames(t_hepar21) = c("cluster", "learn_in_cluster", "combine", "tra")



### Time Variation on the "alarm" dataset

time_comparison1 = function(data, kmax, algo = gs, meth = cluster_BNL, n_run = 100, meth_index = NULL){
  
  if(is.null(meth_index)) meth_index = "Retain all arcs"
  
  if(meth_index == "Tra"){
    time_elapsed = numeric(length = n_run)
    time_elapsed = sapply(time_elapsed, function(x){
      return(system.time(algo(data))[3])
    })
  }else{
    time_elapsed = numeric(length = n_run)
    time_elapsed = sapply(time_elapsed, function(x){
      return(system.time(meth(data, N = kmax, algo = algo))[3])
    })
  }
  
  time_com = data.frame(time = time_elapsed, group = rep(meth_index, n_run), Kindex = rep(kmax, n_run))
  
  return(time_com)
}

data(alarm)
Alarm = as.data.frame(lapply(alarm, as.factor))
t_gs = data.frame()
for(i in 1:ncol(Alarm)){
  if(i == ncol(Alarm)){
    t_gs = rbind(t_gs, time_comparison1(Alarm, kmax = i, algo = gs, n_run = 200, meth_index = "Tra"))
  }else{
    t_gs = rbind(t_gs, time_comparison1(Alarm, kmax = i, algo = gs, n_run = 200, meth = cluster_BNL))
  }
}

t_gs1 = t_gs
t_gs1$group = as.character(t_gs1$group)
t_gs1[t_gs1$Kindex == "37", ]$group = "Traditional Method"
t_gs1$group = as.factor(t_gs1$group)
t_gs1$Kindex = as.factor(t_gs1$Kindex)
ggplot(t_gs1, mapping = aes(x = Kindex, y = time, color = group)) + geom_boxplot() + scale_color_manual(name = "", labels = c("TSCB", "Traditional Method"), values = c("cyan3", "black")) + labs(x = "The Number of Clusters", y = "Time (Sec)")  + theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), legend.position = c(0.8,0.85), axis.title = element_text(size = 45), axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 15), legend.text = element_text(size = 30))


### Elapsed Time Comparison of different traditional methods

modelstring = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                    "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                    "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                    "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
Alarm_actual = model2network(modelstring)

acc_al_op = data.frame(TSCB = rep(0, 6), traditional = rep(0, 6))
acc_al_op$TSCB = sapply(c(gs, iamb, inter.iamb, mmpc, hc, tabu), function(algorithm){
  
  acc1 = numeric(length = ncol(Alarm))
  for(i in 1:ncol(Alarm)){
    acc1[i] = ACC(actual = Alarm_actual$arcs, learning = cluster_BNL(alarm, algo = algorithm, N = i)$arcs, n = ncol(alarm))
  }
  
  n_max = which.max(acc1)
  
  return(ACC(actual = Alarm_actual$arcs, learning = cluster_BNL(alarm, algo = algorithm, N = n_max)$arcs, n = ncol(alarm)))
})

acc_al_op$optimal_K = sapply(c(gs, iamb, inter.iamb, mmpc, hc, tabu), function(algorithm){
  acc1 = numeric(length = ncol(Alarm))
  
  for(i in 1:ncol(Alarm)){
    acc1[i] = ACC(actual = Alarm_actual$arcs, learning = cluster_BNL(alarm, algo = algorithm, N = i)$arcs, n = ncol(alarm))
  }
  
  return(which.max(acc1))
})

acc_al_op$traditional = sapply(c(gs, iamb, inter.iamb, mmpc, hc, tabu), function(algorithm){
  return(ACC(Alarm_actual$arcs, learning = algorithm(alarm)$arcs, n = ncol(alarm)))
})

acc_al_op$Method = c("gs", "iamb", "inter.iamb", "mmpc", "hc", "tabu")



time_comparison1 = function(data, kmax, algo = gs, n_run = 100, meth_index = NULL){
  
  meth = cluster_BNL
  if(is.null(meth_index)) meth_index = "Retain all arcs"
  
  time_elapsed = numeric(length = n_run)
  time_elapsed = sapply(time_elapsed, function(x){
    return(system.time(meth(data, N = kmax, algo = algo))[3])
  })
  
  time_com = data.frame(time = time_elapsed, group = rep(meth_index, n_run), Kindex = rep(kmax, n_run))
  
  return(time_com)
}

result = data.frame()
algorithm = c(gs, iamb, inter.iamb, mmpc, hc, tabu)

for(i in 1:6){
  Time = time_comparison1(alarm, kmax = acc_al_op$optimal_K[i], algo = algorithm[[i]], n_run = 200)$time
  time_tra = sapply(rep(0, 200), function(t){
    return(system.time(algorithm[[i]](Alarm))[3])
  })
  Time = c(Time, time_tra)
  
  result = rbind(result, data.frame(Time, group = c(rep("Retain all arcs", 200), rep("Traditional Method", 200)), Kindex = rep(acc_al_op$optimal_K[i], 200)))
}

result$Method = rep(c("gs", "iamb", "inter.iamb", "mmpc", "hc", "tabu"), each = 400)

ggplot(result) + geom_boxplot(mapping = aes(x = Method, y = Time, color = group)) + scale_color_manual(name = "", labels = c("TSCB", "Traditional Methods"), values = c("cyan3", "black")) + labs(x = "Traditional Methods", y = "Time (Sec)")  + theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), legend.position = c(0.8,0.9), axis.title = element_text(size = 45), axis.text = element_text(size = 28), legend.text = element_text(size = 29))

