## Combine TSCB strategy of constraint-based methods with score-based algorithms

#### Accuracy Experiments

data(alarm)
Alarm = as.data.frame(lapply(alarm, as.factor))

modelstring = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
                    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]",
                    "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]",
                    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]",
                    "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]",
                    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]",
                    "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
Alarm_actual = model2network(modelstring)



library(ggplot2)
alarm_acc2 = numeric(ncol(alarm))
for(i in 1:ncol(alarm)){
  bn.combine = cluster_BNL(Alarm, algo = gs, N = i)
  undirect_arcs = undirected.arcs(bn.combine)
  undirect_arcs1 = apply(undirect_arcs, MARGIN = 1, function(item){
    return(paste(item[1], item[2]))
  })
  arcs1 = apply(bn.combine$arcs, MARGIN = 1, function(item){
    return(paste(item[1], item[2]))
  })
  bn.combine$arcs = bn.combine$arcs[!(arcs1 %in% undirect_arcs1), ]
  bn.combine = hc(alarm, start = bn.combine)
  
  alarm_acc2[i] = ACC(actual = Alarm_actual$arcs, learning = bn.combine$arcs, n = ncol(alarm))
}
alarm_acc22 = data.frame(NG = 1:ncol(alarm), Accuracy = c(alarm_acc2, Alarm_acc1$Accuracy, Alarm_acc2$Accuracy), class = as.factor(rep(c("Combine", "TSCB", "HC"), each = ncol(Alarm))))

ggplot(alarm_acc22, mapping = aes(NG, Accuracy, color = class, shape = class, linetype = class)) + geom_hline(yintercept = Alarm_acc1$Accuracy[nrow(Alarm_acc1)], color = "black", size = 2, linetype = 2) + geom_point(size = 4) + geom_line(size = 1.3) + labs(x = "The Number of Clusters", y = "Accuracy") + scale_shape_discrete(name = "", labels = c("Combined TSCB and HC", "TSCB With GS", "HC With TSCB")) + scale_color_discrete(breaks = NULL) + scale_linetype_discrete(breaks = NULL) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.position = c(0.2, 0.9), legend.text = element_text(size = 20))


accuracy_combine = function(data_actual, algo1 = "gs", algo2 = "hc", n_run = 50, n_samples = 1000, ...){
  require(bnlearn)
  Tra1_acc = numeric(length = n_run)
  Tra2_acc = numeric(length = n_run)
  TSCB_acc = numeric(length = n_run)
  Combine_acc = numeric(length = n_run)
  actual_arcs = arcs(data_actual)
  
  for(i in 1:n_run){
    data = rbn(data_actual, n_samples)
    
    accuracy1 = numeric(ncol(data))
    accuracy2 = numeric(ncol(data))
    flag = 0
    
    for(j in 1:ncol(data)){
      bn_1 = cluster_BNL(data, N = j, algo = get(algo1))
      undi_arcs = undirected.arcs(bn_1)
      if(length(undi_arcs)){
        undi_arcs1 = apply(undi_arcs, MARGIN = 1, function(item){
          return(paste(item[1], item[2]))
        })
        arcs1 = apply(bn_1$arcs, MARGIN = 1, function(item){
          return(paste(item[1], item[2]))
        })
        bn_1$arcs = bn_1$arcs[!(arcs1 %in% undi_arcs1), ]
      }
      re = NULL
      try((re = bn.combine = get(algo2)(data, start = bn_1)), silent = TRUE)
      if(is.null(re)){
        flag = 1
        break
      }
      accuracy1[j] = ACC(actual_arcs, learning = cluster_BNL(data, N = j, algo = get(algo1))$arcs, n = ncol(data))
      accuracy2[j] = ACC(actual_arcs, learning = bn.combine$arcs, n = ncol(data))
    }
    if(flag){
      i = i-1
      break
    }
    TSCB_acc[i] = max(accuracy1)
    Tra1_acc[i] = accuracy1[ncol(data)]
    Tra2_acc[i] = ACC(actual_arcs, learning = get(algo2)(data)$arcs, n = ncol(data))
    Combine_acc[i] = max(accuracy2)
  }
  record = data.frame(Tra1_acc, Tra2_acc, TSCB_acc, Combine_acc)
  names(record) = c(algo1, algo2, "TSCB", "Combine")
  return(record)
}


# Embed GS, HC
# "alarm"
alarm_actual = readRDS("alarm.rds")
alarm_acc_gs_hc = accuracy_combine(alarm_actual, algo1 = "gs", algo2 = "hc")

# "asia"
asia_actual = readRDS("asia.rds")
asia_acc_gs_hc = accuracy_combine(asia_actual, algo1 = "gs", algo2 = "hc")

# "insurance"
insurance_actual = readRDS("insurance.rds")
insurance_acc_gs_hc = accuracy_combine(insurance_actual, algo1 = "gs", algo2 = "hc")

# "hepar2"
hepar2_actual = readRDS("hepar2.rds")
hepar2_acc_gs_hc = accuracy_combine(hepar2_actual, algo1 = "gs", algo2 = "hc")


# Embed IAMB, HC
# "alarm"
alarm_acc_iamb_hc = accuracy_combine(alarm_actual, algo1 = "iamb", algo2 = "hc")

# "asia"
asia_acc_iamb_hc = accuracy_combine(asia_actual, algo1 = "iamb", algo2 = "hc")

# "insurance"
insurance_acc_iamb_hc = accuracy_combine(insurance_actual, algo1 = "iamb", algo2 = "hc")

# "hepar2"
hepar2_acc_iamb_hc = accuracy_combine(hepar2_actual, algo1 = "iamb", algo2 = "hc")


# Embed inter-IAMB, HC
# "alarm"
alarm_acc_inter_IAMB_hc = accuracy_combine(alarm_actual, algo1 = "inter.iamb", algo2 = "hc")

# "asia"
asia_acc_inter_IAMB_hc = accuracy_combine(asia_actual, algo1 = "inter.iamb", algo2 = "hc")

# "insurance"
insurance_acc_inter_IAMB_hc = accuracy_combine(insurance_actual, algo1 = "inter.iamb", algo2 = "hc")

# "hepar2"
hepar2_acc_inter_IAMB_hc = accuracy_combine(hepar2_actual, algo1 = "inter.iamb", algo2 = "hc")


# Embed MMPC, HC
# "alarm"
alarm_acc_mmpc_hc = accuracy_combine(alarm_actual, algo1 = "mmpc", algo2 = "hc")

# "asia"
asia_acc_mmpc_hc = accuracy_combine(asia_actual, algo1 = "mmpc", algo2 = "hc")

# "insurance"
insurance_acc_mmpc_hc = accuracy_combine(insurance_actual, algo1 = "mmpc", algo2 = "hc")

# "hepar2"
hepar2_acc_mmpc_hc = accuracy_combine(hepar2_actual, algo1 = "mmpc", algo2 = "hc")


# Embed GS, TABU
# "alarm"
alarm_acc_gs_tabu = accuracy_combine(alarm_actual, algo1 = "gs", algo2 = "tabu")

# "asia"
asia_acc_gs_tabu = accuracy_combine(asia_actual, algo1 = "gs", algo2 = "tabu")

# "insurance"
insurance_acc_gs_tabu = accuracy_combine(insurance_actual, algo1 = "gs", algo2 = "tabu")

# "hepar2"
hepar2_acc_gs_tabu = accuracy_combine(hepar2_actual, algo1 = "gs", algo2 = "tabu")


# Embed IAMB, TABU
# "alarm"
alarm_acc_iamb_tabu = accuracy_combine(alarm_actual, algo1 = "iamb", algo2 = "tabu")

# "asia"
asia_acc_iamb_tabu = accuracy_combine(asia_actual, algo1 = "iamb", algo2 = "tabu")

# "insurance"
insurance_acc_iamb_tabu = accuracy_combine(insurance_actual, algo1 = "iamb", algo2 = "tabu")

# "hepar2"
hepar2_acc_iamb_tabu = accuracy_combine(hepar2_actual, algo1 = "iamb", algo2 = "tabu")


# Embed inter-IAMB, TABU
# "alarm"
alarm_acc_inter_IAMB_tabu = accuracy_combine(alarm_actual, algo1 = "inter.iamb", algo2 = "tabu")

# "asia"
asia_acc_inter_IAMB_tabu = accuracy_combine(asia_actual, algo1 = "inter.iamb", algo2 = "tabu")

# "insurance"
insurance_acc_inter_IAMB_tabu = accuracy_combine(insurance_actual, algo1 = "inter.iamb", algo2 = "tabu")

# "hepar2"
hepar2_acc_inter_IAMB_tabu = accuracy_combine(hepar2_actual, algo1 = "inter.iamb", algo2 = "tabu")


# Embed MMPC, TABU
# "alarm"
alarm_acc_mmpc_tabu = accuracy_combine(alarm_actual, algo1 = "mmpc", algo2 = "tabu")

# "asia"
asia_acc_mmpc_tabu = accuracy_combine(asia_actual, algo1 = "mmpc", algo2 = "tabu")

# "insurance"
insurance_acc_mmpc_tabu = accuracy_combine(insurance_actual, algo1 = "mmpc", algo2 = "tabu")

# "hepar2"
hepar2_acc_mmpc_tabu = accuracy_combine(hepar2_actual, algo1 = "mmpc", algo2 = "tabu")


#### Time Experiments
time_comparison_com = function(data, kmax, algo1 = "gs", algo2 = "hc", n_run = 100, meth_index = NULL){
  
  if(is.null(meth_index)) meth_index = "Combined"
  
  if(meth_index == "Tra"){
    time_elapsed = numeric(length = n_run)
    time_elapsed = sapply(time_elapsed, function(x){
      return(system.time(get(algo1)(data))[3])
    })
  }else{
    time_elapsed = numeric(length = n_run)
    time_elapsed = sapply(time_elapsed, function(x){
      return(system.time({
        bn_1 = cluster_BNL(data, N = kmax, algo = get(algo1))
        undi_arcs = undirected.arcs(bn_1)
        undi_arcs1 = apply(undi_arcs, MARGIN = 1, function(item){
          return(paste(item[1], item[2]))
        })
        arcs1 = apply(bn_1$arcs, MARGIN = 1, function(item){
          return(paste(item[1], item[2]))
        })
        # Discard the undirected edges
        bn_1$arcs = bn_1$arcs[!(arcs1 %in% undi_arcs1), ]
        bn.combine = get(algo2)(data, start = bn_1)
      })[3])
    })
  }
  
  time_com = data.frame(time = time_elapsed, group = rep(meth_index, n_run), Kindex = rep(kmax, n_run))
  
  return(time_com)
}

data(alarm)
Alarm = as.data.frame(lapply(alarm, as.factor))
t_com = data.frame()
for(i in 1:ncol(Alarm)){
  if(i == ncol(Alarm)){
    t_com = rbind(t_com, time_comparison_com(Alarm, kmax = i, algo1 = "gs", n_run = 100, meth_index = "Tra"))
  }else{
    t_com = rbind(t_com, time_comparison_com(Alarm, kmax = i, algo1 = "gs", algo2 = "hc", n_run = 100))
  }
}

t_com1 = t_com
t_com1$group = as.character(t_com1$group)
t_com1[t_com1$Kindex == "37", ]$group = "GS"
t_com1$group = as.factor(t_com1$group)
t_com1$Kindex = as.factor(t_com1$Kindex)
ggplot(t_com1, mapping = aes(x = Kindex, y = time, color = group)) + geom_boxplot() + scale_color_manual(name = "", labels = c("Combine", "Pure GS"), values = c("cyan3", "black")) + labs(x = "The Number of Clusters", y = "Time (Sec)")  + theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), legend.position = c(0.8,0.85), axis.title = element_text(size = 45), axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 15), legend.text = element_text(size = 30))
