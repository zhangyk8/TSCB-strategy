# Experiments of Accuracies on synthetic benchmark datasets

### Embed Different Methods

accuracy_anal = function(data_actual, algo = gs, n_run = 100, n_samples = 1000, ...){
  require(bnlearn)
  TSCB_acc = numeric(length = n_run)
  Tra_acc = numeric(length = n_run)
  actual_arcs = arcs(data_actual)
  
  for(i in 1:n_run){
    data = rbn(data_actual, n_samples)
    
    accuracy = numeric(ncol(data))
    for(j in 1:ncol(data)){
      re = NULL
      try((re = ACC(actual_arcs, learning = cluster_BNL(data, N = j, algo = algo)$arcs, n = ncol(data))), silent = T)
      if(is.null(re)){
        j = j-1
        break
      }else{
        accuracy[j] = re
      }
    }
    
    TSCB_acc[i] = max(accuracy)
    Tra_acc[i] = accuracy[ncol(data)]
  }
  record = data.frame(TSCB = TSCB_acc, Traditional = Tra_acc)
  return(record)
}

# Embed GS
# "alarm"
alarm_actual = readRDS("alarm.rds")
alarm_acc_gs = accuracy_anal(alarm_actual, algo = gs)

# "asia"
asia_actual = readRDS("asia.rds")
asia_acc_gs = accuracy_anal(asia_actual, algo = gs)

# "insurance" (Run Step-by-step)
insurance_actual = readRDS("insurance.rds")
insurance_acc_gs = accuracy_anal(insurance_actual, algo = gs)

# "hepar2"
hepar2_actual = readRDS("hepar2.rds")
hepar2_acc_gs = accuracy_anal(hepar2_actual, algo = gs)


# Embed IAMB
# "alarm"
alarm_acc_iamb = accuracy_anal(alarm_actual, algo = iamb)

# "asia"
asia_acc_iamb = accuracy_anal(asia_actual, algo = iamb)

# "insurance"
insurance_acc_iamb = accuracy_anal(insurance_actual, algo = iamb)

# "hepar2"
hepar2_acc_iamb = accuracy_anal(hepar2_actual, algo = iamb)


# Embed Iter-IAMB
# "alarm"
alarm_acc_inter_iamb = accuracy_anal(alarm_actual, algo = inter.iamb)

# "asia"
asia_acc_inter_iamb = accuracy_anal(asia_actual, algo = inter.iamb)

# "insurance"
insurance_acc_inter_iamb = accuracy_anal(insurance_actual, algo = inter.iamb)

# "hepar2"
hepar2_acc_inter_iamb = accuracy_anal(hepar2_actual, algo = inter.iamb)


# Embed MMPC
# "alarm"
alarm_acc_mmpc = accuracy_anal(alarm_actual, algo = mmpc)

# "asia"
asia_acc_mmpc = accuracy_anal(asia_actual, algo = mmpc)

# "insurance"
insurance_acc_mmpc = accuracy_anal(insurance_actual, algo = mmpc)

# "hepar2"
hepar2_acc_mmpc = accuracy_anal(hepar2_actual, algo = mmpc)


# Embed HC
# "alarm"
alarm_acc_hc = accuracy_anal(alarm_actual, algo = hc)

# "asia"
asia_acc_hc = accuracy_anal(asia_actual, algo = hc)

# "insurance"
insurance_acc_hc = accuracy_anal(insurance_actual, algo = hc)

# "hepar2"
hepar2_acc_hc = accuracy_anal(hepar2_actual, algo = hc)


# Embed TABU
# "alarm"
alarm_acc_tabu = accuracy_anal(alarm_actual, algo = tabu)

# "asia"
asia_acc_tabu = accuracy_anal(asia_actual, algo = tabu)

# "insurance"
insurance_acc_tabu = accuracy_anal(insurance_actual, algo = tabu)

# "hepar2"
hepar2_acc_tabu = accuracy_anal(hepar2_actual, algo = tabu)



### Accuracy Analysis on "alarm" dataset

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

# GS Algorithm
library(ggplot2)
Alarm_acc1 = numeric(ncol(Alarm))
for(i in 1:ncol(Alarm)){
  Alarm_acc1[i] = ACC(actual = Alarm_actual$arcs, learning = cluster_BNL(Alarm, N = i)$arcs, n = ncol(Alarm))
}
Alarm_acc1 = data.frame(NG = 1:ncol(Alarm), Accuracy = Alarm_acc1)

ggplot(Alarm_acc1, mapping = aes(NG, Accuracy)) + geom_hline(yintercept = Alarm_acc1$Accuracy[nrow(Alarm_acc1)], color = "black", size = 1, linetype = 2) + geom_point(size = 2, shape = 2) + geom_line(color = "red", size = 0.8) + labs(x = "", y = "Accuracy") + theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), axis.title = element_text(size = 25), axis.text = element_text(size = 20))

# hill-climbing algorithm

library(ggplot2)
Alarm_acc2 = numeric(ncol(Alarm))
for(i in 1:ncol(Alarm)){
  Alarm_acc2[i] = ACC(actual = Alarm_actual$arcs, learning = cluster_BNL(Alarm, N = i, algo = hc)$arcs, n = ncol(Alarm))
}
Alarm_acc2 = data.frame(NG = 1:ncol(Alarm), Accuracy = Alarm_acc2)

ggplot(Alarm_acc2, mapping = aes(NG, Accuracy)) + geom_hline(yintercept = Alarm_acc2$Accuracy[nrow(Alarm_acc2)], color = "black", size = 1, linetype = 2) + geom_point(size = 2, shape = 1) + geom_line(color = "blue", size = 0.8) + labs(x = "The Number of Clusters", y = "Accuracy") + theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), axis.title = element_text(size = 25), axis.text = element_text(size = 20))


# Plot the network configuration (Use the fast version)
graphviz.compare(Alarm_actual, cluster_BNL(Alarm, N = which.max(Alarm_acc1$Accuracy)), gs(Alarm), diff.args = list(fp.col = "red", fp.lty = "dotted",fp.lwd = 3, fn.col = "blue", fn.lty = "longdash", fn.lwd = 1.5), layout = "dot", shape = "ellipse")
