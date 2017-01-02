############### Packages Used ###############

install.packages("gridExtra")
install.packages("XLConnect")
install.packages("wavelets")
install.packages("wmtsa")
install.packages("binhf")
install.packages("ggplot2")
install.packages("TSA")
install.packages("fitdistrplus")
install.packages("logspline")
install.packages("actuar")

require(gridExtra)
require(XLConnect)
require(wavelets)
require(wmtsa)
require(binhf)
require(ggplot2)
require(TSA)
require(fitdistrplus)
require(logspline)
require(actuar)



############### 1 Introduction ###############


#### Loading Nitrate Concentration Data into R ####
options(java.parameters = "-Xmx4g" )
wb = loadWorkbook("nitratedata.xlsx")
df = readWorksheet(wb,sheet = "Sheet1",header =TRUE)
nitrate.data = data.frame(Date = df$Date, 
                          Nitrate = as.numeric(df$Nitrate), 
                          I = as.numeric(df$Date), 
                          Index = round(((as.numeric(df$Date) - 
                                            1212555600 + 86400)/86400)-4),0)
n.d = data.frame(nitrate = log(nitrate.data$Nitrate), 
                 date = nitrate.data$Date, 
                 index = nitrate.data$Index)
n.d = data.frame(nitrate = n.d$nitrate[n.d$index >= 1 & n.d$index <= 602], 
                 date = n.d$date[n.d$index >= 1 & n.d$index <= 602], 
                 index = n.d$index[n.d$index >= 1 & n.d$index <= 602])


#### Reproduce Figure 1 ####
png("Nitrate_Concentration_Time_Series.png", 
    units = "in", height = 4, width = 8, res = 300)
par(mfrow=c(1,1),
    oma = c(1, 0, 0, 0),
    mar = c(2, 3, 2, 1.5) + 0.2,
    mgp = c(2, 1, 0))
plot(n.d$index, exp(n.d$nitrate), type = "l",
     ylab = "Nitrate Concentration",
     xlab = "",
     xaxt='n')
axis(side = 1, at = c(1, 25 * 7 + 2, 51 * 7 + 2, 77 * 7 + 3), 
     labels = c("Jun 2008", "Dec 2008", "Jun 2009", "Dec 2009"))
dev.off()



############### 5 Simulation Studies ###############


#### Functions to Acquire Wavelet Coefficients ####

### Haar Wavelet Coefficients Function (for layer 4) ###
jump.dector.MODWT.haar = function(simu.data, layer, index){
  
  n = length(simu.data)
  TY.reflection = 0
  l = 0
  threshold = 0
  my.jump = 0
  
  j_n = log(n,2) - log(n / ((log(n))^2), 2)
  sim.data.reflection = extend.series(simu.data, method = "reflection", 
                                      length = "double")
  sim.wt.reflection = wavMODWT(sim.data.reflection, n.levels = 6, 
                               wavelet = "haar")
  c.char = "sim.wt.reflection$data$d"
  command.char = as.character(paste(c.char,layer, sep = ""))
  TY.reflection = eval(parse(text = command.char))
  
  if (layer == 4) TY.reflection = shift(as.numeric(TY.reflection), places = 7, 
                                        dir = "left")
  l = length(TY.reflection) / 2
  threshold = median(abs(TY.reflection)) * sqrt(2 * log(n)) / 0.6745
  TY = abs(TY.reflection[1:l])
  
  my.frame = data.frame(index = c(1:length(TY)), abs.TY = abs(TY))
  my.jumpTY = my.frame$index[my.frame$abs.TY > threshold]
  
  if (length(my.jumpTY) == 0) my.jump = "NA"
  if (length(my.jumpTY) > 0) my.jump = my.jumpTY
  
  my.jump.fil = rep(NA, length(my.jump))
  if (my.jump[1] == "NA"){my.jump.fil = "NA"}
  if (my.jump[1] != "NA"){
    my.jump.fil[1] = my.jump[1] - 1 + 
      which.max(TY[my.jump[1]:(my.jump[1] + index - 1)])
    i = 1
    j = 1
    for (i in 1:length(my.jump)){
      for (j in 1:length(my.jump) - 1){
        if (is.na(my.jump.fil[i])) {break}
        if ((my.jump[j + 1] - my.jump.fil[i]) >= index) {
          my.jump.fil[i + 1] = my.jump[j + 1] - 1 + 
            which.max(TY[my.jump[j + 1]:(my.jump[j + 1] + index - 1)])
          break
        }
      }
    }
  }
  my.jump = my.jump.fil[!is.na(my.jump.fil)]
  
  if (my.jump[1] == "NA") { 
    my.jump.size = 0 
    jump.adjusted.data = simu.data                     
  }
  if (my.jump[1] != "NA") { 
    my.jump.size = rep(NA,length(my.jump))
    for (j in 1:length(my.jump)){
      my.jump.size[j] = mean(simu.data[my.jump[j]:min((my.jump[j] +
                                                         index / 2 - 1), 
                                                      length(simu.data))]) -
        mean(simu.data[max((my.jump[j] - index / 2),0):(my.jump[j] - 1)])
      my.jump.size[is.nan(my.jump.size)] = 0
    }
    my.jump.part = rep(0,n)
    for (k in 1:length(my.jump)){
      for (i in my.jump[k]:n){
        my.jump.part[i] = my.jump.part[i] + my.jump.size[k]
      }
    }
    jump.adjusted.data = simu.data - my.jump.part
  }
  return(list(detect.jump.loca = my.jump, detect.jump.size = my.jump.size, 
              jump.adj.data = jump.adjusted.data, range = index))
}

### s8 Wavelet Coefficients Function (for layers 1 and 2) ###
jump.dector.MODWT.s8 = function(simu.data, layer, index){
  
  n = length(simu.data)
  TY.reflection = 0
  l = 0
  threshold = 0
  my.jump = 0
  
  sim.data.reflection = extend.series(simu.data,method="reflection", 
                                      length="double")
  sim.wt.reflection = wavMODWT(sim.data.reflection, n.levels=6,wavelet="s8")
  c.char = "sim.wt.reflection$data$d"
  command.char = as.character(paste(c.char,layer, sep = ""))
  TY.reflection = eval(parse(text = command.char))
  
  if (layer == 1) TY.reflection = shift(as.numeric(TY.reflection), 
                                        places = 4, dir = "left")
  if (layer == 2) TY.reflection = shift(as.numeric(TY.reflection), 
                                        places = 11, dir = "left")
  
  l = length(TY.reflection) / 2
  threshold = median(abs(TY.reflection)) * sqrt(2 * log(n)) / 0.6745
  TY=abs(TY.reflection[1:l])
  my.frame = data.frame(index = c(1:length(TY)), abs.TY = abs(TY))
  my.jumpTY = my.frame$index[my.frame$abs.TY > threshold]
  
  if (length(my.jumpTY) == 0) my.jump = "NA"
  if (length(my.jumpTY) > 0) my.jump = my.jumpTY
  
  my.jump.fil = rep(NA,length(my.jump))
  if (my.jump[1] == "NA"){my.jump.fil = "NA"}
  if (my.jump[1] != "NA"){
    my.jump.fil[1] = my.jump[1] - 1 + 
      which.max(TY[my.jump[1]:(my.jump[1] + index - 1)])
    i = 1
    j = 1
    for (i in 1:length(my.jump)){
      for (j in 1:length(my.jump) - 1){
        if (is.na(my.jump.fil[i])) {break}
        if ((my.jump[j + 1] - my.jump.fil[i]) >= index) {
          my.jump.fil[i + 1] = my.jump[j + 1] - 1 + 
            which.max(TY[my.jump[j + 1]:(my.jump[j + 1] + index - 1)])
          break
        }
      }
    }
  }
  my.jump = my.jump.fil[!is.na(my.jump.fil)]
  
  if (my.jump[1] == "NA") { 
    my.jump.size = 0 
    jump.adjusted.data = simu.data                     
  }
  if (my.jump[1] != "NA") { 
    my.jump.size = rep(NA,length(my.jump))
    for (j in 1:length(my.jump)){
      my.jump.size[j] = mean(simu.data[my.jump[j]:min((my.jump[j] +
                                                         index - 1), 
                                                      length(simu.data))]) -
        mean(simu.data[max((my.jump[j] - index), 0):(my.jump[j] - 1)])
      my.jump.size[is.nan(my.jump.size)] = 0
    }
    my.jump.part = rep(0,n)
    for (k in 1:length(my.jump)){
      for (i in my.jump[k]:n){
        my.jump.part[i] = my.jump.part[i] + my.jump.size[k]
      }
    }
    jump.adjusted.data = simu.data - my.jump.part
  }
  return(list(detect.jump.loca = my.jump, detect.jump.size = my.jump.size, 
              jump.adj.data = jump.adjusted.data, range = index))
}

### d4 Wavelet Coefficients function (for layers 2 and 3) ###
jump.dector.MODWT.d4 = function(simu.data, layer, index){
  
  n = length(simu.data)
  TY.reflection = 0
  l = 0
  threshold = 0
  my.jump = 0
  
  sim.data.reflection = extend.series(simu.data,method="reflection", 
                                      length="double")
  sim.wt.reflection = wavMODWT(sim.data.reflection, n.levels = 6, 
                               wavelet = "d4")
  c.char = "sim.wt.reflection$data$d"
  command.char = as.character(paste(c.char,layer, sep = ""))
  TY.reflection = eval(parse(text = command.char))
  
  if (layer == 2) TY.reflection = shift(as.numeric(TY.reflection), 
                                        places = 5, dir = "left")
  if (layer == 3) TY.reflection = shift(as.numeric(TY.reflection), 
                                        places = 12, dir = "left")
  
  l = length(TY.reflection) / 2
  threshold = median(abs(TY.reflection)) * sqrt(2 * log(n)) / 0.6745
  TY=abs(TY.reflection[1:l])
  my.frame = data.frame(index = c(1:length(TY)), abs.TY = abs(TY))
  my.jumpTY = my.frame$index[my.frame$abs.TY > threshold]
  
  if (length(my.jumpTY) == 0) my.jump = "NA"
  if (length(my.jumpTY) > 0) my.jump = my.jumpTY
  
  my.jump.fil = rep(NA,length(my.jump))
  if (my.jump[1] == "NA"){my.jump.fil = "NA"}
  if (my.jump[1] != "NA"){
    my.jump.fil[1] = my.jump[1] - 1 
    + which.max(TY[my.jump[1]:(my.jump[1] + index - 1)])
    i = 1
    j = 1
    for (i in 1:length(my.jump)){
      for (j in 1:length(my.jump) - 1){
        if (is.na(my.jump.fil[i])) {break}
        if ((my.jump[j + 1] - my.jump.fil[i]) >= index) {
          my.jump.fil[i + 1] = my.jump[j + 1] - 1 +
            which.max(TY[my.jump[j + 1]:(my.jump[j + 1] + index - 1)])
          break
        }
      }
    }
  }
  my.jump = my.jump.fil[!is.na(my.jump.fil)]
  
  if (my.jump[1] == "NA") { 
    my.jump.size = 0 
    jump.adjusted.data = simu.data                     
  }
  if (my.jump[1] != "NA") { 
    my.jump.size = rep(NA,length(my.jump))
    for (j in 1:length(my.jump)){
      my.jump.size[j] = mean(simu.data[my.jump[j]:min((my.jump[j] +
                                                         index - 1), 
                                                      length(simu.data))]) - 
        mean(simu.data[max((my.jump[j] - index),0):(my.jump[j] - 1)])
      my.jump.size[is.nan(my.jump.size)] = 0
    }
    my.jump.part = rep(0,n)
    for (k in 1:length(my.jump)){
      for (i in my.jump[k]:n){
        my.jump.part[i] = my.jump.part[i] + my.jump.size[k]
      }
    }
    jump.adjusted.data = simu.data - my.jump.part
  }
  return(list(detect.jump.loca = my.jump, detect.jump.size = my.jump.size, 
              jump.adj.data = jump.adjusted.data, range = index))
}


#### Simulation Data Generator from model (5.1) in the paper ####
data.generator = function(v, j.size, w, noise.level){

  j.n = w
  
  ### Continuous Part Generator ###
  conti.b1 = rnorm(v, 0, sd = sqrt(1 / v))
  conti.b2 = rnorm(v, 0, sd = sqrt(1 / v))
  conti.b3 = -0.62 * conti.b1 + sqrt(1 - 0.62 ^ 2) * conti.b2

  ### Euler Scheme ###
  conti.log.sigma.sq = rep(0, v)
  conti.log.sigma.sq.0 = rnorm(1, mean = -6.802, sd = sqrt(0.3125))
  for (i in 0:v - 1) {
    if (i == 0) conti.log.sigma.sq[i + 1] = conti.log.sigma.sq.0 
    + (-0.6802 - 0.10 * conti.log.sigma.sq.0) / v 
    + 0.25 * conti.b1[i + 1]
    if (i >= 1) conti.log.sigma.sq[i + 1] = conti.log.sigma.sq[i] 
    + (-0.6802 - 0.10 * conti.log.sigma.sq[i]) / v 
    + 0.25 * conti.b1[i + 1]
  }
  conti.sigma.sq = exp(conti.log.sigma.sq)
  conti.sigma = sqrt(exp(conti.log.sigma.sq))
  conti.diff = conti.sigma * conti.b3
  conti.part = rep(NA, v)
  for (i in 1:v) {
    conti.part[i] = sum(conti.diff[1:i])
  }
  
  ### Jump Part Generator ###
  if (w >= 2){
    jump.loca = sort(round(runif(w, min = 0, max = 1) * v, digits = 0))
    while (min(diff(jump.loca)) <= 16) {
      jump.loca = sort(round(runif(w, min = 0, max = 1) * v, digits = 0))
    }
  }
  if (w <= 1){
    jump.loca = round(runif(w, min = 0, max = 1) * v, digits = 0)
  }
  jump.size = j.size
  jump.part = rep(0,v)
  for (j in 1:length(jump.loca)){
    for (i in jump.loca[j]:v){
      jump.part[i] = jump.part[i] + jump.size[j]
    }
  }
  
  ### Noise Part Generator ###
  noise.part = rnorm(v, mean = 0, sd = noise.level)
  
  ### Simulation Data ###
  simu.data = conti.part + jump.part + noise.part
  simu.data.nojump = conti.part + noise.part
  simu.data.conti = conti.part
  
  return(list(simu.data = simu.data, simu.data.nojump = simu.data.nojump, 
              simu.data.conti = simu.data.conti, jump.loca = jump.loca, 
              sigma.sq = conti.sigma.sq))
}


#### Integrated Volatilities Analysis Functions ####

### RV ###
RV = function(jump.adjusted.data){
  data = jump.adjusted.data
  RV = sum((diff(data) ^ 2))
  RV
}

### ASRV ###
ASRV.calculator = function(jump.adjusted.data, K){
  data = jump.adjusted.data
  n = length(data)
  sum.ASRV = 0
  for (i in 1:(n - K)){
    sum.ASRV = sum.ASRV + (data[i + K] - data[i]) ^ 2  
  }
  ASRV = sum.ASRV / K
  return(ASRV)
}

### JTSRV ###
JTSRV = function(jump.adjusted.data,K){
  data = jump.adjusted.data
  JTSRV = ASRV.calculator(data,K) - 1 / K * ASRV.calculator(data,1)
  JTSRV
}

### JMSRV ###
JMSRV = function(jump.adjusted.data,M,C){ 
  data = jump.adjusted.data 
  n = length(data)
  v = vector(mode="numeric", length = 1)
  a = vector(mode="numeric", length = M)
  theta = 0
  for (m in 1:M){
    a[m] = 12*(m + C) * (m - M / 2 - 1 / 2) / (M * (M ^ 2 - 1))
  }
  zeta = (M + C) * (C + 1) / ((n + 1) * (M - 1))
  for (m in 1:M){
    Km = m + C
    v = v + a[m] * ASRV.calculator(data, Km)
  }
  theta = v + zeta*(ASRV.calculator(data, C + 1) -
                      ASRV.calculator(data, C + M))
  JMSRV = theta
  return(JMSRV)
}

### RBPV ###
RBPV.vol = function(simu.data){
  data = simu.data
  n = length(data)
  RBPV = 0
  i = 3
  for (i in 3:n){
    RBPV = RBPV + 1 / (0.79788) ^ 2 * abs(data[i] - data[i - 1]) * 
      abs(data[i - 1] - data[i - 2])
  }
  return(RBPV)
}

RBPV.vol.subsamp = function(simu.data, K){
  data = simu.data
  n = length(data)
  RBPV = 0
  i = 2 * K + 1
  for (i in (2 * K + 1):n){
    RBPV = RBPV + 1 / (0.79788) ^ 2 / K * abs(data[i]-data[i - K]) *
      abs(data[i - K] - data[i - 2 * K])
  }
  return(RBPV)
}


#### Reproduce of Figure 2 and Figure 3 Contents ####
set.seed(54321)

RV.MSE = c()
RV.adj.MSE = c()
RBPV.MSE = c()
RBPV.5.MSE = c()
RBPV.15.MSE = c()
RBPV.jump.MSE = c()
RBPV.5.jump.MSE = c()
RBPV.15.jump.MSE = c()
JMSRV.MSE = c()
#JTSRV.MSE = c()
JV.MSE = c()
jump.exact.pct = c()
s1.JMSRV.MSE = c()
s1.JV.MSE = c()
s1.jump.exact.pct = c()
s2.JMSRV.MSE = c()
s2.JV.MSE = c()
s2.jump.exact.pct = c()
d2.JMSRV.MSE = c()
d2.JV.MSE = c()
d2.jump.exact.pct = c()
d3.JMSRV.MSE = c()
d3.JV.MSE = c()
d3.jump.exact.pct = c()

j = 1
i = 1
t = 20
m = 3000
max.noise = 0.005
for (j in 1:t) {
  RV.record = c()
  RV.adj.record = c()
  RBPV.record = c()
  RBPV.5.record = c()
  RBPV.15.record = c()
  RBPV.jump.record = c()
  JMSRV.record = c()
  JV.MSE.record = c()
  jump.num.detect = c()
  jump.acc = c()
  s1.JMSRV.record = c()
  s1.JV.MSE.record = c()
  s1.jump.num.detect = c()
  s1.jump.acc = c()
  s2.JMSRV.record = c()
  s2.JV.MSE.record = c()
  s2.jump.num.detect = c()
  s2.jump.acc = c()
  d2.JMSRV.record = c()
  d2.JV.MSE.record = c()
  d2.jump.num.detect = c()
  d2.jump.acc = c()
  d3.JMSRV.record = c()
  d3.JV.MSE.record = c()
  d3.jump.num.detect = c()
  d3.jump.acc = c()
  noise.level = max.noise * (j - 1) / (t - 1)
  
  for (i in 1: m){
    w = round(runif(1, -0.5, 6.5))
    
    if (w != 0) {
      j.s = rnorm(w, 0.03, 0.01)
      while (min(j.s) < 0.025) {
        j.s = rnorm(w, 0.03, 0.01)
      }
      r.n = runif(w)
      for (k in 1: w){
        if (r.n[k] < 0.5) j.s[k] = -1 * j.s[k]
      }
    }
    
    if (w == 0){
      j.s = 0
      w = 1
    }
    
    n = 672
    result = data.generator(n, j.s, w, noise.level)
    if (j.s[1] == 0 & w == 1){w = 0}
    
    RV.record[i] = (RV(result$simu.data) - 
                      sum(result$sigma.sq) / n) ^ 2
    RBPV.record[i] = (RBPV.vol(result$simu.data) - 
                        sum(result$sigma.sq) / n) ^ 2
    RBPV.5.record[i] = (RBPV.vol.subsamp(result$simu.data, 5) - 
                          sum(result$sigma.sq) / n) ^ 2
    RBPV.15.record[i] = (RBPV.vol.subsamp(result$simu.data, 15) - 
                           sum(result$sigma.sq) / n) ^ 2
    RBPV.jump.record[i] = (RV(result$simu.data) -
                             RBPV.vol.subsamp(result$simu.data, 1) - 
                             sum(j.s ^ 2)) ^ 2
    
    Info = jump.dector.MODWT.haar(result$simu.data, 4, 16)
    RV.adj.record[i] = (RV(Info$jump.adj.data) - 
                          sum(result$sigma.sq) / n) ^ 2
    JMSRV.record[i] = (JMSRV(Info$jump.adj.data, 26, 26) - 
                         sum(result$sigma.sq) / n) ^ 2 
    JV.MSE.record[i] = (sum(Info$detect.jump.size ^ 2) - 
                          sum(j.s ^ 2)) ^ 2
    if (Info$detect.jump.loca[1] != "NA") {
      jump.num.detect[i] = length(Info$detect.jump.loca)
      }
    if (Info$detect.jump.loca[1] =="NA") {jump.num.detect[i] = 0}
    if (jump.num.detect[i] == w) {jump.acc[i] = 1}
    if (jump.num.detect[i] != w) {jump.acc[i] = 0}
    
    Info.s1 = jump.dector.MODWT.s8(result$simu.data, 1, 16)
    s1.JMSRV.record[i] = (JMSRV(Info.s1$jump.adj.data, 26, 26) - 
                            sum(result$sigma.sq) / n) ^ 2 
    s1.JV.MSE.record[i] = (sum(Info.s1$detect.jump.size^2) - 
                             sum(j.s ^ 2)) ^ 2
    if (Info.s1$detect.jump.loca[1] != "NA") {
      s1.jump.num.detect[i] = length(Info.s1$detect.jump.loca)
      }
    if (Info.s1$detect.jump.loca[1] =="NA") {s1.jump.num.detect[i] = 0}
    if (s1.jump.num.detect[i] == w) {s1.jump.acc[i] = 1}
    if (s1.jump.num.detect[i] != w) {s1.jump.acc[i] = 0}
    
    Info.s2 = jump.dector.MODWT.s8(result$simu.data, 2, 16)
    s2.JMSRV.record[i] = (JMSRV(Info.s2$jump.adj.data, 26, 26) - 
                            sum(result$sigma.sq) / n) ^ 2 
    s2.JV.MSE.record[i] = (sum(Info.s2$detect.jump.size ^ 2) - 
                             sum(j.s ^ 2)) ^ 2
    if (Info.s2$detect.jump.loca[1] != "NA") {
      s2.jump.num.detect[i] = length(Info.s2$detect.jump.loca)
      }
    if (Info.s2$detect.jump.loca[1] =="NA") {s2.jump.num.detect[i] = 0}
    if (s2.jump.num.detect[i] == w) {s2.jump.acc[i] = 1}
    if (s2.jump.num.detect[i] != w) {s2.jump.acc[i] = 0}
    
    Info.d2 = jump.dector.MODWT.d4(result$simu.data, 2, 16)
    d2.JMSRV.record[i] = (JMSRV(Info.d2$jump.adj.data, 26, 26) - 
                            sum(result$sigma.sq) / n) ^ 2 
    d2.JV.MSE.record[i] = (sum(Info.d2$detect.jump.size ^ 2) - 
                             sum(j.s ^ 2)) ^ 2
    if (Info.d2$detect.jump.loca[1] != "NA") {
      d2.jump.num.detect[i] = length(Info.d2$detect.jump.loca)
      }
    if (Info.d2$detect.jump.loca[1] =="NA") {d2.jump.num.detect[i] = 0}
    if (d2.jump.num.detect[i] == w) {d2.jump.acc[i] = 1}
    if (d2.jump.num.detect[i] != w) {d2.jump.acc[i] = 0}
    
    Info.d3 = jump.dector.MODWT.d4(result$simu.data, 3, 16)
    d3.JMSRV.record[i] = (JMSRV(Info.d3$jump.adj.data, 26, 26) - 
                            sum(result$sigma.sq) / n) ^ 2 
    d3.JV.MSE.record[i] = (sum(Info.d3$detect.jump.size ^ 2) - 
                             sum(j.s ^ 2)) ^ 2
    if (Info.d3$detect.jump.loca[1] != "NA") {
      d3.jump.num.detect[i] = length(Info.d3$detect.jump.loca)
      }
    if (Info.d3$detect.jump.loca[1] =="NA") {d3.jump.num.detect[i] = 0}
    if (d3.jump.num.detect[i] == w) {d3.jump.acc[i] = 1}
    if (d3.jump.num.detect[i] != w) {d3.jump.acc[i] = 0}
  }
  
  RV.MSE[j] = sum(RV.record) / m
  RV.adj.MSE[j] = sum(RV.adj.record) / m
  RBPV.MSE[j] = sum(RBPV.record) / m
  RBPV.5.MSE[j] = sum(RBPV.5.record) / m
  RBPV.15.MSE[j] = sum(RBPV.15.record) / m
  RBPV.jump.MSE[j] = sum(RBPV.jump.record) / m
  JMSRV.MSE[j] = sum(JMSRV.record) / m
  JV.MSE[j] = sum(JV.MSE.record) / m
  jump.exact.pct[j] = sum(jump.acc) / m * 100
  s1.JMSRV.MSE[j] = sum(s1.JMSRV.record) / m
  s1.JV.MSE[j] = sum(s1.JV.MSE.record) / m
  s1.jump.exact.pct[j] = sum(s1.jump.acc) / m * 100
  s2.JMSRV.MSE[j] = sum(s2.JMSRV.record) / m
  s2.JV.MSE[j] = sum(s2.JV.MSE.record) / m
  s2.jump.exact.pct[j] = sum(s2.jump.acc) / m * 100
  d2.JMSRV.MSE[j] = sum(d2.JMSRV.record) / m
  d2.JV.MSE[j] = sum(d2.JV.MSE.record) / m
  d2.jump.exact.pct[j] = sum(d2.jump.acc) / m * 100
  d3.JMSRV.MSE[j] = sum(d3.JMSRV.record) / m
  d3.JV.MSE[j] = sum(d3.JV.MSE.record) / m
  d3.jump.exact.pct[j] = sum(d3.jump.acc) / m * 100
}

### Jump Number Accuracy ###
detection.result.2 = c()
jump.num.2 = c()
jump.acc.2 = c()
detection.result.2 = c()
s1.detection.result.2 = c()
s1.jump.num.2 = c()
s1.jump.acc.2 = c()
s1.detection.result.2 = c()
s2.detection.result.2 = c()
s2.jump.num.2 = c()
s2.jump.acc.2 = c()
s2.detection.result.2 = c()
d2.detection.result.2 = c()
d2.jump.num.2 = c()
d2.jump.acc.2 = c()
d2.detection.result.2 = c()
d3.detection.result.2 = c()
d3.jump.num.2 = c()
d3.jump.acc.2 = c()
d3.detection.result.2 = c()

N = 3000
for (k in 1:20){
  j.size = 0.005 * (k)
  n = 672
  
  for (i in 1:N) {
    
    w = round(runif(1, 0.500001, 6.5))
    j.s = rep(0, w)
    r.n = runif(w)
    for (l in 1:w){
      if (r.n[l] < 0.5) j.s[l] = -j.size
      if (r.n[l] >= 0.5) j.s[l] = j.size
    }
    noise.level = max(rnorm(1, 0.0025, 0.001), 0)
    result = data.generator(n, j.s, w, noise.level)
    Info = jump.dector.MODWT.haar(result$simu.data, 4, 16)
    if (Info$detect.jump.loca[1] != "NA") {
      detection.result.2[i] = length(Info$detect.jump.loca)}
    if (Info$detect.jump.loca[1] =="NA") {detection.result.2[i] = 0
    }
    if (detection.result.2[i] == w) {jump.acc.2[i] = 1}
    if (detection.result.2[i] != w) {jump.acc.2[i] = 0}
    
    Info.s1 = jump.dector.MODWT.s8(result$simu.data, 1, 16)
    if (Info.s1$detect.jump.loca[1] != "NA") {
      s1.detection.result.2[i] = length(Info.s1$detect.jump.loca)
      }
    if (Info.s1$detect.jump.loca[1] =="NA") {s1.detection.result.2[i] = 0}
    if (s1.detection.result.2[i] == w) {s1.jump.acc.2[i] = 1}
    if (s1.detection.result.2[i] != w) {s1.jump.acc.2[i] = 0}
    
    Info.s2 = jump.dector.MODWT.s8(result$simu.data, 2, 16)
    if (Info.s2$detect.jump.loca[1] != "NA") {
      s2.detection.result.2[i] = length(Info.s2$detect.jump.loca)
      }
    if (Info.s2$detect.jump.loca[1] =="NA") {s2.detection.result.2[i] = 0}
    if (s2.detection.result.2[i] == w) {s2.jump.acc.2[i] = 1}
    if (s2.detection.result.2[i] != w) {s2.jump.acc.2[i] = 0}
    
    Info.d2 = jump.dector.MODWT.d4(result$simu.data, 2, 16)
    if (Info.d2$detect.jump.loca[1] != "NA") {
      d2.detection.result.2[i] = length(Info.d2$detect.jump.loca)
      }
    if (Info.d2$detect.jump.loca[1] =="NA") {d2.detection.result.2[i] = 0}
    if (d2.detection.result.2[i] == w) {d2.jump.acc.2[i] = 1}
    if (d2.detection.result.2[i] != w) {d2.jump.acc.2[i] = 0}
    
    Info.d3 = jump.dector.MODWT.d4(result$simu.data, 3, 16)
    if (Info.d3$detect.jump.loca[1] != "NA") {
      d3.detection.result.2[i] = length(Info.d3$detect.jump.loca)
      }
    if (Info.d3$detect.jump.loca[1] =="NA") {d3.detection.result.2[i] = 0}
    if (d3.detection.result.2[i] == w) {d3.jump.acc.2[i] = 1}
    if (d3.detection.result.2[i] != w) {d3.jump.acc.2[i] = 0}
  }
  
  jump.num.2[k] = sum(jump.acc.2) / N * 100
  s1.jump.num.2[k] = sum(s1.jump.acc.2) / N * 100
  s2.jump.num.2[k] = sum(s2.jump.acc.2) / N * 100
  d2.jump.num.2[k] = sum(d2.jump.acc.2) / N * 100
  d3.jump.num.2[k] = sum(d3.jump.acc.2) / N * 100
}


#### Reproduce Figure 2 ####
png("Wavelet_Method_Evaluation.png", 
    units = "in", height = 5.1, width = 9.7, res = 300)

par(mfrow = c(1,2), oma = c(1.5, 0, 0, 0), 
    mar = c(3, 3, 3, 1.5) + 0.2, mgp = c(2, 1, 0))

x = seq(0,max.noise,max.noise/19)
plot(x, JV.MSE, type = "l", lty = 1, 
     col = "black", 
     cex.lab = 0.85, cex.axis = 0.85, 
     cex.main = 1, cex.sub = 0.85,
     ylim = c(0, 0.06 / 10000), ylab = "MSE", 
     xlab = "Noise Level", main = "(a) Jump Variation Accuracy")
lines(x, RBPV.jump.MSE, type = "l", lty = 3, col = "grey")

plot(x, JMSRV.MSE, type = "l", lty = 1, 
     col = "black",
     cex.lab = 0.85, cex.axis = 0.85, 
     cex.main = 1, cex.sub = 0.85,
     ylim = c(0, 0.4 / 10000),
     ylab = "MSE", xlab = "Noise Level", 
     main = "(b) Integrated Volatility Accuracy")
lines(x, RV.MSE, type = "l", lty = 2, col = "grey50")
lines(x, RV.adj.MSE, type = "l", lty = 5, col = "grey80")
lines(x, RBPV.MSE, type = "l", lty =3, col = "grey")
lines(x, RBPV.5.MSE, type = "l", lty = 4, col = "grey50")
lines(x, RBPV.15.MSE, type = "l", lty = 6, col = "grey80")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
    mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Haar Level 4", "RV", "RV.adj", "RBPV", 
                   "RBPV-5", "RBPV-15"), 
       col = c("black", "grey50", "grey80", "grey", "grey50", "grey80"),
       xpd = FALSE, 
       horiz = TRUE, 
       inset = c(0, 0),
       lty = c(1, 2, 5, 3, 4, 6),
       bty = "n",
       cex = 0.85,
       text.width = c(0.20, 0.155, 0.18, 0.155, 0.20, 0.20),
       x.intersp = 0.8)

dev.off()


#### Reproduce Figure 3 ####
png("Wavelet_Filter_Comparison.png", 
    units = "in", height = 9, width = 9, res = 300)

par(mfrow = c(2, 2), oma = c(2, 0, 0, 0), 
    mar = c(3, 3, 3, 1.5) + 0.2, mgp = c(2, 1, 0))

x = seq(0, max.noise, max.noise / 19)
plot(x, JV.MSE, type = "l", lty = 1, col = "black",
     ylim = c(0,0.08 / 10000), ylab = "MSE", 
     xlab = "Noise Level", main = "(a) Jump Variation Accuracy")
lines(x, s1.JV.MSE, type = "l", lty = 2, col = "grey50")
lines(x, s2.JV.MSE, type = "l", lty = 5, col = "grey80")
lines(x, d2.JV.MSE, type = "l", lty = 4, col = "grey50")
lines(x, d3.JV.MSE, type = "l", lty = 6, col = "grey80")
lines(x, RBPV.jump.MSE, type="l", lty = 3, col = "black")

plot(x, JMSRV.MSE, type = "l", lty = 1,  col = "black",
     ylim = c(0,0.08 / 10000),
     ylab = "MSE", xlab = "Noise Level", 
     main = "(b) Integrated Volatility Accuracy")
lines(x, s1.JMSRV.MSE, type = "l", lty = 2, col = "grey50")
lines(x, s2.JMSRV.MSE, type = "l", lty = 5, col = "grey80")
lines(x, d2.JMSRV.MSE, type = "l", lty = 4, col = "grey50")
lines(x, d3.JMSRV.MSE, type = "l", lty = 6, col = "grey80")
lines(x, RBPV.MSE, type="l", lty = 3, col = "black")

plot(x, jump.exact.pct, type = "l", col = "black", 
     ylim = c(0,100), ylab = "Percentage (%)", 
     xlab = "Noise Level", 
     main = "(c) Jump Number Accuracy vs. Noise Level")
lines(x, s1.jump.exact.pct, type = "l", lty = 2, col = "grey50")
lines(x, s2.jump.exact.pct, type = "l", lty = 5, col = "grey80")
lines(x, d2.jump.exact.pct, type = "l", lty = 4, col = "grey50")
lines(x, d3.jump.exact.pct, type = "l", lty = 6, col = "grey80")

u = seq(0.005, 0.1, 0.005)
plot(u, jump.num.2, type = "l", col = "black",
     ylim = c(0, 100), ylab="Percentage (%)", 
     xlab = "Jump Magnitude", 
     main = "(d) Jump Number Accuracy vs. Jump Magnitude")
lines(u, s1.jump.num.2, type = "l", lty = 2, col = "grey50")
lines(u, s2.jump.num.2, type = "l", lty = 5, col = "grey80")
lines(u, d2.jump.num.2, type = "l", lty = 4, col = "grey50")
lines(u, d3.jump.num.2, type = "l", lty = 6, col = "grey80")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
    mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Haar Level 4", "s8 Level 1", "s8 Level 2", 
                   "d4 Level 2", "d4 Level 3", "RBPV"), 
       col = c("black", "grey50", "grey80", "grey50", "grey80", "grey"),
       xpd = FALSE, 
       horiz = TRUE, 
       inset = c(0, 0),
       lty = c(1, 2, 5, 4, 6, 3),
       bty = "n",
       cex = 1,
       text.width = c(0.24, rep(0.225, 5)),
       x.intersp = 0.8,
       xjust = 0.5)

dev.off()



############### 6 Nitrate Concentration Analysis ###############

########## Reproduce 6.1 Process Variation Analysis ##########

#### Reproduce Process Variation Analysis Information ####

### Weekly jump adjusted data vectors and jump information Estimation ###

for (i in 1 : 86){
  data = as.vector(n.d[n.d$index > 7*(i - 1) & n.d$index <= 7 * i, 1])
  assign(paste("data.week", i, sep = "."), data) 
  data23 = as.vector(
    n.d[max(min(which(n.d$index > 7 * (i - 1) & n.d$index <= 7 * i) - 16), 0): 
          min(max(which(n.d$index > 7 * (i - 1) & n.d$index <= 7 * i) + 16), 
              length(n.d$index)), 1])
  Info = jump.dector.MODWT.haar(data23,layer = 4, index = 16)
  if (i == 1) {
    jump.adj.data = Info$jump.adj.data[1 : (length(Info$jump.adj.data) - 16)]
    }
  if (i > 1 & i <= 86) {
    jump.adj.data = Info$jump.adj.data[17 : 
                                         (length(Info$jump.adj.data) - 16)] + 
      sum(Info$detect.jump.size[which(Info$detect.jump.loca <= 16)])
    }
  assign(paste("jump.adj.week", i, sep = "."), jump.adj.data)
  
  if (i > 1 & i <= 86 & Info$detect.jump.loca[1] != "NA") {
    Info23 = data.frame(detect.jump.loca = Info$detect.jump.loca - 16, 
                        detect.jump.size = Info$detect.jump.size)
    Info23 = Info23[which(Info23$detect.jump.loca >= 1 & 
                            Info23$detect.jump.loca <= (length(
                              Info$jump.adj.data) - 32)), ]
  }
  if (i == 1 & Info$detect.jump.loca[1] != "NA") {
    Info23 = data.frame(detect.jump.loca = Info$detect.jump.loca, 
                        detect.jump.size = Info$detect.jump.size)
    Info23 = Info23[which(Info23$detect.jump.loca >= 1 & 
                            Info23$detect.jump.loca <= (length(
                              Info$jump.adj.data) - 16)), ]
  }
 if (Info$detect.jump.loca[1] == "NA"){
    Info23 = data.frame(detect.jump.loca = Info$detect.jump.loca, 
                        detect.jump.size = Info$detect.jump.size)
  }
  if (nrow(Info23)==0){
    Info23[1,1] = "NA"
    Info23[1,2] = 0
  }
  assign(paste("jump.loca.week", i, sep = "."), 
         as.vector(Info23$detect.jump.loca))
  assign(paste("jump.posi.loca.week", i, sep = "."), 
         as.vector(Info23$detect.jump.loca[which(
           Info23$detect.jump.size > 0)])) 
  assign(paste("jump.nega.loca.week", i, sep = "."), 
         as.vector(Info23$detect.jump.loca[which(
           Info23$detect.jump.size < 0)]))
  assign(paste("jump.size.week", i, sep = "."), Info23$detect.jump.size)
  assign(paste("jump.posi.size.week", i, sep = "."), 
         Info23$detect.jump.size[Info23$detect.jump.size > 0])
  assign(paste("jump.nega.size.week", i, sep = "."), 
         Info23$detect.jump.size[Info23$detect.jump.size < 0])
  assign(paste("length.week", i, sep = "."), length(jump.adj.data)) 
}

### Jump Variation Generator ###
jump.var.generator.1 = function(jump.size){
  est.jump.var = sum(jump.size ^ 2)
  return(est.jump.var)  
}

### Weekly data length vectors ###
all.length.week = c()
for (i in 1 : 86){
  char = as.character(paste("length.week",i, sep = "."))
  all.length.week = c(all.length.week, eval(parse(text = char)))
}

### Weekly data cut points vector ###
week.cut = c()
for (i in 1 : 86){
  week.cut[i] = sum(all.length.week[1 : i])
}

### Jump adjusted Nitrate concentration data vector for the whole period ###
all.week.jump.adj.data = c()
for (i in 1 : 86){
  char = as.character(paste("jump.adj.week",i, sep = "."))
  all.week.jump.adj.data = c(all.week.jump.adj.data, 
                             eval(parse(text = char)))
}

### Jump locations vector ###
all.week.jump.loca.data = c()
test = c()
for (i in 2:86){
  char = as.character(paste("jump.loca.week", i, sep = "."))
  if (eval(parse(text = char))[1] == "NA") {week.jump.loca.data.i = "NA"}
  if ((eval(parse(text = char))[1] != "NA")) {
    week.jump.loca.data.i = eval(parse(text = char)) + week.cut[i - 1]
    }
  all.week.jump.loca.data = c(all.week.jump.loca.data, week.jump.loca.data.i)
}
all.week.jump.loca.data = c(jump.loca.week.1, 
                            as.numeric(all.week.jump.loca.data[
                              all.week.jump.loca.data != "NA"]))

### Jump Sizes vector ###
all.week.jump.size.data = c()
for (i in 2 : 86){
  char = as.character(paste("jump.size.week",i, sep = "."))
  all.week.jump.size.data = c(all.week.jump.size.data, 
                              eval(parse(text = char)))
}
all.week.jump.size.data = c(jump.size.week.1, 
                            as.numeric(all.week.jump.size.data[
                              all.week.jump.size.data != 0]))
jump.info23 = data.frame(jump.loca = all.week.jump.loca.data, 
                         jump.size = all.week.jump.size.data)

### Weekly jump numbers ###
week.jump.num = c()
week.posi.jump.num = c()
week.nega.jump.num = c()
for (i in 1:86){
  char = as.character(paste("jump.loca.week",i, sep = "."))
  if (eval(parse(text = char))[1] == "NA") {week.jump.num[i] = 0}
  if (eval(parse(text = char))[1] != "NA") {
    week.jump.num[i] = length(eval(parse(text = char)))
    }
  char = as.character(paste("jump.posi.size.week",i, sep = "."))
  if (is.na(eval(parse(text = char))[1])) {week.posi.jump.num[i] = 0}
  if (!is.na(eval(parse(text = char))[1])) {
    week.posi.jump.num[i] = length(eval(parse(text = char)))
    }
  char = as.character(paste("jump.nega.size.week",i, sep = "."))
  if (is.na(eval(parse(text = char))[1])) {week.nega.jump.num[i] = 0}
  if (!is.na(eval(parse(text = char))[1])) {
    week.nega.jump.num[i] = length(eval(parse(text = char)))
    } 
}

### Weekly jump variation ###
week.jump.var = c()
week.posi.jump.var = c()
week.nega.jump.var = c()
for (i in 1:86){
  char = as.character(paste("jump.size.week",i, sep = "."))
  if (is.na(eval(parse(text = char))[1])) {week.jump.var[i] = 0}
  if (!is.na(eval(parse(text = char))[1])) {
    week.jump.var[i] = jump.var.generator.1(eval(parse(text = char)))
    } 
  char = as.character(paste("jump.posi.size.week",i, sep = "."))
  if (is.na(eval(parse(text = char))[1])) {week.posi.jump.var[i] = 0}
  if (!is.na(eval(parse(text = char))[1])) {
    week.posi.jump.var[i] = jump.var.generator.1(eval(parse(text = char)))
    } 
  char = as.character(paste("jump.nega.size.week",i, sep = "."))
  if (is.na(eval(parse(text = char))[1])) {week.nega.jump.var[i] = 0}
  if (!is.na(eval(parse(text = char))[1])) {
    week.nega.jump.var[i] = jump.var.generator.1(eval(parse(text = char)))
    } 
}

### Integrated Volatility Estimation by JMSRV ###
ni.JMSRV = c()

for (i in 1:86){
  char = as.character(paste("jump.adj.week",i, sep = "."))
  data = eval(parse(text = char))
  ni.JMSRV[i] = JMSRV(data, 26, 26)
}


#### Reproduce Figure 4 ####
jump.number.detail = data.frame(jump.index = all.week.index,
                                jump.loca = all.week.jump.loca2.data,
                                jump.sign = sign(all.week.jump.size.data),
                                jump.size = exp(n.d$nitrate[
                                  all.week.jump.loca.data - 1]) * 
                                  (exp(all.week.jump.size.data) - 1))

p6 = ggplot(jump.number.detail, 
            aes(x = jump.loca, y = jump.size)) + 
  geom_point(aes(colour = factor(jump.sign), shape = factor(jump.sign))) + 
  scale_color_manual("Jump Category: ", 
                     labels = c("Negative Jump ", "Positive Jump "), 
                     values = c("grey70", "black")) +
  scale_shape_manual("Jump Category: ", 
                     labels = c("Negative Jump ", "Positive Jump "), 
                     values = c(19, 17)) +
  facet_wrap(~jump.index, ncol = 8, as.table = T) + 
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, color="grey60", 
             linetype="solid", size = 0.25) +
  theme(panel.background = element_rect(fill = 'grey96')) +
  geom_segment(aes(xend = jump.loca, yend = 0, 
                   color = factor(jump.sign)), size = 0.5) +
  theme(strip.background = element_rect(fill = "white")) + 
  theme(strip.text.y = element_text(angle = 180)) +
  labs(x = "Weekly  Jump  Location", y = "Jump   Size") +
  theme(plot.margin = unit(c(0.25,0.25,-0.25,0.25), "cm")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylim(round(-max(abs(jump.number.detail$jump.size))) - 0.5, 
       round(max(abs(jump.number.detail$jump.size))) + 0.5)

png("Original_Scale_Weekly_Nitrate_Jump_Detail_Haar_4.png", 
    units = "in", height = 8, width = 8, res = 300)
p6
dev.off()


#### Reproduce Figure 5 ####
png("Log_Scale_Nitrate_Summary_Haar_4.png", 
    units = "in", height = 8, width = 8, res = 300)

par(mfrow=c(2, 2), oma = c(2.5, 0, 0, 0), 
    mar = c(2, 3, 3, 1.5) + 0.2, mgp = c(2, 1, 0))

nega.num = week.nega.jump.num
posi.num = week.posi.jump.num
bound1 = max(nega.num)
bound2 = max(posi.num)
plot(posi.num,type = "h", las = 2, xaxt = 'n', yaxt = 'n', 
     ylab = "Jump Numbers", xlab = "", ylim = c(-bound1, bound2),
     main = "(a) Detected Weekly Jump Numbers")
lines(-nega.num, type = "h", las = 2, lty = 3)
axis(side = 1, at = c(1, 26, 52, 78), 
     labels = c("Jun 2008", "Dec 2008", "Jun 2009", "Dec 2009"))
axis(side = 2, at = seq(-bound1, bound2, 3), 
     labels = c(abs(seq(-bound1, bound2, 3))))

nega.var = week.nega.jump.var
posi.var = week.posi.jump.var
bound = round(max(week.nega.jump.var, posi.var) + 5, digits = -1)
bound1 = round(max(nega.var)+1, digits = 0)
bound2 = round(max(posi.var)+1, digits = 0)
plot(posi.var, type = "h", las = 2, yaxt = 'n', ylab = "Jump Variation", 
     ylim = c(-bound1, bound2), xaxt = 'n', xlab="", 
     main="(b) Estimated Weekly Jump Variation (JV)")
lines(-nega.var, type = "h", las = 2, lty = 3)
axis(side = 1, at = c(1, 26, 52, 78), 
     labels = c("Jun 2008", "Dec 2008", "Jun 2009", "Dec 2009"))
axis(side = 2, at = seq(-bound1, bound2, 4), 
     labels = abs(seq(-bound1, bound2, 4)))

plot(ni.JMSRV, type = "h", las = 2, xaxt = 'n', xlab = "", 
     ylab = "Integrated Volatility", 
     main = "(c) Estimated Weekly Integrated Volatility (IV)")
axis(side = 1, at = c(1, 26, 52, 78), 
     labels = c("Jun 2008", "Dec 2008", "Jun 2009", "Dec 2009"))

posi.ratio = (week.posi.jump.var / ni.JMSRV)
nega.ratio = (week.nega.jump.var / ni.JMSRV)
bound1 = ceiling(max(nega.ratio) + 1)
bound2 = ceiling(max(posi.ratio))
plot(posi.ratio, type = "h", las = 2, xaxt = 'n', xlab = "", 
     yaxt = 'n', ylab = "Ratio", ylim = c(-bound1, bound2), 
     main = "(d) Ratios of Estimated JV to IV")
lines(-nega.ratio, type = "h", las = 2, lty = 3)
axis(side = 1, at = c(1, 26, 52, 78), 
     labels=c("Jun 2008", "Dec 2008", "Jun 2009", "Dec 2009"))
axis(side = 2, at = seq(-bound1, bound2, 3), 
     labels = abs(seq(-bound1, bound2, 3)))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
    mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Positive Jump Related Quantity", 
                   "Negative Jump Related Quantity"), 
       xpd = FALSE, horiz = TRUE, lty = c(1, 3), bty = "n", 
       cex = 1, x.intersp = 1)

dev.off()


########## Reproduce 6.2 Jump Distribution Analysis ##########

#### Reproduce Jump Distribution Analysis Information ####

### Jump Week Index Vector ###

all.week.index = c()
for (i in 1:86){
  char = as.character(paste("jump.size.week", i, sep = "."))
  if ((eval(parse(text = char))[1] != 0)) {
    week.jump.size.data.i = eval(parse(text = char))
    if (i <= 9){char.i = as.character(paste("Week 0", i, sep = ""))}
    if (i >= 10){char.i = as.character(paste("Week", i, sep = " "))}
    all.week.index = c(all.week.index, 
                       rep(char.i, length(week.jump.size.data.i)))
  }
}

### Detected Jump Time and Quarter Information ###
date.time = data.frame(Date = df$Date, 
                       Time = format(df$Time, format="%H:%M:%S"), 
                       Value = df$Nitrate)
date.time = date.time[353 : length(date.time[, 1]), 1 : 2]
jump.date.time = date.time[all.week.jump.loca.data, ]
getSeason = function(DATES) {
  WS = as.Date("2012-10-1", format = "%Y-%m-%d") # Q1
  SE = as.Date("2012-1-1",  format = "%Y-%m-%d") # Q2
  SS = as.Date("2012-4-1",  format = "%Y-%m-%d") # Q3
  FE = as.Date("2012-7-1",  format = "%Y-%m-%d") # Q4
  d = as.Date(strftime(DATES, format="2012-%m-%d"))
  ifelse (d >= WS | d < SE, "Q4",
          ifelse (d >= SE & d < SS, "Q1",
                  ifelse (d >= SS & d < FE, "Q2", "Q3")))
}
jump.season = getSeason(jump.date.time[ , 1]) 

### Detected Jump Location (By Weekly Index) ###
all.week.jump.loca2.data=c()
for (i in 2:86){
  char = as.character(paste("jump.loca.week",i, sep = "."))
  if (eval(parse(text = char))[1] == "NA") {week.jump.loca2.data.i = "NA"}
  if ((eval(parse(text = char))[1] != "NA")) {
    week.jump.loca2.data.i = eval(parse(text = char))
  }
  all.week.jump.loca2.data = c(all.week.jump.loca2.data, 
                               week.jump.loca2.data.i)
}
all.week.jump.loca2.data = c(jump.loca.week.1, 
                             as.numeric(all.week.jump.loca2.data[
                               all.week.jump.loca2.data != "NA"]))

jump.detail.powerful = data.frame(jump.date = jump.date.time[,1], 
                                  jump.time = jump.date.time[,2], 
                                  jump.season = jump.season, 
                                  jump.size = exp(n.d$nitrate[
                                    all.week.jump.loca.data - 1]) * 
                                    (exp(all.week.jump.size.data) - 1))

jump.data = jump.detail.powerful
all.jump = jump.data$jump.size
posi.jump = jump.data$jump.size[which(jump.data$jump.size > 0)]
nega.jump = jump.data$jump.size[which(jump.data$jump.size < 0)]
nega.jump.abs = abs(nega.jump)

Q2.jump = jump.data$jump.size[which(jump.data$jump.season == "Q2")]
NonQ2.jump = jump.data$jump.size[-which(jump.data$jump.season == "Q2")]

### Runs test on Runs above or below 0 for All the Jumps ###
runs(all.jump)

### Two-sample Kolmogorov-Smirnov Test ### 
### for Positive Jump Sizes and Negative Jump Sizes ###
ks.test(posi.jump, nega.jump.abs)

### Two-sample Kolmogorov-Smirnov Test ### 
### for Quarter-2 Jump Sizes and Non-quarter-2 Jump Sizes ###
ks.test(abs(Q2.jump), abs(NonQ2.jump))


#### Figure 6 Producing ####
png("ACF.png", units = "in", height = 3, width = 10, res = 300)
acf(all.jump, lag.max = 20, main = "")
dev.off()


#### Figure 7 Producing ####
mybreaks1 <- seq(from = -9, to = 9, by = 0.3)
p11 = hist(posi.jump, breaks = mybreaks1, plot = FALSE)
p21 = hist(nega.jump, breaks = mybreaks1, plot = FALSE)

mybreaks <- seq(from = -9, to = 9, by = 0.30)
p1 = hist(c(abs(Q2.jump)), breaks = mybreaks, plot = FALSE)
p2 = hist(-c(abs(NonQ2.jump)), breaks = mybreaks, plot = FALSE)

png("Jump_Estimates_Comparison.png", 
    units = "in", height = 5.1, width = 9.7, res = 300)
par(mfrow = c(1, 2))
plot(p21, col = "white", border = "grey50", 
     xlim = c(-10, 10), ylim = c(0,100), 
     xlab = "Jump Estimate", ylab = "Frequency", 
     main = "(a) Jump Estimates Distribution")
plot(p11, col = "grey", xlim=c(-10, 10), border = "grey50", add = T)
plot(p2, col = "white", border = "grey50", 
     xlim=c(-10,10), ylim = c(0, 100), 
     xlab = "Jump Size", ylab = "Frequencey", xaxt = "n", 
     main = "(b) Quarter 2 (Grey) vs Other Quarters (White)")
plot(p1, xlim = c(-10,10), col = "grey", border = "grey50", 
     add = T)
axis(side = 1, at = c(-10, -5, 0, 5, 10), labels = c(10, 5, 0, 5, 10))
dev.off()


#### Table 1 Producing ####
### Quarter-2 ###
summary(abs(Q2.jump))
descdist(abs(Q2.jump), discrete = FALSE, graph = FALSE)
var(abs(Q2.jump))

### Non-Quarter-2 ###
summary(abs(NonQ2.jump))
descdist(abs(NonQ2.jump), discrete = FALSE, graph = FALSE)
var(abs(NonQ2.jump))


#### Table 2 Producing ####
abs.Q2.jump = abs(Q2.jump)
fit.gamma = fitdist(abs.Q2.jump, "gamma")
fit.ln = fitdist(abs.Q2.jump, "lnorm")
fit.weibull = fitdist(abs.Q2.jump, "weibull")
fit.burr = fitdist(abs.Q2.jump, "burr", 
                   start = list(shape1 = 1, shape2 = 1, rate = 2))
gofstat(list(fit.gamma, fit.ln, fit.weibull, fit.burr), 
        fitnames = c("gamma", "lnorm", "weibull", "burr"))
ks.test(abs.Q2.jump, "pgamma", shape = fit.gamma$estimate[1], 
        rate = fit.gamma$estimate[2])
ks.test(abs.Q2.jump, "plnorm", fit.ln$estimate[1], fit.ln$estimate[2])
ks.test(abs.Q2.jump, "pweibull", shape = fit.weibull$estimate[1], 
        scale = fit.weibull$estimate[2])
ks.test(abs.Q2.jump, "pburr", shape1 = fit.burr$estimate[1], 
        shape2 = fit.burr$estimate[2], rate = fit.burr$estimate[3])


#### Figure 8 Producing ####
Q2.fit = fit.ln
png("Quarter-2_Diagnostic.png", 
    units = "in", height = 5.1, width = 9.7, res = 300)
par(mfrow = c(1, 2))
denscomp(Q2.fit, addlegend = FALSE, fitcol = "black", fitlty = 2,
         fitpch = 6, lwd=1, border = "grey50",
         main = "(a) Fitted Density vs Empirical Histogram", 
         probability = TRUE, breaks = 75, xlab = "Jump Size",
         xlim = c(0, 9))
ppcomp(Q2.fit, fitcol = "grey50", addlegend = FALSE, main = "(b) P-P Plot",
       line01col = "black", line01lty = 2, 
       xlab = "Fitted Probability", ylab = "Empirical Probability")
dev.off()


#### Table 3 Producing ####
abs.NonQ2.jump = abs(NonQ2.jump)
fit.gamma = fitdist(abs.NonQ2.jump, "gamma")
fit.ln = fitdist(abs.NonQ2.jump, "lnorm")
fit.weibull = fitdist(abs.NonQ2.jump, "weibull")
fit.burr = fitdist(abs.NonQ2.jump, "burr", 
                   start = list(shape1 = 1, shape2 = 2, rate = 5))
gofstat(list(fit.gamma, fit.ln, fit.weibull, fit.burr), 
        fitnames = c("gamma", "lnorm", "weibull", "burr"))
ks.test(abs.NonQ2.jump, "pgamma", shape = fit.gamma$estimate[1], 
        rate = fit.gamma$estimate[2])
ks.test(abs.NonQ2.jump, "plnorm", fit.ln$estimate[1], fit.ln$estimate[2])
ks.test(abs.NonQ2.jump, "pweibull", shape = fit.weibull$estimate[1], 
        scale = fit.weibull$estimate[2])
ks.test(abs.NonQ2.jump, "pburr", shape1 = fit.burr$estimate[1], 
        shape2 = fit.burr$estimate[2], rate = fit.burr$estimate[3])


#### Figure 9 Producing ####
NonQ2.fit = fit.burr
png("Non-Quarter-2_Diagnostic.png", 
    units = "in", height = 5.1, width = 9.7, res = 300)

par(mfrow = c(1, 2))
denscomp(NonQ2.fit, addlegend = FALSE, fitcol = "black", fitlty = 2,
         fitpch = 6, lwd=1, border = "grey50",
         main = "(a) Fitted Density vs Empirical Histogram", 
         probability = TRUE, breaks = 75, xlab = "Jump Size",
         xlim = c(0, 9))
ppcomp(NonQ2.fit, fitcol = "grey50", addlegend = FALSE, 
       main = "(b) P-P Plot", 
       line01col = "black", line01lty = 2, 
       xlab = "Fitted Probability", ylab = "Empirical Probability")

dev.off()


#### Table 4 Producing ####
qlnorm(c(0.90, 0.95, 0.99), Q2.fit$estimate[1], Q2.fit$estimate[2])

qburr(c(0.90, 0.95, 0.99), NonQ2.fit$estimate[1], NonQ2.fit$estimate[2], 
      NonQ2.fit$estimate[3])
