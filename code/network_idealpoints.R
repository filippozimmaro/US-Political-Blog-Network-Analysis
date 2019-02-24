library(compiler)
enableJIT(3)
library(maxLik)
library(ggplot2)

setwd("/Users/shichaoma/Google Drive/Job Market/Non-academic/InsightInterview/")

rm(list=ls())
load("polblogsSmall.Rdata")

LogitNetworkIdealPointS.Two = function(Y, Tol = 1e-4){
  # Two step MLE
  # Y: sociomatrix
  llik.cost = function(b, Y, theta0, n){
    # n = nrow(Y)
    # cost = b
    costM = matrix(b,n,n)
    # c(0,theta0) # ideal points
    thetaM = matrix(c(0,theta0),n,n)
    utility = -(thetaM - t(thetaM))^2 - costM
    ll = - Y * log(1+exp(-utility)) - (1-Y) * log(1+exp(utility))
    diag(ll) = NA
    return(na.omit(c(t(ll))))
  }
  grad.cost = function(b,Y, theta0, n){
    costM = matrix(b,n,n)
    # c(0,theta0) # ideal points
    thetaM = matrix(c(0,theta0),n,n)
    utility = -(thetaM - t(thetaM))^2 - costM
    dC = -Y * exp(-utility) / (1+exp(-utility)) + (1-Y) * exp(utility) / (1+exp(utility))
    diag(dC) = 0
    return(apply(dC,1,"sum"))
  }
  hess.cost = function(b,Y, theta0, n){
    costM = matrix(b,n,n)
    # c(0,theta0) # ideal points
    thetaM = matrix(c(0,theta0),n,n)
    utility = -(thetaM - t(thetaM))^2 - costM
    dC2 = -Y * exp(-utility) / (1+exp(-utility))^2 - (1-Y) * exp(utility) / (1+exp(utility))^2
    diag(dC2) = 0
    return(dC2)
    # return(apply(dC2,1,"sum"))
  }
  llik.theta = function(b,Y,cost,n){
    # n = nrow(Y)
    # c(0,b) # ideal points
    thetaM = matrix(c(0,b),n,n)
    costM = matrix(cost,n,n)    
    utility = -(thetaM - t(thetaM))^2 - costM
    ll = - Y * log(1+exp(-utility)) - (1-Y) * log(1+exp(utility))
    diag(ll) = NA
    return(na.omit(c(t(ll))))
  }
  grad.theta = function(b,Y,cost,n){
    # n = nrow(Y)
    # c(0,b) # ideal points
    thetaM = matrix(c(0,b),n,n)
    costM = matrix(cost,n,n)    
    utility = -(thetaM - t(thetaM))^2 - costM
    dtheta = -Y * 2 * exp(-utility) * (thetaM - t(thetaM)) / (1+exp(-utility)) + 
      (1-Y) * 2 * exp(utility) * (thetaM - t(thetaM)) / (1+exp(utility))
    diag(dtheta) = 0
    return((apply(dtheta,1,"sum") + apply(-dtheta,2,"sum"))[-1])
  }
  hess.theta = function(b,Y,cost,n){
    # n = nrow(Y)
    # c(0,b) # ideal points
    thetaM = matrix(c(0,b),n,n)
    costM = matrix(cost,n,n)    
    utility = -(thetaM - t(thetaM))^2 - costM
    dtheta2 = -Y * (4 * exp(-utility) * (thetaM - t(thetaM))^2 + 2 * exp(-utility) * (1+exp(-utility))) / (1+exp(-utility))^2 + 
      (1-Y) * (-4 * exp(utility) * (thetaM - t(thetaM))^2 + 2 * exp(utility) * (1+exp(utility))) / (1+exp(utility))^2
    diag(dtheta2) = 0
    return(dtheta2[-1,-1])
    #return((apply(dtheta2,1,"sum") + apply(dtheta2,2,"sum"))[-1])
  }
  n = nrow(Y)
  StartValue.cost = runif(n)
  StartValue.theta = runif(n-1)
  costGood = F
  thetaGood = F
  thetall = 100 # initial likelihood
  costll = 100 # initial likelihood
  loop = 1
  while( (costGood == F) || (thetaGood == F)){
    # cost loop
    costMax = maxLik(llik.cost, grad = grad.cost, 
                     #hess = hess.cost,  #Hessian is not used by BFGS
                     start = StartValue.cost, method="BFGS", 
                     iterlim = 500, tol = Tol, finalHessian = F,
                     Y = Y, theta0 = StartValue.theta, n = n)
    if(abs(costMax$maximum - costll) < Tol){
      costGood = T
    }else{
      costGood = F
    }
    costll = costMax$maximum
    StartValue.cost = costMax$estimate
    cat("Cost done at loop",loop, "!\n")
    
    # theta loop
    thetaMax = maxLik(llik.theta, grad = grad.theta, 
                      # hess = hess.theta,  #Hessian is not used by BFGS
                      start=StartValue.theta, method="BFGS",
                      iterlim = 500, tol = Tol, finalHessian = F,
                      Y = Y, cost = StartValue.cost, n = n)
    if(abs(thetaMax$maximum - thetall) < Tol){
      thetaGood = T
    }else{
      thetaGood = F
    }
    thetall = thetaMax$maximum
    StartValue.theta = thetaMax$estimate
    cat("Theta done at loop",loop, "!\n")
    
    cat("loop =", loop, "has done!\n\n")
    loop = loop + 1
  }
  ## One step at the end to get Hessian
  print("Finalizing...\n")
  costMax = maxLik(llik.cost, grad = grad.cost, 
                   start = StartValue.cost, method="BFGS", 
                   iterlim = 500,
                   Y = Y, theta0 = StartValue.theta, n = n)
  thetaMax = maxLik(llik.theta, grad = grad.theta, 
                    start=StartValue.theta, method="BFGS",
                    iterlim = 500,
                    Y = Y, cost = StartValue.cost, n = n)
  return(list(thetaMax,costMax))
}

# do the estimation
blog.result = LogitNetworkIdealPointS.Two(Y = matrix.polblogs.small)
save.image(file="bolgResult.RData")

load("bolgResult.RData")

print("Sumamry of theta:\n")
summary(blog.result[[1]]$estimate)
print("Sumamry of cost:\n")
summary(blog.result[[2]]$estimate)

## calibrate the middle point
Calibrate = function(input,label){
  m.R = median(input[as.logical(label)])
  m.D = median(input[!as.logical(label)])
  midPoint = (m.R + m.D)/2
  output = input - midPoint
  if(m.R < midPoint && m.D > midPoint){
    output = -output
  }
  return(output)
}

idealPoints = Calibrate(c(0,blog.result[[1]]$estimate), LeftRight.small)

# evaluate the model using self-identification
Score = function(input, label){
  correct = (input>=0) == as.logical(label)
  return(sum(correct)/length(correct))
}
Score(idealPoints, LeftRight.small)

# plot distribution of ideal points
histPlot = ggplot(data.frame(idealPoints, LeftRight.small)) + 
  geom_histogram(aes(x=idealPoints, fill = factor(1 - LeftRight.small), 
                     y = ..density..), 
                 binwidth=.1,alpha=.5, position="identity") +
  geom_density(aes(x=idealPoints)) +
  xlab("Ideal Point") +
  ylab("Density") + 
  scale_fill_discrete(name = "Political Leaning",
                      breaks=c("0", "1"),
                      labels=c("Conversitive", "Liberal"))
histPlot
ggsave(filename = "hist_density.pdf", plot = histPlot)

# plot the standard deviation of ideal points
standard.dev = sqrt(diag(vcov(blog.result[[1]])))
sdPlot = ggplot(data.frame(idealPoints[-1], standard.dev, LeftRight.small[-1])) + 
  geom_point(aes(x=idealPoints..1., y = standard.dev,
                 colour = factor(1 - LeftRight.small..1.)
  )) +
  xlab("Ideal Point") +
  ylab("Standard Error") + 
  scale_colour_discrete(name = "Political Leaning",
                        breaks=c("0", "1"),
                        labels=c("Conversitive", "Liberal"))
sdPlot
ggsave(filename = "sd.pdf", plot = sdPlot)

outgoing = apply(matrix.polblogs.small, 1, "sum")
outPlot = ggplot(data.frame(idealPoints, outgoing, LeftRight.small)) + 
  geom_point(aes(x=idealPoints, y = outgoing,
                 colour = factor(1 - LeftRight.small))) +
  xlab("Ideal Point") +
  ylab("# Outgoing Relations") + 
  scale_colour_discrete(name = "Political Leaning",
                        breaks=c("0", "1"),
                        labels=c("Conversitive", "Liberal"))
outPlot
ggsave(filename = "outgoing.pdf", plot = outPlot)

incoming = apply(matrix.polblogs.small, 2, "sum")
inPlot = ggplot(data.frame(idealPoints, incoming, LeftRight.small)) + 
  geom_point(aes(x=idealPoints, y = incoming,
                 colour = factor(1 - LeftRight.small))) +
  xlab("Ideal Point") +
  ylab("# Incoming Relations") + 
  scale_colour_discrete(name = "Political Leaning",
                        breaks=c("0", "1"),
                        labels=c("Conversitive", "Liberal"))
inPlot
ggsave(filename = "incoming.pdf", plot = inPlot)


AverageIdeology = function(idealPoints, Y){
  averageOutgoingIdeology = rep(NA,length(idealPoints))
  averageIncomingIdeology = rep(NA,length(idealPoints))
  for(i in 1:length(idealPoints)){
    averageOutgoingIdeology[i] = mean(ifelse(Y[i, ] == 1, idealPoints, NA),
                                      na.rm = T)
    averageIncomingIdeology[i] = mean(ifelse(Y[, i] == 1, idealPoints, NA), 
                                      na.rm = T)
  }
  return(cbind(averageOutgoingIdeology, averageIncomingIdeology))
}
avgIdeology = AverageIdeology(idealPoints, matrix.polblogs.small)
avgOut = avgIdeology[, 1]
avgIn = avgIdeology[, 2]

avgOutPlot = ggplot(data.frame(idealPoints, avgOut, LeftRight.small)) + 
  geom_point(aes(x=idealPoints, y = avgOut,
                 colour = factor(1 - LeftRight.small))) +
  geom_abline(intercept=0, slope=1) +
  xlab("Ideal Point") +
  ylab("Average Ideal Points of Outgoing Relations") + 
  scale_colour_discrete(name = "Political Leaning",
                        breaks=c("0", "1"),
                        labels=c("Conversitive", "Liberal"))
avgOutPlot
ggsave(filename = "avg_out.pdf", plot = avgOutPlot)

avgInPlot = ggplot(data.frame(idealPoints, avgIn, LeftRight.small)) + 
  geom_point(aes(x=idealPoints, y = avgIn,
                 colour = factor(1 - LeftRight.small))) +
  geom_abline(intercept=0, slope=1) +
  xlab("Ideal Point") +
  ylab("Average Ideal Points of Incoming Relations") + 
  scale_colour_discrete(name = "Political Leaning",
                        breaks=c("0", "1"),
                        labels=c("Conversitive", "Liberal"))
avgInPlot
ggsave(filename = "avg_in.pdf", plot = avgInPlot)

cost = blog.result[[2]]$estimate
costPlot = ggplot(data.frame(idealPoints, cost, LeftRight.small)) + 
  geom_point(aes(x=idealPoints, y = cost,
                 colour = factor(1 - LeftRight.small))) +
  xlab("Ideal Point") +
  ylab("Blog Fixed Effect (Net Cost)") + 
  scale_colour_discrete(name = "Political Leaning",
                        breaks=c("0", "1"),
                        labels=c("Conversitive", "Liberal"))
costPlot
ggsave(filename = "cost.pdf", plot = costPlot)

save.image(file="bolgResult.RData")
