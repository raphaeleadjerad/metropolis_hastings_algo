# Metropolis Hastings Algorithm ----------------------------
# Adjerad R. 


# packages ---------------------------------
library(ggplot2)
library(boot)
library(invgamma)
library(cumstats)
library(gridExtra)
library(bayesmeta)

# Import dataset ---------------------------
dataset <- cars
speed <- dataset$speed
speed.2 <- speed^2
dist <- dataset$dist
dataset.2 <- cbind(dataset,speed.2)

# linear regression ----------------
ggplot(data = dataset.2)+
  geom_point(mapping = aes(x = speed, y = dist), alpha = 0.6)+
  geom_smooth(aes(x = speed, y = dist))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Distance to stop increases with speed",
          caption = "Computational Statistics", 
          x = "Speed",
          y = 'Distance to stop')
ggplot(data = dataset.2)+
  geom_point(mapping = aes(x = speed.2, y = dist), alpha = 0.6)+
  geom_smooth(aes(x = speed.2, y = dist))+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Distance to stop increases with speed^2",
          caption = "Computational Statistics", 
          x = "Speed squared",
          y = 'Distance to stop')

model <- dist ~ speed + speed.2
fit.1 <- lm(model)
summary(fit.1)

coeff = summary(fit.1)[["coefficients"]]                                 # return coefficient vector
coeff.fit.1 <- c("a.hat" = coeff[1,1], "b.hat" = coeff[2,1], "c.hat" = coeff[3,1],
                 "sigma.2" = (summary(fit.1)$sigma)^2,"sd.a" = coeff[1,2], "sd.b" = coeff[2,2], "sd.c" = coeff[3,2])
coeff.fit.1 <- as.data.frame(t(coeff.fit.1))
coeff.fit.1
# Histogram of posterior distribution -------------------------------

# In order to get empirical histograms of our coefficients, we bootstrap them.
# Parameters
R <- 10000
set.seed(134)
# Bootstrap
boot.coef <- function(data, indices){
  data <- data[indices,]                                                 
  # select obs. in bootstrap sample
  mod <- lm(dist ~ speed + speed.2, data= data)
  coeff = summary(mod)[["coefficients"]]                                 
  
  return(c("a.hat" = coeff[1,1], "b.hat" = coeff[2,1], "c.hat" = coeff[3,1],
           "sigma.2" = (summary(mod)$sigma)^2,"sd.a" = coeff[1,2], "sd.b" = 
             coeff[2,2], "sd.c" = coeff[3,2]))
}


boot.fit.1 <- boot(data = dataset.2, statistic = boot.coef, R = R)
rep.boot <- as.data.frame(boot.fit.1['t'])
colnames(rep.boot) <- c("a.hat.b" , "b.hat.b" , "c.hat.b" ,
                        "sigma.2.b" ,"sd.a.b" , "sd.b.b" , "sd.c.b")


# plot ----------------------------------------------------------------
# a.hat.bootstrapped
ggplot(data = rep.boot,mapping = aes(x = a.hat.b))+
  geom_histogram(bins = 15, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = a.hat.b, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',size=1) + 
  stat_function(fun = dnorm, args = list(mean = coeff.fit.1$a.hat, 
                                         sd =coeff.fit.1$sd.a), 
                aes(colour = 'Theoretical normal density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Empirical distribution of coef a",
          caption = "Computational Statistics", 
          x = "a",
          y = 'Empirical density')

# b.hat.bootstrapped
ggplot(data = rep.boot,mapping = aes(x = b.hat.b))+
  geom_histogram(bins = 15, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = b.hat.b, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',size=1) + 
  stat_function(fun = dnorm, args = list(mean = coeff.fit.1$b.hat, 
                                         sd =coeff.fit.1$sd.b), 
                aes(colour = 'Theoretical normal density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Empirical distribution of coef b",
          caption = "Computational Statistics", 
          x = "b",
          y = 'Empirical density')

# c.hat.bootstrapped
ggplot(data = rep.boot,mapping = aes(x = c.hat.b))+
  geom_histogram(bins = 15, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = c.hat.b, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',size=1) + 
  stat_function(fun = dnorm, args = list(mean = coeff.fit.1$c.hat, 
                                         sd =coeff.fit.1$sd.c), 
                aes(colour = 'Theoretical normal density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Empirical distribution of coef c",
          caption = "Computational Statistics", 
          x = "c",
          y = 'Empirical density')

# sigma.hat.bootstrapped
args.1 <-  1+ length(dist)/2
args.2 <-  1 + 1/2*sum(fit.1$residuals^2)

ggplot(data = rep.boot,mapping = aes(x = sigma.2.b))+
  geom_histogram(bins = 15, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = sigma.2.b, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',size=1) + 
  stat_function(fun = dinvgamma, args =list(args.1,args.2), 
                aes(colour = 'Theoretical inverse gamma density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"))+ labs(
          title = "Empirical distribution of coef sigma2",
          caption = "Computational Statistics", 
          x = "sigma2",
          y = 'Empirical density')

# Metropolis Hastings --------------------------

# Parameters
param.init <- as.matrix(coeff.fit.1[1,])
iteration  <-  100000
start.value <-  param.init[1:4]
ssr <- sum(fit.1$residuals^2)
set.seed(139)

# Posterior distribution 
l.likelihood <- function(theta){
  ll <- (-length(dist)/2)*log(theta[4]) - 
    (2*theta[4])^(-1)*sum((dist - theta[1]-
                             theta[2]*speed - 
                             theta[3]*speed.2)^2)
  return(ll)
}

# Proposal density : We pick the product of the marginals 
l.P.dens <-  function(theta){
  
  l.a.dens <-  dnorm(theta[1],param.init[1],param.init[5],log=T)
  l.b.dens <- dnorm(theta[2],param.init[2],param.init[6],log=T)
  l.c.dens <-  dnorm(theta[3],param.init[3],param.init[7],log=T)
  l.sigma.2.dens <-  dinvgamma(theta[4],1+length(dist)/2,1+ssr/2,log=T)
  
  l.P.density <- l.a.dens + l.b.dens + l.c.dens + l.sigma.2.dens
  return(l.P.density)
}


# Proposal P(x,) : Note that it does not depend on where the chain is at 
P <-  function(){
  a.p <- rnorm(1,param.init[1],param.init[5])
  b.p <- rnorm(1,param.init[2],param.init[6])
  c.p <- rnorm(1,param.init[3],param.init[7])
  sigma.2.p <- rinvgamma(1,1+length(dist)/2,1+ssr/2)
  return(c(a.p, b.p, c.p, sigma.2.p))
} 

# Metropolis algorithm
algo.metropolis <- function(iteration, start.value = start.value, mu, proposal, P){
  
  chain <- matrix(NA,
                  ncol = length(start.value), nrow = iteration, byrow = T)
  chain[1,] <- start.value
  
  for(i in 1:(iteration-1)){
    x <- as.vector(chain[i,])
    y <- P()
    ratio <- exp(mu(y) + proposal(x) - mu(x)- proposal(y))
    rho <- min(1,ratio)
    accept <- runif(1) < rho
    if(accept){
      chain[i+1,] <- y
    }else{
      chain[i+1,]<- chain[i,]
    }
  }
  return(chain)
}

algo.metro.1 <-  algo.metropolis(iteration, start.value, l.likelihood, l.P.dens, P)

# Burnin
burnIn <- 10000
algo.metro.1 <- data.frame(algo.metro.1)
chain.metro <- cbind(N.iter = 1:iteration,algo.metro.1)
chain.metro.burnin <- chain.metro[-(1:burnIn),]
# We burnin the first 10000 observations, too dependent on the beginning
colnames(chain.metro.burnin) <- c("N.iter","a","b","c","sigma.2")

# plot ----------------------------------------------------
ggplot(data = chain.metro.burnin, mapping = aes(x = a))+
  geom_histogram(bins = 25, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = a, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',adjust = 3,size=1) + 
  stat_function(fun = dnorm,args = list(param.init[1],param.init[5]), 
                aes(colour = 'Theoretical normal density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Histogram of a",
          caption = "Computational Statistics", 
          x = "a")

ggplot(data = chain.metro.burnin, mapping = aes(x = b))+
  geom_histogram(bins = 25, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = b, y = ..density.., colour = 'Empirical density'),
            stat = 'density',adjust = 3,size=1) + 
  stat_function(fun = dnorm,args = list(param.init[2],param.init[6]), 
                aes(colour = 'Theoretical normal density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Histogram of b",
          caption = "Computational Statistics", 
          x = "b")

ggplot(data = chain.metro.burnin, mapping = aes(x = c))+
  geom_histogram(bins = 25, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = c, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',adjust = 3,size=1) + 
  stat_function(fun = dnorm,args = list(param.init[3],param.init[7]), 
                aes(colour = 'Theoretical normal density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Histogram of c",
          caption = "Computational Statistics", 
          x = "c")

ggplot(data = chain.metro.burnin, mapping = aes(x = sigma.2))+
  geom_histogram(bins = 25, fill = "white", colour = "black", aes(y =..density..))+
  geom_line(aes(x = sigma.2, y = ..density.., colour = 'Empirical density'), 
            stat = 'density',adjust = 3,size=1) + 
  stat_function(fun = dinvgamma,args = list(1+length(dist)/2,1+ssr/2), 
                aes(colour = 'Theoretical invgamma density'), size = 1) +
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Histogram of sigma.2",
          caption = "Computational Statistics", 
          x = "sigma.2")

ggplot(data = chain.metro.burnin, mapping = aes(x = N.iter,y=a))+
  geom_line(colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Chain of a",
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Value of a")

ggplot(data = chain.metro.burnin, mapping = aes(x = N.iter,y=b))+
  geom_line(colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Chain of b",
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Value of b")

ggplot(data = chain.metro.burnin, mapping = aes(x = N.iter,y=c))+
  geom_line(colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Chain of c",
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Value of c")

ggplot(data = chain.metro.burnin, mapping = aes(x = N.iter,y=sigma.2))+
  geom_line(colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = "Chain of sigma.2",
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Value of sigma.2")


# Monitoring convergence -------------------------------


m.chain.metro <- chain.metro.burnin
for (i in 1:4){
  m.chain.metro[,i+1] <-  cummean(m.chain.metro[,i+1])
}

# Cumulative mean a
g1 <- ggplot(data = m.chain.metro)+
  geom_line(mapping = aes(x = N.iter,y = a),colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = paste0("Cumulative mean a"),
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Cum mean a")

# Cumulative mean b
g2 <- ggplot(data = m.chain.metro)+
  geom_line(mapping = aes(x = N.iter,y = b),colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = paste0("Cumulative mean b"),
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Cum mean b")

# Cumulative mean c
g3 <- ggplot(data = m.chain.metro)+
  geom_line(mapping = aes(x = N.iter,y = a),colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = paste0("Cumulative mean c"),
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Cum mean c")

# Cumulative mean sigma
g4 <- ggplot(data = m.chain.metro)+
  geom_line(mapping = aes(x = N.iter,y = a),colour = "black",alpha = 0.6)+
  theme(text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(face ="italic"),
        legend.position = "bottom")+ labs(
          title = paste0("Cumulative mean sigma2"),
          caption = "Computational Statistics", 
          x = "Iterations",
          y ="Cum mean sigma2")

grid.arrange(g1,g2,g3,g4, nrow = 2, ncol = 2)

# Metropolis Hastings with Student distribution -------------------
# Posterior distribution 
nu <- 4
l.likelihood.t <- function(theta){
  ll <- (-length(dist)/2)*log(theta[4]) -((nu+1)/2)*
    sum(log(1+((dist - theta[1]-theta[2]*speed -theta[3]*speed.2)^2)/(nu*theta[4])))
  return(ll)
}
l.likelihood.t(param.init)

# Proposal density : We pick the product of the marginal
l.P.t.dens <-  function(theta){
  # Marginals for each parameter
  l.t <- function(dist,nu,param,sd.2) -1/2*log(sd.2) - 
    (nu+1)/2*log(1+(dist-param)^2/(nu*sd.2))
  l.a.dens <-  l.t(theta[1],nu,param.init[1],param.init[5])
  l.b.dens <- l.t(theta[2],nu,param.init[2],param.init[6])
  l.c.dens <-  l.t(theta[3],nu,param.init[3],param.init[7])
  l.sigma.2.dens <-  dhalft(theta[4],scale = (ssr/(length(dist)-3)), nu,log=T)
  
  l.P.density <- l.a.dens + l.b.dens + l.c.dens + l.sigma.2.dens
  return(l.P.density)
}
l.P.t.dens(param.init)

# Proposal 
P.t <-  function(){
  a.p <- sqrt(param.init[5])*rt(1,df = nu, ncp = param.init[1])
  b.p <- sqrt(param.init[6])*rt(1,df = nu, ncp = param.init[2])
  c.p <- sqrt(param.init[7])*rt(1,df = nu, ncp = param.init[3])
  sigma.2.p <- rhalft(1,scale=(ssr/(length(dist)-3)), nu)
  return(c(a.p, b.p, c.p, sigma.2.p))
} 
P.t()

# Output of algorithm ------------------------------------

# takes a minute !
algo.metro.1 <-  algo.metropolis(iteration, start.value, l.likelihood.t, l.P.t.dens, P.t)

# Burnin
burnIn <- 10000
algo.metro.1 <- data.frame(algo.metro.1)
chain.metro <- cbind(N.iter = 1:iteration,algo.metro.1)
chain.metro.burnin <- chain.metro[-(1:burnIn),]
# We burnin the first 10000 observations, too dependent on the beginning
colnames(chain.metro.burnin) <- c("N.iter","a","b","c","sigma.2")


