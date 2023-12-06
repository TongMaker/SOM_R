library(rstan)

means<-matrix(c(10,50,10,60,20,20,80,10,40,20,20,10,30,40,40,30,10,30,50,50,50,40,60,10,80,20,60), ncol = 9)
covs<-array(rep(20*diag(3),9), dim = c(3,3,9))
weights<-c(0.08,0.1,0.12,0.1,0.15,0.1,0.2,0.1,0.05)

sizek<-1000
componente<-sample(1:9, size = sizek, replace = T, prob = weights)

samples<-do.call(rbind,
                 lapply(1:9,function(k){
                   mues_k<-MASS::mvrnorm(n = sum(componente==k), means[,k], covs[,,k])
                   cbind(mues_k,k)
                 }
                 )
)[sample(sizek),][,1:3] 

fitstan = rstan::stan(file = "SoftKmean.stan",
                      data = list(
                        N=1000,
                        D=3,
                        K=9,
                        y=samples),
                      iter=1100, warmup=100, chains=1)
