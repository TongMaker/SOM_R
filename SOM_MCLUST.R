library(MASS)

source("/Users/tongzhou/Library/Mobile Documents/com~apple~CloudDocs/Desktop/PHD/RCODE/Aux_Program/Mixture_Parameter_Generator.R", echo=TRUE)
#function to split data
sample_split<-function(samples, num_split){
  sample_size<-nrow(samples)
  size_miniSample<-sample_size/num_split
  if(size_miniSample-floor(size_miniSample)!=0){
    stop("num_split is not a multiple of samples")
  }
  samples_split<-list()
  cont<-1:size_miniSample
  
  for (L in 1:num_split){
    samples_split[[L]]<-samples[cont,]
    cont<-cont+size_miniSample
  }
  samples_split
}



#function of SOM Algorithm
SOM_BASIC<-function(samples, gridsom, cluster, eta, sigma0, weight=NULL){
  gridnum<-ncol(gridsom)
  dimension<-ncol(samples)
  if(is.null(weight)){weight<-matrix(runif(gridnum*dimension),ncol = dimension)}
  
  tmax<-nrow(samples)
  ri<-array(dim = c(tmax,nrow(gridsom)))
  i<-c()
  for(t in 1:tmax){
    dist_WeightAndData<-apply(weight, 1,function(x) sqrt(sum((samples[t,]-x)^2)))
    i_star<-which.min(dist_WeightAndData) #the selected weight is at position i_star, 
    i[t]<-i_star
    sigma<-sigma0*exp(-2*sigma0*t/tmax)
    ri_star<-gridsom[, i_star]
    lambda<-apply(gridsom, 2, function(x) exp(-sum((x-ri_star)^2)/(2*sigma^2)))
    Dweight<-matrix(rep(eta*lambda, dimension),cluster)*t(apply(weight, 1,function(x) samples[t,]-x))
    weight<-weight+Dweight
    ri[t,]<-ri_star
  }
  list(weight=weight, classification=cbind(samples, ri,i))
}


means<-Gpoints(num_point=6, dimension=3, dmin=100, pinit=5)
covs<-create_param(D=3,K=6, type="cov", covdiag = 100)


means<-matrix(c(10,10,10,30,20,20,10,10,40,20,20,10,30,40,40,30,10,30), ncol = 6)
covs<-array(rep(20*diag(3),6), dim = c(3,3,6))
weights<-c(0.18,0.12,0.25,0.1,0.2,0.15)
sizek<-1000000
componente<-sample(1:6, size = sizek, replace = T, prob = weights)

samples<-do.call(rbind,
                 lapply(1:6,function(k){
                   mues_k<-mvrnorm(n = sum(componente==k), means[,k], covs[,,k])
                   cbind(mues_k,k)
                 }
                 )
)[sample(sizek),] 

samples_k<-sample_split(samples,1000)
SOM_res<-lapply(samples_k, \(x) SOM_BASIC(x[,c(1,2,3)],gridsom, cluster = 6, eta = 0.5, sigma0 = 0.5))
centroid<-do.call(rbind,lapply(SOM_res, \(x) x[[1]]))

library(mclust)

mclust_est<-mclust::Mclust(centroid, G=6)
mclust_est$parameters$mean



chatgpt::ask_chatgpt()

library(chatgpt)

gpt_get_completions(
  prompt,
  openai_api_key = Sys.getenv("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"),
  messages = NULL
)
Sys.setenv(OPENAI_API_KEY = "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")