source("/Users/tongzhou/Library/Mobile Documents/com~apple~CloudDocs/Desktop/PHD/RCODE/Aux_Program/Mixture_Parameter_Generator.R", echo=TRUE)
library(tidyverse)
library(kohonen)
library(MASS)
#
ProbByNode_AllSom_K<-list()
ProbByNode_Total_K<-list()
Prob_AllSom_K<-list()
Prob_Total_K<-list()
Welldone_K<-list()

for(K in 1:1000){
print(K)
tini<-Sys.time()  
sample_size<-1000000
#data dimension
dd=3
#number cluster
kk=6
#weights<-diff(c(0,sort(sample(seq(0,1,1e-5),kk-1)),1))  
weights<-rep(1/kk,kk)
sum(weights)

#means<-create_param(D=dd,K=kk, type="mean")
covs<-create_param(D=dd,K=kk, type="cov", covdiag = 100)

means<-Gpoints(num_point=6, dimension=3, dmin=100, pinit=123)

componente<-sample(1:kk, size = sample_size, replace = T, prob = weights)
#Main sample
samples<-do.call(rbind,
                 lapply(1:kk,function(k){
                   mues_k<-mvrnorm(n = sum(componente==k), means[,k], covs[,,k])
                   cbind(mues_k,k)
                 }
                 )
)[sample(sample_size),] 


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


#function to verify if two matrix are similar no matter the rows ordering
SOM_SUBSET_FIT2<-function(m,s,threshold_n){
  rp<-(runif(ncol(m))-0.5)*1e5
  matrp<-t(array(rp, dim = c(ncol(m), nrow(m))))
  
  distrp_m<-sqrt(rowSums((m-matrp)^2))
  distrp_s<-sqrt(rowSums((s-matrp)^2))
  
  m_ordered<-m[order(distrp_m),]
  s_ordered<-s[order(distrp_s),]
  
  dist_m_s<-sqrt(rowSums((m_ordered-s_ordered)^2))
  max(dist_m_s)<threshold_n
}


samples_split <- sample_split(samples, 1000)


df<-data.frame(samples_split[[1]] )
names(df)<-c("X","Y","Z","k")
plotly::plot_ly(df, x=~X, y=~Y, z=~Z,size=0.001, type = "scatter3d", mode = "markers" )


gridsom<- matrix(c(1,1,1,2,2,1,2,2,3,1,3,2), nrow = 2)
#gridsom<- matrix(c(1,1,1.5,2,2,1,2.5,2,3,1,3.5,2), nrow = 2)
plot(t(gridsom),asp=1)

#the relationship between grid and weight is the order of matrix
res1<-SOM_BASIC(samples_split[[1]][,c(1,2,3)],gridsom, cluster = 6, eta = 0.1, sigma0 = 0.5)

SOM_res<-lapply(samples_split, \(x) SOM_BASIC(x[,c(1,2,3)],gridsom, cluster = 6, eta = 0.5, sigma0 = 0.5))

WelldoneSOM<-lapply(SOM_res, \(x) SOM_SUBSET_FIT2(t(means),x$weight,20)) %>% unlist()
table(WelldoneSOM)



vector_join<-function(long_vector, join_vector){
short_table<-data.frame(a=1:length(join_vector),b=join_vector) 
Long_table<-data.frame(a=long_vector,ord=1:length(long_vector))
cc<-merge(short_table, Long_table, by="a")
dd<-cc[order(cc$ord),]$b
dd
}

#load("/Users/tongzhou/Library/Mobile Documents/com~apple~CloudDocs/Desktop/PHD/RCODE/Rdata/confusion_mat.RData")

#loop of confusion matrices
confusion_mat<-list()
for(i in 1:length(SOM_res)){
  res<-SOM_res[[i]]$classification
  ws<-SOM_res[[i]]$weight
  data<-samples_split[[i]]
  wm<-t(means)
  
  # 5,5,1,3,6,4 estimacion som 1 <=> mean 5, som 2 <=> mean 5, som 3 <=> mean 1, som 4 <=> mean 3, som 5 <=> mean 6, som 6 <=> mean 4
  wsorder<-apply(ws,1,\(x) which.min(sqrt(rowSums(t(x-t(wm))^2))))

  
  #match: replace the number of vector: value of position ->> position number
  dataposi<-data[,4]
  orderedposi<-vector_join(res[,6], wsorder)             
  no_exist_node<-setdiff(1:nrow(wm), unique(orderedposi)) #puede ser mas de un elemento
  
  #Add the non-existent estimated node in order to place the column correctly in the confusion matrix.
  confusion_mat[[i]]<-table(c(no_exist_node,dataposi),c(no_exist_node,orderedposi))
  #substract the 1 unit of the non-existent column(s)
  confusion_mat[[i]][no_exist_node,no_exist_node]<-confusion_mat[[i]][no_exist_node,no_exist_node]-diag(length(no_exist_node))
  
}



ProbByNode_AllSom<-do.call(rbind,lapply(1:length(confusion_mat), \(i) diag(confusion_mat[[i]])/rowSums(confusion_mat[[i]])))
ProbByNode_Total<-colSums(ProbByNode_AllSom)


Prob_AllSom<-unlist(lapply(1:length(confusion_mat), \(i) sum(diag(confusion_mat[[i]]))/sum(confusion_mat[[i]])))
Prob_Total<-sum(Prob_AllSom)

quantile(Prob_AllSom, 0.05)

#result store
ProbByNode_AllSom_K[[K]]<-ProbByNode_AllSom
ProbByNode_Total_K[[K]]<-ProbByNode_Total
Prob_AllSom_K[[K]]<-Prob_AllSom
Prob_Total_K[[K]]<-Prob_Total
Welldone_K[[K]]<-table(WelldoneSOM)
tfin<-Sys.time()

print(paste0("Time consumption of iteration ",K,": ",floor(tfin-tini)," seconds"))
}

load("/Users/tongzhou/Library/Mobile Documents/com~apple~CloudDocs/Desktop/PHD/RCODE/Rdata/Som_uncertainty1000_1000_1000.RData")


which(Prob_AllSom<0.74)
confusion_mat[[723]]
Welldone_K[[4]]
hist(ProbByNode_AllSom_K[[1]][,5], breaks = 200)
quantile(ProbByNode_AllSom_K[[1]][,3], 0.1)

hist(Prob_AllSom_K[[1000]], breaks = 100)
plot(density(Prob_AllSom_K[[1000]]))
plot(density(unlist(Prob_Total_K)))


hist(unlist(Prob_Total_K),breaks = 100)
plot(density(unlist(Prob_Total_K)))
WD<-do.call(rbind,Welldone_K)
hist(WD[,1] ,breaks = 50)


#====================================================================================================================================

ggplot(data.frame(Precision=Prob_AllSom_K[[145]]), aes(x=Precision)) + geom_density(adjust=2, fill="#5d8aa8") +guides(fill="none", colour="none")+theme_light()
summary(Prob_AllSom_K[[1000]])
summary(ProbByNode_AllSom_K[[146]])->s
ratios<-ProbByNode_AllSom_K[[146]]
mean_ratios<-rbind(reduce(ProbByNode_AllSom_K, `+`)/1000, t(rep(0.8,6)))
all_ratios<-do.call(rbind,ProbByNode_AllSom_K)
apply(mean_ratios, 2, function(x)quantile(x,0.005))
apply(mean_ratios, 2, function(x)quantile(x,0.01))


round(t(apply(all_ratios, 2, function(x)quantile(x,c(0.05,0.1,0.2,0.3,0.4,0.5,0.75,0.9,1)))),4)


#paste0("Cluster 6"," & ",paste(t(s)[6,] %>% str_remove_all("[:alpha:]+\\s*[:punct:]*\\s*[:punct:]*") %>% str_trim(), collapse = " & "), "\\\\")
#PrecisionByNode<-as.data.frame(ProbByNode_AllSom_K[[146]])%>%gather(key="cluster", value = "Precision") %>% mutate(cluster=paste0("Cluster ",cluster))
PrecisionByNode<-as.data.frame(mean_ratios)%>%gather(key="cluster", value = "Precision") %>% mutate(cluster=paste0("Cluster ",cluster))
PrecisionByNode<-as.data.frame(all_ratios)%>%gather(key="cluster", value = "Precision") %>% mutate(cluster=paste0("Cluster ",cluster))

options(scipen = 999)
ggplot(PrecisionByNode, aes(x=Precision)) + 
  geom_histogram(binwidth=.008, alpha=.5, position="identity")+
  guides(fill="none", colour="none")+
  facet_wrap(~cluster, ncol=3, scales = "free_y")+
  scale_y_continuous(labels = scales::label_number_si())+
  #scale_x_continuous(labels = scales::percent)+
  theme_light()

ggplot(PrecisionByNode, aes(x=Precision)) + 
  geom_density(fill="gray", adjust=1, alpha=.5)+
  guides(fill="none", colour="none")+
  facet_wrap(~cluster, ncol=3, scales = "fixed")+theme_light()
#we can use the average of 1000 times SOMs weight for the initial weight of the next som calculation.


#the worst
ggplot(data.frame(Precision=Prob_AllSom_K[[603]]), aes(x=Precision)) + geom_density(adjust=2, fill="#2c8ba1") +guides(fill="none", colour="none")+theme_light()
#the best
ggplot(data.frame(Precision=Prob_AllSom_K[[146]]), aes(x=Precision)) + geom_density(adjust=2, fill="#5d8aa8") +guides(fill="none", colour="none")+theme_light()
#together
all_dist<-do.call(rbind,lapply(c(603,146), \(x) data.frame(prob=Prob_AllSom_K[[x]], num=as.character(x))))
ggplot(all_dist, aes(x=prob, fill=num)) + geom_density(adjust=2, alpha=.2) +guides(fill="none", colour="none")+theme_light()

which.min(unlist(lapply(Prob_AllSom_K,mean))) #max:146, min 603




all_dist<-do.call(rbind,lapply(1:30, \(x) data.frame(Precision=Prob_AllSom_K[[x]], num=x))) %>% arrange(num)
all_dist$Precision[all_dist$Precision<0.55]<-1

ggplot(all_dist, aes(x=Precision)) + 
  geom_density(fill="#5d8aa8", adjust=2, alpha=.5)+
  guides(fill="none", colour="none")+
  facet_wrap(~num, ncol=5)+
  theme_light()#+ theme(strip.background = element_blank(), strip.placement = "outside")


all_dist<-do.call(rbind,lapply(1:1000, \(x) data.frame(Precision=Prob_AllSom_K[[x]], num=as.character(x))))

ggplot(all_dist, aes(x=Precision, fill=num)) + geom_density(adjust=2, alpha=.2) +guides(fill="none", colour="none")+theme_light()



ggplot(all_dist, aes(x=prob, fill=num)) + geom_density(adjust=2, alpha=.2) +guides(fill="none", colour="none")

ggplot(all_dist, aes(x=prob, fill=num)) + geom_histogram(binwidth=.0005, alpha=.5, position="identity")+guides(fill="none", colour="none")

indiv<-dist<-data.frame(prob=Prob_AllSom_K[[19]])
ggplot(indiv, aes(x=prob)) + geom_histogram(binwidth=.001, alpha=.5, position="identity")+guides(fill="none", colour="none")



#=============================
#latex
m=paste0(t(cbind(confusion_mat[[723]],confusion_mat[[965]])))
for(i in 1:length(m)){
  if(i%%12!=0){
    m[i]<-paste0(m[i]," & ")
  }
  else{
  m[i]<-paste0(m[i]," \\\\ \n")
  }
}
mm<-paste(m, collapse = "")
write_file(mm, "/Volumes/DiscB/PHD/RCODE/Routputs/mat.txt")

#=============================
#latex
m=paste0(t(round(t(apply(all_ratios, 2, function(x)quantile(x,c(0.05,0.1,0.2,0.3,0.4,0.5,0.75,0.9,1)))),4)))
for(i in 1:length(m)){
  if(i%%9!=0){
    m[i]<-paste0(m[i]," & ")
  }
  else{
    m[i]<-paste0(m[i]," \\\\ \n")
  }
}
mm<-paste(m, collapse = "")
write_file(mm, "/Volumes/DiscB/PHD/RCODE/Routputs/matporc.txt")
utils::browseURL("/Volumes/DiscB/PHD/RCODE/Routputs/matporc.txt")
