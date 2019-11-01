# 1. For each selection method p and clustering method, produce a sequence of statistics values for each cluster size for.
# Show this on a table and graph.
# 2. Produce a Summary table of optimal number of clusters.
require(cluster)
# max clusters
K_max=20;
# data
df <- USArrests
df <- na.omit(df)
df <- scale(df)
d <- dist(df, method = "euclidean")
# Hierarchical clustering using ward.D
hc1 <- hclust(d, method = "ward.D" )
hclust1 <- function(data,k){ list(cluster=cutree(hclust(dist(data),method="ward.D"),k))     }

# selection function 
# 1.gap
# clusGap(df, FUN = hclust1, K.max = K_max, B = 60)

# 2.second derivative
# x is an array;dd number of differences
clust_sel = function(x,jrange=3:K_max,dd=2) {
  df=x;
  df <- na.omit(df);
  df <- scale(df);
  d <- dist(df, method = "euclidean");
  y <-hclust(d, method = "ward.D" );
  win_ss = function(x,y,w = rep(1, length(y))) sum(lm(x~factor(y),weights = w)$resid^2*w)
  ### win_ss calculates within cluster sum of squares
  sm1 = NULL
  for(i in jrange) sm1[i] = win_ss(x,cutree(y,i)) 
  sm1=sm1[jrange]
  k = if(dd==1) sm1[-1] else -diff(sm1) 
  plot( jrange[-length(k)+1:0] , -diff(k)/k[-length(k)]*100, xlab = "nclust",ylab = "change")  #choose max.
  
  jrange [sort.list(diff(k)/k[-length(k)]*100)[1]]
  cat("second derivative statistics (3,4,5..K): ",-diff(k)/k[-length(k)]*100)
  cat("\n optimal sec_der_cluster",jrange [sort.list(diff(k)/k[-length(k)]*100)[1]])
}

# 3.silhouette 
# return ("cluster", "neighbor", "sil_width"); close to 1 means that the data is appropriately clustered
sil_sel=function(data,jrange=1:K_max){
  d <- dist(data, method = "euclidean")
  hc1 <- hclust(d, method = "ward.D" )
  Y <- array(0,dim=c(1,K_max))
  for(i in jrange){ 
    sil_all<-silhouette(cutree(hc1, h=i) ,d,FUN=hclust1(data,i), title=title(main = 'test'))
    summary(sil_all)$avg.width
    Y[i]=summary(sil_all)$avg.width ;
  }
  optimal=max(Y);optimal_ind=which.max(Y);
    
  sil_all<-silhouette(cutree(hc1, h=optimal_ind) ,d,FUN=hclust1(df,optimal_ind), title=title(main = 'test'))
  plot(sil_all)
  cat("optimal_sil_sel:",optimal," ");
  cat(" optimal cluster:",optimal_ind);
}

# Plot the obtained dendrogram
par(mfrow=c(3,3))
plot(hc1, cex = 0.6, hang = -1,xlab = "Ward's Method")
# clustering...
clusGap(df, FUN = hclust1, K.max = K_max, B = 60)
clust_sel(df,jrange=3:K_max,dd=2)
sil_sel(df,jrange=1:K_max)
