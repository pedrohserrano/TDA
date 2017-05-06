#CLUSTERING

#----------------------------- LIBRARIES -----------------------------
#install.packages("fpc")
#install.packages("kernlab")
#install.packages("dbscan")
#install.packages("jsonlite")
#install.pachkages("qrage")
library(fpc)
library(dbscan)
library(kernlab)
library(jsonlite)
library(igraph)
library(RColorBrewer)
library(qrage)

#----------------------------- LOAD DATA -----------------------------
rm(list = ls())
rutawork = ('/home/denny/itam/topologia/proyecto_clustering/')
datos <- read.csv(paste(rutawork,'ecobici_preprocessed.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
str(datos)

#por lo pronto, me quedo con las variables int
datos <- datos[c("Edad_Usuario", "Distancia_km", "Duracion_viaje")]

#----------------------------- PCA ----------------------------------
#Esta parte nos ayuda a seleccionar la variable mas significativa que explica la dinamica de los datos

#PCA
# datos.pca <- prcomp(datos ,center = TRUE,scale. = TRUE)
# summary(datos.pca)
# str(datos.pca)
# datos.pca$rotation
#en este caso observarmos que la variable mas importante es distancia_km


#KERNEL PCA
variables <- length(datos)
sigma <- 1
kres <- kpca(~., data=datos,features=variables,kernel="rbfdot",kpar = list(sigma = sigma))
data_kpca <- as.data.frame(kres@rotated)
#ordena de mayor a menor (por la que tiene mayor varianza)
datos_prueba <- data_kpca[c("V1","V2")]

#----------------------------- GENERATE INTERVALS----------------------------------

df <- datos_prueba      #choose a dataset

#----------------------------- NECESSARY PARAMETERS -----------------------------
#var_o <- data$x1    #variable we will use to make the overlapping subsets
var_o <- datos_prueba$V1   #if we want to use kernel pca variable to cut
n_int <- 6  #number of intervals we want
p <- 0.2  #proportion of each interval that should overlap with the next

#----------------------------- CREATING THE INTERVALS -----------------------------

#this section will create a data frame in which we will construct overlapping intervals
intervals_centers <- seq(min(var_o),max(var_o),length=n_int)  #basic partition = centers
interval_length <- intervals_centers[2]-intervals_centers[1]  #to create the overlaps of p% of this length
intervals <- data.frame(centers=intervals_centers)            #create a data frame
#create the overlapping intervals  
intervals$min <- intervals_centers - (0.5+p)*interval_length                     
intervals$max <- intervals_centers + (0.5+p)*interval_length
intervals$interval <- seq(1,n_int)
intervals$name <- with(intervals, sprintf("[%.2f;%.2f)",min,max))

#function that will split the variable according to the invervals
res <- lapply(split(intervals,intervals$interval), function(x){   
  return(df[var_o> x$min & var_o <= x$max,])     #res will be a list with each element res[i]
})                                               #being the points on the i'th subset

#res

#----------------------------- CLUSTERING FOR EACH INTERVAL -----------------------------

###-------------------------FUNCTIONS------------------------------------------------
## We obtain a clustering that attemps to not depend on the eps parameter
## we give a maximum eps for attempting clustering and evaluate on the percentage of 
## noise obtained
noClust<-function(data, eps=0.7, eps_cl=6.5, np=.1){
  #Default parameters for dbscan
  p_noise <- 0.05       # we use this as a rule of thumb
  ##Number of clusters detected
  numClust<-0
  ##Noise percentage
  noise_perc<-1
  MinPts <- p_noise*dim(data)[1]
  # We iterate eps through an geometric function starting on the given value
  for(j in 0:10){
    eps<-eps+j*eps
    ## We iterate also on the eps_cl parameter with an exponential function
    for(i in 0:3 ){
      result<-optics(data,eps=eps,minPts=MinPts,eps_cl = eps_cl*10**-i)
      noise_perc=length(result$cluster[result$cluster==0])/length(result$cluster[result$cluster!=0])
      if (noise_perc < np) {
        numClust<-max(result$cluster)
        return(list(cluster=result$cluster, noise_perc=noise_perc, num_clust=numClust))
      }
    }
  }
  list(cluster=result$cluster, noise_perc=noise_perc)
}
p_noise <- 0.05
#ITERATE EVERY ELEMENT OF THE LIST (res[i]) AND CLUSTERIZE INSIDE
ints<-list()
counter1<-1;counter2<-1

for(i in 1:(n_int-1)){
  df1<-as.data.frame(res[[i]])
  df2<-as.data.frame(res[[i+1]])
  
  if(i==1){
    MinPts <- p_noise*dim(df1)[1]
    result1<-(noClust(df1))
    df1$cluster1 <- result1$cluster
    
    #create columns in the original matrix to show which cluster they belong to
    df[dim(df)[2]+i]<-rep(0,dim(df)[1])
    df[row.names(df1),dim(df)[2]]<-result1$cluster
    
  }else{result1 <- result2              #use the results for the last iteration
  df1$cluster1 <- result1$cluster #this ensures that the cluster labels will be correct for the adj. matrix
  }
  
  MinPts <- p_noise*dim(df2)[1]
  result2<-(noClust(df2))
  df2$cluster2 <- result2$cluster
  
  #create columns in the original matrix to show which cluster they belong to
  df[dim(df)[2]+1]<-rep(0,dim(df)[1])
  df[row.names(df2),dim(df)[2]]<-result2$cluster
  
  intersection <- merge(df1,df2,all=TRUE)            #points in the intersection
  intersection[is.na(intersection)] <- 0
  ints[[i]]<-as.data.frame(unique(intersection[3:4]))               #list of all the clusters that intersect
  
}

#plot(df$V1,df$V2)


#---------------------DISTANCE BETWEEN CLUSTERS AND ADJ_MATRIX -----------------------------

#Bajo el supuesto que leemos una base que tiene columnas para cada intervalo
#Leemos la base de datos con los intervalos
base <- df
#Creamos una columna para los clusters
base$clusters<-0

#Columna en donde empieza el intervalo 1:
int_ini <- 3
#Columna en donde se ubica el ultimo intervalo:
int_fin <- 8
#Columna donde se creo la columna de "clusters":
col_cluster <- 9


for(i in seq(nrow(base[,int_ini:int_fin]))){
  temp<-c()
  for(m in seq(int_ini,int_fin)){
    if (base[i,m] > 0){
      
      temp<-c(paste0(temp,base[i,m],sep = ","))
    }
  }
  if(length((temp))>0){
    aux<-unlist(strsplit(temp,","))
    aux2<-unique(aux)
    aux3<-paste(aux2,collapse=",")
    base[i,col_cluster]<-aux3
  }
}

base <- data.frame(base$clusters)
names(base) <- c("clusters")

#Creamos una variable para enumerar la observacion
base$obs <- paste("obs",seq(1,length(base$clusters)))

#Detectamos los clusters
num_clusters <- sort(unique(strsplit(paste0(base$clusters, collapse=","),",")[[1]]))
clusters <- length((num_clusters))

#Creamos una columna con ceros para cada cluster
for(x in num_clusters){
  base[[paste("c",x,sep="_")]] <- rep(0,nrow(base))
}

#Para cada columna que creamos agregamos un 1 según el cluster al que pertenece la obs
base$clusters<- as.character(base$clusters)


#Ojo es x+3, porque en la columna 3 en adelante es donde va vaciar los "1" de cada cluster
for(i in seq(nrow(base))){
  vector <- strsplit(base$clusters[i], ",")[[1]]
  vector <- sort(as.numeric(vector))
  for(x in vector){
    base[i,(x+3)] <- 1
  }
}


dummy_mat<-base[,3:ncol(base)]
n_clusters<-ncol(dummy_mat)
cluster_list<-lapply(1:n_clusters,function (col) {which(dummy_mat[,col]==1)})

adj_matrix<-matrix(0,nrow=n_clusters,ncol=n_clusters)

for(i in 1:(n_clusters-1)){
  for(j in (i+1):n_clusters){
    distancia<-setdiff(cluster_list[[i]], cluster_list[[j]])
    cercania<-length(cluster_list[[i]])-length(distancia)
    adj_matrix[i,j]<-round(cercania/min(length(cluster_list[[i]]),length(cluster_list[[j]])),2)
    adj_matrix[j,i]<-adj_matrix[i,j]
  }
}

summary_cluster<-matrix(0,nrow=1,ncol=n_clusters)
for(i in 1:n_clusters){
  summary_cluster[1,i]<-length(cluster_list[[i]])
}


#KEPLER
nodes.n <- clusters
nodes.size<- as.numeric(summary_cluster)/100
nodes.tooltips <- paste("Grupo:", 1:nodes.n)
nodes.names <- 1:nodes.n
nodes.color <- as.character(1:nodes.n)

# ------- AHORA TENEMOS QUE CREAR UN JSON DE ESO -----------------------------
adj.matrix <- adj_matrix
aux_mat <- data.frame()
for(i in 1:nodes.n) for(j in 1:nodes.n) if(adj.matrix[i, j]!=0) aux_mat <- rbind(aux_mat, data.frame(source=i-1, target=j-1, value=adj.matrix[i, j]))
linksJSON <- toJSON(aux_mat)
nodesJSON <- toJSON(data.frame(color=nodes.color, group=nodes.size, name=nodes.names, tooltip=nodes.tooltips))
graphJSON <- sprintf("{\"nodes\": %s, \"links\": %s}", nodesJSON, linksJSON)
#head(graphJSON)

# ------------  CREAMOS EL HTML ----------------------------------------------------------
htmlFile <- readLines('/home/denny/itam/topologia/ManifoldLearning/www/index.html')
#htmlFile <- readLines("www/index.html")
graph_def_line <- which(grepl("graph =", htmlFile))
#htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
#writeLines(htmlFile, "www/index.html")
writeLines(htmlFile, '/home/denny/itam/topologia/ManifoldLearning/www/index.html')
browseURL("file:////home/denny/itam/topologia/ManifoldLearning/www/index.html")


#----------------------------OTRA VISUALIZACION IGRAPH ---------------------------------------

colores <- brewer.pal(12,"Set3")
grafica <- function(mat,colores,tam=c(1),lay='kamada') {
  if(tam==c(1)) {
    tam <- rep(1,dim(mat)[1])
  }
  g <- graph.adjacency(mat,mode='undirected',weighted=T)
  g <- simplify(g)
  if(lay=='kamada') {
    plot(g,vertex.color=colores,vertex.size=(tam/sum(mat)),edge.arrow.size=.3,rescale=F, layout=layout.kamada.kawai)
  }
}
grafica(adj_matrix,colores,summary_cluster,'kamada')


#------------------------OTRA VISUALIZACION QRAGE------------------------

tofdg <- function(matriz){
  a <- as.data.frame(matrix(0,nrow=dim(matriz)[1]**2/2,ncol=3))
  contador <- 0
  for(i in 1:dim(matriz)[2]){
    col <- matriz[,i]
    for(j in 1:i){
      a[contador+1,3] <- as.numeric(matriz[i,j])
      a[contador+1,2] <- j
      a[contador+1,1] <- i
      contador <- contador +1
    }
  }
  a
}

x <- tofdg(adj.matrix)
z <- as.data.frame(cbind(seq(1,clusters,1),t(summary_cluster)))
colores2 <- as.data.frame(cbind(seq(1,clusters,1),t(colores)[1:clusters]))


qrage(links=x, width = 1000, height = 800,distance=8000,nodeValue=z
      ,nodeColor=colores2,linkColor='#00f',arrowColor='#f00'
      ,cut=0.01,textSize=12
      ,linkWidth=c(1,8),linkOpacity=c(0.6,1))


