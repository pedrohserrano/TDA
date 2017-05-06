library(fpc)
library(dbscan)
library(kernlab)
library(jsonlite)
library(RColorBrewer)

#----------------------------- LOAD DATA -----------------------------
#we define 3 num_variables to begin, and loading data from s3
rm(list = ls())
rutawork <- ('/Users/pedrohserrano/TDA_proyectos/proyecto_clustering/')
#consultas <- read.csv(paste(rutawork,'p7-consultas-yza-integrado.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
consultas <- read.csv(paste(rutawork,'p7-med-yza-join-cronicos.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
str(consultas)
#consultas<-consultas[c("CONSULTAS", "SERVICIOS", "REDENCION")]
consultas<-consultas[c("ID_USUARIO", "DOSIS", "FRECUENCIA")]

#----------------------------- CLUSTERING METHOD-----------------------------
#At the very begining we define the clustering method we gonna use later, this time breath fisrt search
clusterBFS <- function(distmat, neigh){
  mat <- distmat
  dimnames(mat)=list(1:nrow(distmat),1:nrow(distmat))
  clustList <- list()
  while(nrow(mat)!=0){
    indices <- which(bfs_density_component(mat, nbs=neigh))
    if((nrow(mat)-length(indices))!=1){
      individuos <- row.names(mat)[indices]
      mat <- mat[-indices, -indices]
      clustList <- c(clustList, list(individuos))  
    } else{
      el_we_solito <- row.names(mat)[!(row.names(mat) %in% row.names(mat)[indices])]
      individuos <- row.names(mat)[indices]
      clustList <- c(clustList, list(individuos), list(el_we_solito))
      break
    }
  }
  return(clustList)
}

#----------------------------- KERNEL PCA ----------------------------------
#We are going to show what's the most important variable to understand data
num_variables <- length(consultas)
sigma <- 1
kernel_response <- kpca(~., data=consultas,features=num_variables,kernel="rbfdot",kpar = list(sigma = sigma))
kpca_vectors <- as.data.frame(kernel_response@rotated) #desc by variance
kpca_vectors <- kpca_vectors[c("V1","V2")]

#----------------------------- CREATING THE INTERVALS -----------------------------
#We will create a data frame in which we will construct overlapping intervals, then over intervals build clusters
vec_1 <- kpca_vectors$V1   #if we want to use kernel pca variable to cut
n_int <- 6  #number of intervals we want
p <- 0.2  #proportion of each interval that should overlap with the next

intervals_centers <- seq(min(vec_1),max(vec_1),length=n_int)  #basic partition = centers
interval_length <- intervals_centers[2]-intervals_centers[1]  #to create the overlaps of p% of this length
intervals <- data.frame(centers=intervals_centers)            #create a data frame

#create the overlapping intervals  
intervals$min <- intervals_centers - (1-p)*interval_length                     
intervals$max <- intervals_centers + (1-p)*interval_length
intervals$interval <- seq(1,n_int)
intervals$name <- with(intervals, sprintf("[%.2f;%.2f)",min,max))

#function that will split the variable according to the invervals
res <- lapply(split(intervals,intervals$interval), function(x){   
  return(kpca_vectors[vec_1> x$min & vec_1 <= x$max,])     #res will be a list with each element res[i]
})                                               #being the points on the i'th subset

kpca_vectors$Clusters<-numeric(nrow(kpca_vectors))
for(i in 1:length(res)) {
  if(nrow(res[[i]])!=0){
    interval_temp<-res[[i]]
    distmat <- as.matrix(dist(interval_temp))
    neigh <- 6
    clusters_i<-clusterBFS(distmat, neigh) #Here lays the clustering method!!
    interval_temp$Cluster<-numeric(nrow(interval_temp))
  
  for(j in 1:length(clusters_i)){
    for(z in 1:length(clusters_i[[j]])){
      interval_temp$Cluster[as.numeric(clusters_i[[j]][z])]<-j
    }
  }
  
  res[[i]]<-interval_temp
  }
}

for(i in 1:length(res)){
  
  kpca_vectors[[ncol(kpca_vectors)+1]]<-numeric(nrow(kpca_vectors))
  interval_temp<-res[[i]]
  interval_temp$filas<-rownames(interval_temp)
  for(j in 1:nrow(kpca_vectors)){
    
    if(rownames(kpca_vectors)[j] %in% rownames(interval_temp)){
      kpca_vectors[[ncol(kpca_vectors)]][j]<-interval_temp$Cluster[which(interval_temp$filas==rownames(kpca_vectors)[j])]
    }#else{
#       kpca_vectors[[ncol(kpca_vectors)+1]][j]<-0
#     }
    
  }
  
}

#----------------------------- CLUSTER TAGGING -----------------------------
#We need to know on what cluster every element belong to
for(i in 1:length(res)){
 
  ncol<-ncol(kpca_vectors)
  kpca_vectors[[ncol+1]]<-numeric(nrow(kpca_vectors))
  ncol<-ncol(kpca_vectors)
  
  if(i==1){
    clusters_counter<-0
    kpca_vectors[[ncol]]<-kpca_vectors[[3+i]]
  }else{
    
    for(j in 1:nrow(kpca_vectors)){
      
      if(kpca_vectors[[3+i]][j]!=0){
        
        kpca_vectors[[ncol]][j]<-kpca_vectors[[3+i]][j]+clusters_counter
        
        }
      
      }
    
    }
  clusters_counter<-clusters_counter+max(kpca_vectors[[3+i]])
}

kpca_vectors<-kpca_vectors[,-c(3:8)]


#---------------------DISTANCE BETWEEN CLUSTERS AND ADJ_MATRIX -----------------------------
#We read a dataset and it have columns for each interval, then we read kpca_vectors de consultas con los intervalos

kpca_vectors$clusters<-0 #Creamos una columna para los clusters
int_fin <- ncol(kpca_vectors)-1 #column where the last interval belongs
int_ini <- int_fin-n_int+1 #column where the interval 1 starts
col_cluster <- ncol(kpca_vectors) #column where we create the "clusters"

for(i in seq(nrow(kpca_vectors[,int_ini:int_fin]))){
  temp<-c()
  for(m in seq(int_ini,int_fin)){
    if (kpca_vectors[i,m] > 0){
      
      temp<-c(paste0(temp,kpca_vectors[i,m],sep = ","))
    }
  }
  if(length((temp))>0){
    aux<-unlist(strsplit(temp,","))
    aux2<-unique(aux)
    aux3<-paste(aux2,collapse=",")
    kpca_vectors[i,col_cluster]<-aux3
  }
}

kpca_vectors <- data.frame(kpca_vectors$clusters)
names(kpca_vectors) <- c("clusters")

#Creamos una variable para enumerar la observacion
kpca_vectors$obs <- paste("obs",seq(1,length(kpca_vectors$clusters)))

#Detectamos los clusters
num_clusters <- sort(as.numeric(unique(strsplit(paste0(kpca_vectors$clusters, collapse=","),",")[[1]])))
a<-strsplit(paste0(kpca_vectors$clusters, collapse=","),",")
clusters <- length((num_clusters))

#Creamos una columna con ceros para cada cluster
for(x in num_clusters){
  kpca_vectors[[paste("c",x,sep="_")]] <- rep(0,nrow(kpca_vectors))
}

#Para cada columna que creamos agregamos un 1 seg??n el cluster al que pertenece la obs
kpca_vectors$clusters<- as.character(kpca_vectors$clusters)

#Ojo es x+3 (int_ini), porque en la columna 3 en adelante es donde va vaciar los "1" de cada cluster
for(i in seq(nrow(kpca_vectors))){
  vector <- strsplit(kpca_vectors$clusters[i], ",")[[1]]
  vector <- sort(as.numeric(vector))
  for(x in vector){
    kpca_vectors[i,(x+int_ini)] <- 1
  }
}

if(sum(kpca_vectors[[3]])==0){
  dummy_mat<-kpca_vectors[,4:ncol(kpca_vectors)]
}else{dummy_mat<-kpca_vectors[,3:ncol(kpca_vectors)]}

n_clusters<-ncol(dummy_mat)
cluster_list<-lapply(1:n_clusters,function (col) {which(dummy_mat[,col]==1)})

adj_matrix<-matrix(0,nrow=n_clusters,ncol=n_clusters)

for(i in 1:(n_clusters-1)){
  # i<-1 
  for(j in (i+1):n_clusters){
    # j<-i+1
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


#---------------------TOOLTIPS-----------------------------
library(dplyr)
master <- cbind(consultas,dummy_mat)
ww <- paste0("c_",seq(2,clusters,1))

bb <- data.frame()
for(x in 1:length(ww)){
  myCols <- (c("Consultas","Servicios", ww[x]))
  colNums <- match(myCols,names(master))
  a  <- (master  %>% select(colNums) %>% filter(master[3+x]=='1') %>% 
           summarise(dist_prom=mean(CONSULTAS),duracion_prom=mean(SERVICIOS)))
  bb <- rbind.data.frame(bb,a)
}
#print(bb)


#install.packages('gplots')
library('gplots')
#vamos a ordenar los colores segun la variable "var_color"
var_color <- bb$dist_prom
a <- col2hex(colorRampPalette(c("blue", "red"))(110))
a <- a[order(var_color)]
c <- paste(shQuote(a,type="cmd"), collapse=", ")
b <- paste(shQuote(order(var_color),type="cmd"), collapse=", ")
##############################


#---------------------KEPLER MAPPER-----------------------------
nodes.n <- clusters
nodes.size<- as.numeric(summary_cluster[,-1])#/100
nodes.tooltips <- paste("Grupo:", 1:nodes.n)
nodes.names <- 1:nodes.n
nodes.color <- as.character(1:nodes.n)


# ------- AHORA TENEMOS QUE CREAR UN JSON -----------------------------
adj.matrix <- adj_matrix
aux_mat <- data.frame()
for(i in 1:nodes.n) {
  
  for(j in 1:nodes.n){
    
  if(adj.matrix[i, j]!=0) {

        aux_mat <- rbind(aux_mat, data.frame(source=i-1, target=j-1, value=adj.matrix[i, j]))
    }
  }
}

linksJSON <- toJSON(aux_mat)
nodesJSON <- toJSON(data.frame(color=nodes.color, group=nodes.size, name=nodes.names, tooltip=nodes.tooltips))
graphJSON <- sprintf("{\"nodes\": %s, \"links\": %s}", nodesJSON, linksJSON)

# ------------  CREAMOS EL HTML ----------------------------------------------------------
htmlFile <- readLines('/Users/pedrohserrano/TDA_proyectos/ManifoldLearning/www/index.html')
graph_def_line <- which(grepl("graph =", htmlFile))
htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
writeLines(htmlFile, '/Users/pedrohserrano/TDA_proyectos/ManifoldLearning/www/index.html')



browseURL("file:////Users/pedrohserrano/TDA_proyectos/ManifoldLearning/www/index.html")


# ------------  CREAMOS EL HTML ----------------------------------------------------------
# htmlFile <- readLines('/home/denny/itam/topologia/ManifoldLearning/www/index.html')
# #htmlFile <- readLines("www/index.html")
# graph_def_line <- which(grepl("graph =", htmlFile))
# #htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
# htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
# #writeLines(htmlFile, "www/index.html")
# writeLines(htmlFile, '/home/denny/itam/topologia/ManifoldLearning/www/index.html')
# ######colores FALTA CAMBIAR COLOR POR VARIABLE
# #modificar el index
# setwd('..')
# setwd("ManifoldLearning/www/")
# htmlFile <- readLines('index.html')
# 
# numbers_line <- which(grepl("domain", htmlFile))
# colors_line <- which(grepl("range", htmlFile))
# htmlFile[numbers_line] <- paste0('.domain([',b,']);')
# htmlFile[numbers_line+1] <- paste0('.range([',c,']);')
# writeLines(htmlFile,'index.html')
# ###
# browseURL("file:////home/denny/itam/topologia/ManifoldLearning/www/index.html")


