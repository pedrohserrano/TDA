library(fpc)
library(dbscan)
library(kernlab)
library(jsonlite)
library(igraph)
library(RColorBrewer)
library(qrage)

#----------------------------- LOAD DATA -----------------------------
rm(list = ls())
#rutawork = ('/home/denny/itam/topologia/proyecto_clustering/')
rutawork = ('/Users/pedrohserrano/TDA_proyectos/proyecto_clustering/')
#datos <- read.csv(paste(rutawork,'ecobici_preprocessed.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
datos <- read.csv(paste(rutawork,'consultas-p7-yza-integrado.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
datos_origin<-datos
str(datos)

# proporcion_entrena<-1
# inTraining <- createDataPartition(datos_origin$Genero_Usuario, p = proporcion_entrena, list = FALSE)
# datos <- datos_origin[ inTraining,]
# 
# #por lo pronto, me quedo con las variables int
# datos <- datos[c("Edad_Usuario", "Distancia_km", "Duracion_viaje")]
datos<-datos_origin[c("CONSULTAS", "SERVICIOS", "REDENCION")]

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
print("WEEEEEEEEEY  YA QUED?? EL PINCHE PCA")
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

df$Clusters<-numeric(nrow(df))
for(i in 1:length(res)) {
  # i<-1
  if(nrow(res[[i]])!=0){
    interval_temp<-res[[i]]
    distmat <- as.matrix(dist(interval_temp))
    neigh <- 6
    clusters_i<-clustITAM(distmat, neigh)
    interval_temp$Cluster<-numeric(nrow(interval_temp))
  
  for(j in 1:length(clusters_i)){
    #j<-1
    
    for(z in 1:length(clusters_i[[j]])){
    #  z<-1
      interval_temp$Cluster[as.numeric(clusters_i[[j]][z])]<-j
      #df$Clusters[as.numeric(rownames(interval_temp)[nrow(interval_temp)])]<-j
      
    }
  
  }
  
  res[[i]]<-interval_temp
  }
}


for(i in 1:length(res)){
  
  df[[ncol(df)+1]]<-numeric(nrow(df))
  interval_temp<-res[[i]]
  interval_temp$filas<-rownames(interval_temp)
  for(j in 1:nrow(df)){
    
    if(rownames(df)[j] %in% rownames(interval_temp)){
      df[[ncol(df)]][j]<-interval_temp$Cluster[which(interval_temp$filas==rownames(df)[j])]
    }#else{
#       df[[ncol(df)+1]][j]<-0
#     }
    
  }
  
}

####---------Empieza reetiquetamiento de Clusters----------####
df_<-df
# df_<-df_salvado

for(i in 1:length(res)){
 
  ncol<-ncol(df_)
  df_[[ncol+1]]<-numeric(nrow(df_))
  ncol<-ncol(df_)
  
  if(i==1){
    
    contador_clusters<-0
    df_[[ncol]]<-df_[[3+i]]
    
  }else{
    
    for(j in 1:nrow(df_)){
      
      if(df_[[3+i]][j]!=0){
        
        df_[[ncol]][j]<-df_[[3+i]][j]+contador_clusters
        
        }
      
      }
    
    }
  contador_clusters<-contador_clusters+max(df_[[3+i]])
  
}

df_salvado<-df
df<-df_[,-c(3:8)]#hasta 9
#df<-df_[,-c(3:9)]

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#---------------------DISTANCE BETWEEN CLUSTERS AND ADJ_MATRIX -----------------------------

#Bajo el supuesto que leemos una base que tiene columnas para cada intervalo
#Leemos la base de datos con los intervalos
base <- df
#Creamos una columna para los clusters
base$clusters<-0


#Columna en donde se ubica el ultimo intervalo:
int_fin <- ncol(base)-1
#Columna en donde empieza el intervalo 1:
int_ini <- int_fin-n_int+1
#Columna donde se creo la columna de "clusters":
col_cluster <- ncol(base)


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
num_clusters <- sort(as.numeric(unique(strsplit(paste0(base$clusters, collapse=","),",")[[1]])))
a<-strsplit(paste0(base$clusters, collapse=","),",")
clusters <- length((num_clusters))

#Creamos una columna con ceros para cada cluster
for(x in num_clusters){
  base[[paste("c",x,sep="_")]] <- rep(0,nrow(base))
}

#Para cada columna que creamos agregamos un 1 seg??n el cluster al que pertenece la obs
base$clusters<- as.character(base$clusters)


#Ojo es x+3 (int_ini), porque en la columna 3 en adelante es donde va vaciar los "1" de cada cluster
for(i in seq(nrow(base))){
  vector <- strsplit(base$clusters[i], ",")[[1]]
  vector <- sort(as.numeric(vector))
  for(x in vector){
    base[i,(x+int_ini)] <- 1
  }
}

if(sum(base[[3]])==0){
  dummy_mat<-base[,4:ncol(base)]
}else{dummy_mat<-base[,3:ncol(base)]}

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


#######################
library(dplyr)
master <- cbind(datos,dummy_mat)
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


#KEPLER
nodes.n <- clusters
nodes.size<- as.numeric(summary_cluster[,-1])#/100
nodes.tooltips <- paste("Grupo:", 1:nodes.n)
nodes.names <- 1:nodes.n
nodes.color <- as.character(1:nodes.n)

# ------- AHORA TENEMOS QUE CREAR UN JSON DE ESO -----------------------------
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
#head(graphJSON)

# ------------  CREAMOS EL HTML ----------------------------------------------------------
# htmlFile <- readLines('/home/denny/itam/topologia/ManifoldLearning/www/index.html')
# htmlFile <- readLines('/Users/jalfredomb/Dropbox/01_ITAM_Ciencia_de_Datos/2do_semestre/Topologia/TDA_proyectos//ManifoldLearning/www/index.html')
htmlFile <- readLines('/Users/pedrohserrano/TDA_proyectos/ManifoldLearning/www/index.html')
#htmlFile <- readLines("www/index.html")
graph_def_line <- which(grepl("graph =", htmlFile))
#htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
#writeLines(htmlFile, "www/index.html")
#writeLines(htmlFile, '/home/denny/itam/topologia/ManifoldLearning/www/index.html')
#writeLines(htmlFile,'/Users/jalfredomb/Dropbox/01_ITAM_Ciencia_de_Datos/2do_semestre/Topologia/TDA_proyectos//ManifoldLearning/www/index.html')
writeLines(htmlFile, '/Users/pedrohserrano/TDA_proyectos/ManifoldLearning/www/index.html')
#browseURL("file:////Users/jalfredomb/Dropbox/01_ITAM_Ciencia_de_Datos/2do_semestre/Topologia/TDA_proyectos//ManifoldLearning/www/index.html")


browseURL("file:////Users/pedrohserrano/TDA_proyectos/ManifoldLearning/www/index.html")


# ------------  CREAMOS EL HTML ----------------------------------------------------------
htmlFile <- readLines('/home/denny/itam/topologia/ManifoldLearning/www/index.html')
#htmlFile <- readLines("www/index.html")
graph_def_line <- which(grepl("graph =", htmlFile))
#htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
#writeLines(htmlFile, "www/index.html")
writeLines(htmlFile, '/home/denny/itam/topologia/ManifoldLearning/www/index.html')
######colores FALTA CAMBIAR COLOR POR VARIABLE
#modificar el index
setwd('..')
setwd("ManifoldLearning/www/")
htmlFile <- readLines('index.html')

numbers_line <- which(grepl("domain", htmlFile))
colors_line <- which(grepl("range", htmlFile))
htmlFile[numbers_line] <- paste0('.domain([',b,']);')
htmlFile[numbers_line+1] <- paste0('.range([',c,']);')
writeLines(htmlFile,'index.html')
###
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
