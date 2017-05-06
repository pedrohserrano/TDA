# CLUSTERING ALGORITH!!!!
neigh <- 6
distmat <- as.matrix(dist(iris[ ,1:4]))

clustITAM <- function(distmat, neigh){
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

clustITAM(distmat, neigh)
