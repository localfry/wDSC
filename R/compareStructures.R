compareStructures <- function(structures, method=NULL, hausdorff.method=NULL, verbose=TRUE, plot=TRUE, pixels=100, dose) {
  if (class(structures) != "structure.list") {
    warning("Input 'structures' must be of class 'structure.list'")
    return()
  }
  empty <- unlist(lapply(structures, function(struct) {return(dim(struct)[1] <= 0)}))
  if (any(empty)) {
    warning(paste("Skipping empty structure(s): ", paste(names(structures[empty]), collapse=", ", sep=""), sep=""))
    structures <- structures[!empty]
  }
  N <- length(structures)
  if (N < 2) {
    warning("Need at least 2 structures to perform comparison")
    return()
  }
  method <- match.arg(method, choices=c("axial", "surface", "hausdorff", "grid", "DSC", "EMD", "wDSC"))
  hausdorff.method <- match.arg(hausdorff.method, choices=c("mean", "median", "absolute"))
  switch(method,
         axial = contours <- compareStructures.axial(structures, pixels=pixels),
         grid = {
           warning("method='grid' no longer supported, use method='axial' instead")
           contours <- compareStructures.axial(structures, pixels=pixels)
         },
         surface = contours <- compareStructures.surface(structures),
         hausdorff = return(compareStructures.hausdorff(structures, method=hausdorff.method, verbose=verbose)),
         DSC = {
           contours <- compareStructures.axial(structures, pixels=pixels)
           N <- dim(contours)[2]-3
           results <- matrix(0, nrow=N, ncol=N, dimnames=list(names(structures), names(structures)))
           for (i in 1:N) {
             for (j in 1:N) {
               if (i == j) {
                 results[i, j] <- 1#
                 next
               }
               results[i, j] <- 2*sum((contours[,i+3]>0) & (contours[,j+3]>0))/(sum(contours[,i+3]>0)+sum(contours[,j+3]>0))
             }
           }
           return(results)
         },
         EMD = return(compareStructures.EMD(structures)),
         wDSC = return(compareStructures.wDSC(structures, dose, pixels))
  )
  if ((plot) & (method %in% c("axial", "grid"))) {
    mar.old <- par()$mar
    par(mar=c(0, 0, 0, 0))
    z.unique <- sort(unique(contours[,3]))
    N.z <- length(z.unique)
    layout(matrix(c(1:N.z*2, 1:N.z*2-1), nrow=N.z, ncol=2), widths=c(1, 10), heights=1)
    levels <- 0:N
    for (z.i in z.unique) {
      contours.i <- contours[which(contours[,3] == z.i), ]
      sum.i <- apply(contours.i[, 1:N+3], 1, sum, na.rm=TRUE)
      x <- unique(contours.i[, 1])
      y <- unique(contours.i[, 2])
      lvl.i <- matrix(sum.i, nrow=length(x), ncol=length(y))
      plot(range(x), range(y), type="n", xaxt="n", yaxt="n")
      graphics::.filled.contour(x, y, z=lvl.i, levels=levels, col=c(NA,rev(heat.colors(N-1))))
      contour(x, y, z=lvl.i, levels=levels, col="black", add=TRUE, drawlabels=FALSE, lwd=0.25)
      plot(1,type="n",xaxt="n",yaxt="n")
      text(1, labels=paste("z=", z.i, sep=""))
    }
    par(mar=mar.old)
  }
  return(contours)
}	

compareStructures.surface <- function (structures) {	
  N <- length(structures)
  z <- as.list(rep(NA, N))
  pts <- matrix(nrow=0, ncol=3, dimnames=list(NULL, c("X", "Y", "Z")))
  for (i in 1:N) {
    if (length(structures[[i]]$vertices) < 1) {
      next
    }
    z[[i]] <- unlist(lapply(structures[[i]]$closed.polys, function(closed.poly) {return(unique(closed.poly[,3]))}))
    pts <- rbind(pts, structures[[i]]$vertices)
  }
  results <- matrix(0, nrow=dim(pts)[1], ncol=N, dimnames=list(NULL, names(structures)))
  for (i in 1:N) {
    #		plot3d(structures[[i]]$vertices,col="gray",cex=0.2)
    for (j in unique(z[[i]])) {
      pts.j <- pts[which(pts[, 3]== j), 1:2]
      results.j <- rep(0, dim(pts.j)[1])
      z.j <- which(z[[i]] == j)
      ## THIS LOOP ACCOUNTS FOR AXIAL SLICES WITH MULTIPLE SEPARATE CLOSED POLYGONS (e.g. 3 ROOTS FOR SINGLE TOOTH)
      ## THIS LOOP DOES NOT(!!!) ACCOUNT FOR DONUTS (E.G. STRUCTURES WITH HOLE IN THEM -- NEED TO FIGURE OUT HOW THOSE ARE STORED FIRST) -- IF STRUCTURE HAS A HOLE, ALL BETS ARE OFF AT THE MOMENT... SOLUTION WILL BE TO DO LOGICAL SUBTRACTION RATHER THAN ADDITION OF RESULTS
      for (k in 1:length(z.j)) {
        results.j <- results.j + as.numeric(pointInPoly2D(pts.j[,1:2], structures[[i]]$closed.polys[[z.j[k]]][,1:2]))
      }
      results[which(pts[, 3]== j), i] <- results[which(pts[, 3]== j), i]+results.j
    }
  }
  #	points3d(pts, col=rainbow(n=3)[apply(results,1,sum)])
  #	points3d(pts[which(apply(results,1,sum)==3),],col="black",cex=2)
  return(cbind(pts, results))
}


compareStructures.axial <- function (structures, pixels=100) {	
  N <- length(structures)
  z <- as.list(rep(NA, N))
  bounds <- range(structures, na.rm=TRUE)
  x.coords <- seq(from=bounds[1,1], to=bounds[2,1], length.out=pixels)
  y.coords <- seq(from=bounds[1,2], to=bounds[2,2], length.out=pixels)
  for (i in 1:N) {
    if (length(structures[[i]]$vertices) < 1) {
      next
    }
    z[[i]] <- unlist(lapply(structures[[i]]$closed.polys, function(closed.poly) {return(unique(closed.poly[,3]))}))
  }
  z.coords <- unique(unlist(z))
  pts <- matrix(nrow=length(x.coords)*length(y.coords)*length(z.coords), ncol=3, dimnames=list(NULL, c("X", "Y", "Z")))
  pts <- matrix(c(rep(x.coords, each=length(y.coords)*length(z.coords)), rep(rep(y.coords, each=length(z.coords)), length(x.coords)), rep(z.coords, length(x.coords)*length(y.coords))), nrow=length(x.coords)*length(y.coords)*length(z.coords), ncol=3, dimnames=list(NULL, c("X", "Y", "Z"))) 
  results <- matrix(0, nrow=dim(pts)[1], ncol=N, dimnames=list(NULL, names(structures)))
  for (i in 1:N) {
    for (j in unique(z[[i]])) {
      pts.j <- pts[which(pts[, 3]== j), 1:2]
      results.j <- rep(0, dim(pts.j)[1])
      z.j <- which(z[[i]] == j)
      ## THIS LOOP ACCOUNTS FOR AXIAL SLICES WITH MULTIPLE SEPARATE CLOSED POLYGONS (e.g. 3 ROOTS FOR SINGLE TOOTH)
      ## IF CLOSED POLYGONS ARE NESTED, THEY WILL BE INTERPRETED AS HOLES, SUCH THAT POINTS BETWEEN TWO POLYGONS MAY BE INTERPRETED AS EXTERIOR TO THE POLYGONS THEMSELVES (NOTE THAT THIS ASSUMES THE POLYGONS DO NOT CROSS EACH OTHER AT ANY POINT)
      for (k in 1:length(z.j)) {
        results.j <- results.j + as.numeric(pointInPoly2D(pts.j[,1:2], structures[[i]]$closed.polys[[z.j[k]]][,1:2]))
      }
      results[which(pts[, 3]== j), i] <- results[which(pts[, 3]== j), i]+(results.j %% 2 != 0)	
    }
  }
  return(cbind(pts, results))
}


compareStructures.wDSC <- function (structures, dose, pixels = 100){
  print("Testing if Voxels Contain Any Structures")
  N <- length(structures) # how many structures are there
  bounds <- range(structures, na.rm=TRUE) 
  voxelDimensions <- c((bounds[2,1]-bounds[1,1])/(pixels-1), (bounds[2,2]-bounds[1,2])/(pixels-1), (bounds[2,3]-bounds[1,3])/(pixels-1))
  halfVoxels <- voxelDimensions*1/2 #defines distance from midpoint where dose will be approximated
  x.coords <- seq(from=bounds[1,1], to=bounds[2,1], length.out=pixels) 
  y.coords <- seq(from=bounds[1,2], to=bounds[2,2], length.out=pixels) 
  toSelectFrom <- compareStructures.axial(structures, pixels)
  columns <- length(structures) + 3 #allows function to work with any number of structures
  toFill <- matrix(rep(NA, times = columns*nrow(toSelectFrom)), ncol= columns) #this will be filled with rows from toSelectFrom that are in at least one structure, saves a lot of computation time
  print("Selecting Voxels Containing Structures")
  for (p in 1:nrow(toSelectFrom)) {
    iszero <- c(0)
    for (u in 4:columns) {
      iszero[u] <- toSelectFrom[p,u] == 0
    }
    if(sum(iszero, na.rm = TRUE) == length(structures)){ #iszero = true when row contains coordintes not in any structure, so we leave this row in toFill with NAs 
      next()
    }else{
      for (s in 1:columns) {
        toFill[[p,s]] <- toSelectFrom[[p,s]] #if iszero >0, then the row contains a coordinate in at least one structure and needs to be used in calculating the wDSC
      }
    }
  }
  newMat <- toFill[apply(toFill, 1, function(x)!any(is.na(x))), , drop=F] #removes rows with NAs (that weren't in any structure)
  xStep <- halfVoxels[1]  # distance from voxel midpoint where dose will be aproximated 
  yStep <- halfVoxels[2]
  zStep <- halfVoxels[3]
  magnitudesMatrix <- matrix(0, ncol = 1, nrow = nrow(newMat)) #this will be be the magnitudes of the net radiation dose gradient vector 
  dhms <- function(t){  #this function will estimate time that the calculation will take
    paste(t %/% (60*60*24) 
          ,paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0")
                 ,formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
                 ,formatC(t %% 60, width = 2, format = "d", flag = "0")
                 ,sep = ":"
          )
    )
  }
  ETAseconds <- nrow(newMat)/150
  print("Calculating Radiation Dose Gradient At All Voxels")
  print("Estimated time for full weighted DSC calculation in D H:M:S") 
  print(dhms(ETAseconds))
  pb <- txtProgressBar(min = 0, max = nrow(newMat), style = 3)
  for (q in 1:nrow(newMat)) {
    setTxtProgressBar(pb, q)
    midpointX <- newMat[q,1]   #this is the x, y, z coordinate of the voxel center, we will move in halfVoxels steps in all 14 directions
    midpointY <- newMat[q,2]
    midpointZ <- newMat[q,3]
    doseAtMidpoint <- approx3D(data = dose, x = midpointX , y = midpointY, z = midpointZ)
    #the below vectorMatrix is a matrix that has columns i, j, k, deltaDose. The first three represent the direction of the vector and the deltaDose is the difference in radiation dose at the end of that vector and the midpoint of the voxel
    vectorMatrix <- matrix(data = c(0,0,0,0,-1,1,-1,-1,-1,-1,1,1,1,1,  #i column
                                    0,0,1,-1,0,0,1,-1,1,-1,1,-1,1,-1,  #j column
                                    1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,  #k column
                                    approx3D(data = dose, x = midpointX, y = midpointY, z = midpointZ + zStep) - doseAtMidpoint, #up
                                    approx3D(data = dose, x = midpointX, y = midpointY, z = midpointZ - zStep) - doseAtMidpoint, #down
                                    approx3D(data = dose, x = midpointX, y = midpointY + yStep, z = midpointZ) - doseAtMidpoint, #right
                                    approx3D(data = dose, x = midpointX, y = midpointY - yStep, z = midpointZ) - doseAtMidpoint, #left
                                    approx3D(data = dose, x = midpointX - xStep, y = midpointY, z = midpointZ) - doseAtMidpoint, #forward
                                    approx3D(data = dose, x = midpointX + xStep, y = midpointY, z = midpointZ) - doseAtMidpoint, #backward
                                    approx3D(data = dose, x = midpointX - xStep, y = midpointY + yStep, z = midpointZ + zStep) - doseAtMidpoint, #forward up right
                                    approx3D(data = dose, x = midpointX - xStep, y = midpointY - yStep, z = midpointZ + zStep) - doseAtMidpoint, #forward up left
                                    approx3D(data = dose, x = midpointX - xStep, y = midpointY + yStep, z = midpointZ - zStep) - doseAtMidpoint, #forward down right
                                    approx3D(data = dose, x = midpointX - xStep, y = midpointY - yStep, z = midpointZ - zStep) - doseAtMidpoint, #forward down left
                                    approx3D(data = dose, x = midpointX + xStep, y = midpointY + yStep, z = midpointZ + zStep) - doseAtMidpoint, #backward up right
                                    approx3D(data = dose, x = midpointX + xStep, y = midpointY - yStep, z = midpointZ + zStep) - doseAtMidpoint, #backward up left
                                    approx3D(data = dose, x = midpointX + xStep, y = midpointY + yStep, z = midpointZ - zStep) - doseAtMidpoint, #backward down right
                                    approx3D(data = dose, x = midpointX + xStep, y = midpointY - yStep, z = midpointZ - zStep) - doseAtMidpoint #backward down left
    ),
    
    ncol = 4)
    for (w in 1:14) { #in this loop, the i, j, and k columns are multiplied by the deltaDose 
      vectorMatrix[w, 1] <- vectorMatrix[w, 1] * vectorMatrix[w, 4]
      vectorMatrix[w, 2] <- vectorMatrix[w, 2] * vectorMatrix[w, 4]
      vectorMatrix[w, 3] <- vectorMatrix[w, 3] * vectorMatrix[w, 4]
    }
    netVector <- colSums(vectorMatrix[,c(1,2,3)]) #the columns are summed to get the net vector, the mangnitude of the vector is calculated on the next line
    magnitudesMatrix[q,] <- sqrt(netVector[1]^2 + netVector[2]^2 + netVector[3]^2)
  }
  close(pb)
  multipliedMatrix <- newMat[ , c(4:columns)] * c(magnitudesMatrix) #takes the matrix with 1 and 0 representing whether or not a voxel belongs to a structure and multiplies the 1 by the dose gradient magnitude
  DSCtable <- matrix(0, nrow = columns - 3, ncol = columns - 3, dimnames=list(names(structures), names(structures)))
  for (i in 1:(columns-3)) {
    for (j in 1:(columns-3)) {
      if (i == j) {
        DSCtable[i, j] <- 1      #DSC for identical structures is 1
        next()
      }
      alpha <- colSums(multipliedMatrix)[i] #represents |A| sum all gradient magnitudes for structure A
      beta <-  colSums(multipliedMatrix)[j] #represents |B| sum all gradient magnitudes for structure B
      gamma <- matrix(0, ncol = 1, nrow = nrow(multipliedMatrix))
      for (k in 1:nrow(multipliedMatrix)) {
        if(multipliedMatrix[k,i] == multipliedMatrix[k,j]){
          gamma[k] <- multipliedMatrix[k,i] #if voxel is in both structures, add the gradient magnitude to gamma
        }else{
          gamma[k] <- 0
        }
      }
      gamma <- sum(gamma, na.rm = TRUE) #represnts |A ∩ B|
      DSCtable[i, j] <- 2*gamma/(alpha+beta)
    }
  }
  print(DSCtable)
}

compareStructures.hausdorff <- function (structures, verbose=TRUE, method=NULL) {
  
  hausdorff.dist <- function (A, B, method) {
    if (ncol(A) != ncol(B)){
      warning("Dimensionality of A and B must be the same")
      return(NA)
    }
    compute.dist = function (a0, B0){
      C0 <- matrix(rep(a0, each=nrow(B0)),byrow=F, ncol=ncol(B0))
      return(min(apply(C0-B0, 1, function(x) {sqrt(sum(t(x)*x))}), na.rm=TRUE))
    }
    
    if (method == "mean") {
      d1 <- apply(A, 1, compute.dist, B0=B)
      d2 <- apply(B, 1, compute.dist, B0=A)
      return(mean(c(d1, d2), na.rm=TRUE))
    }
    else if (method == "median") {
      d1 <- apply(A, 1, compute.dist, B0=B)
      d2 <- apply(B, 1, compute.dist, B0=A)
      return(median(c(d1, d2), na.rm=TRUE))
    }
    else if (method == "absolute") {
      d1 <- max(apply(A, 1, compute.dist, B0=B))
      d2 <- max(apply(B, 1, compute.dist, B0=A))
      return(max(d1, d2, na.rm=TRUE))
    }
    else {
      warning("Invalid 'method' argument; must be one of 'mean', 'median', or 'absolute'")
      return(NA)
    }
  }
  
  method <- match.arg(method, choices=c("mean", "median", "absolute"))
  N <- length(structures)
  results <- matrix(0, nrow=N, ncol=N, dimnames=list(names(structures), names(structures)))
  for (i in 1:N) {
    if (verbose) {
      cat("Analyzing structure ", i, "/", N, " (", structures[[i]]$name, ") ... ", sep="")
    }
    for (j in i:N) {
      results[i, j] <- hausdorff.dist(structures[[i]]$vertices, structures[[j]]$vertices, method=method)
      results[j, i] <- results[i, j]
    }
    if (verbose) {
      cat("FINISHED\n")
    }
  }
  return(results)
}

compareStructures.EMD <- function (structures) {	
  cat("EMD method currently unavailable (in development) -- will be available in a future release\n")
  return()
}

pointInPoly2D <- function (points, poly) {
  poly <- matrix(unique(poly), ncol=2)
  n <- dim(poly)[1]
  x <- diff(poly[c(1:n,1),1])
  y <- poly[,2] + poly[c(2:n,1),2]
  if (sum(x*y/2) >= 0) {
    #clockwise poly
    return(pip2d(poly[n:1,], points) >= 0)
  }
  else {
    #anti-clockwise poly
    return(pip2d(poly, points) >= 0)
  }
}
