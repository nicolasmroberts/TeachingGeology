plot2DVelocityFieldFromL = function(a, b, c, d, gridSpacing = 2) {
     
     #Create a grid of locations
     locations = createPointGrid(gridSpacing = gridSpacing)
     
     L = rbind(c(a, b),
               c(c, d))
     
     velocityDataFrame = velocitiesAndLocationsFromLocationsAndL(L, locations)
     
     plotField(velocityDataFrame)

}


plot2DDisplacementFieldFromF = function(a, b, c, d, gridSpacing = 2, radius = 8, plotPsi = FALSE, plotStrain = FALSE) {
     locs = createPointGrid(gridSpacing = gridSpacing)
     
     positionGradientTensor = rbind(c(a, b),
                                   c(c, d))
     
     displacementsDataFrame = displacementsFromLocationsAndF(positionGradientTensor, locs)
     
     
     
     offsetMarker1 = computeVerticalOffset(positionGradientTensor)
     
     
     
     finiteStrain = computeFiniteStrainSemiAxes(positionGradientTensor, radius)
     
     maxStrain = data.frame(t(finiteStrain[[1]]))
     minStrain = data.frame(t(finiteStrain[[2]]))
     
     if (plotPsi == TRUE && plotStrain == TRUE){
          plotFieldWithOffsetAndMaxStrain(displacementsDataFrame, offsetMarker = offsetMarker1, maxS = maxStrain, minS = minStrain)
     } else if (plotPsi == TRUE && plotStrain == FALSE) {
          plotFieldWithOffsetMarker(displacementsDataFrame, offsetMarker1)
     } else if (plotPsi == FALSE && plotStrain == TRUE) {
          plotFieldWithMaxStrain(displacementsDataFrame, maxS = maxStrain, minS = minStrain)
     } else if(plotPsi == FALSE && plotStrain == FALSE) {
          plotField(displacementsDataFrame)
     }
}


plot2DCircleDisplacementFieldFromF = function(a,b,c,d, radius = 8, plotPsi = FALSE, plotStrain = FALSE) {
     locs = createPointCircle(radius = radius)
     
     positionGradientTensor = rbind(c(a, b),
                                    c(c, d))
     
     offsetMarker1 = computeVerticalOffset(positionGradientTensor)
     
     
     finiteStrain = computeFiniteStrainSemiAxes(positionGradientTensor, radius)
     
     maxStrain = data.frame(t(finiteStrain[[1]]))
     minStrain = data.frame(t(finiteStrain[[2]]))
     
     displacementsDataFrame = displacementsFromLocationsAndF(positionGradientTensor, locs)
     
     if (plotPsi == TRUE && plotStrain == TRUE){
          plotFieldWithOffsetAndMaxStrain(displacementsDataFrame, offsetMarker =  offsetMarker1, maxS = maxStrain, minS = minStrain)
     } else if (plotPsi == TRUE && plotStrain == FALSE) {
          plotFieldWithOffsetMarker(displacementsDataFrame, offsetMarker1)
     } else if (plotPsi == FALSE && plotStrain == TRUE) {
          plotFieldWithMaxStrain(displacementsDataFrame, maxS = maxStrain, minS = minStrain)
     } else if(plotPsi == FALSE && plotStrain == FALSE) {
          plotField(displacementsDataFrame)
     }
 
     #plotFieldWithMaxStrain(displacementsDataFrame, maxStrainSlope)
     
         
}

# field headings in the data frame should have xStart, yStart, xEnd, yEnd whether the field is velocity or position gradients
plotField = function(dataFrameWithStartsEnds) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 3) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 3, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          theme_void()
}

plotFieldWithMaxStrain = function(dataFrameWithStartsEnds, maxS, minS) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 3) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 3, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          geom_segment(data = maxS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 3) +
          geom_segment(data = minS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 3) +
          theme_void()
}


plotFieldWithOffsetMarker = function(dataFrameWithStartsEnds, offsetMarker) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 3) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 3, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          geom_segment(data = offsetMarker, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "black", size = 2) +
          theme_void()
}


plotFieldWithOffsetAndMaxStrain = function(dataFrameWithStartsEnds, maxS, minS, offsetMarker) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 3) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 3, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          geom_segment(data = maxS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 3) +
          geom_segment(data = minS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 3) +
          geom_segment(data = offsetMarker, aes(x = -X1, y = -X2, xend = X1, yend = X2), color = "black", size = 2) +
          theme_void()
}


# Creates a grid of points that you can use to draw velocity or displacement fields
createPointGrid = function(xmin = -10, ymin = -10, xmax = 10, ymax = 10, gridSpacing = 0.5){
     locations = list()
     counter = 1
     for(i in seq(xmin, xmax, gridSpacing)){
          for(n in seq(ymin,ymax, gridSpacing)){
               locations[[counter]] <- c(i, n)
               counter = counter + 1
          }
     }
     return(locations)
}

createPointCircle = function(radius) {
     locations = list()
     counter = 1
     for(i in seq(0,2*pi, 2*pi/32)){
          locations[[counter]] = c(radius * cos(i), radius * sin(i)) 
          counter = counter + 1
     }
     return(locations)
}

# Creates a data frame from a list of locations and the velocity gradient tensor
velocitiesAndLocationsFromLocationsAndL = function(L, locations) {
     velocities = t(sapply(locations, function(s) L %*% s))
     locations = t(sapply(locations,function(s) s))
     locationsAndVelocities = data.frame(cbind(locations,velocities))
     names(locationsAndVelocities) = c("xStart", "yStart", "xVel", "yVel")
     locationsAndVelocities$xEnd = locationsAndVelocities$xStart + locationsAndVelocities$xVel
     locationsAndVelocities$yEnd = locationsAndVelocities$yStart + locationsAndVelocities$yVel
     return(locationsAndVelocities)
}


# Creates a data frame from a list of locations and the position gradient tensor
displacementsFromLocationsAndF = function(positionGradientTensor, locations) {
     endPositions = t(sapply(locations, function(S) positionGradientTensor %*% S))
     locations = t(sapply(locations, function(S) S))
     locationsAndEndPositions = data.frame(cbind(locations, endPositions))
     names(locationsAndEndPositions) = c("xStart", "yStart", "xEnd", "yEnd")
     return(locationsAndEndPositions)
}


computeFiniteStrainSemiAxes = function(positionGradientTensor, radius) {     
     
     finiteStrainEigens = computeFiniteStrainAxes(positionGradientTensor)
     maxVector = finiteStrainEigens$vectors[,1]
     maxValue = finiteStrainEigens$values[1]^(1/2)
     
     minVector = finiteStrainEigens$vectors[,2]
     minValue = finiteStrainEigens$values[2]^(1/2)
     
     maxStrain = c(maxValue * radius * maxVector[1], maxValue * radius * maxVector[2])
     minStrain = c(minValue * radius * minVector[1], minValue * radius * minVector[2])

     return(list(maxStrain, minStrain))

}

computeFiniteStrainAxes = function(positionGradientTensor) {
     finiteStrain = positionGradientTensor %*% t(positionGradientTensor)
     finiteStrainEigens = eigen(finiteStrain)
     return(finiteStrainEigens)
}

computeVerticalOffset = function(positionGradientTensor) {
     x = 0
     y = 10
     xy = rbind(c(x), c(y))
     XYend = positionGradientTensor %*% xy
     Xend = XYend[1,1]
     Yend = XYend[2,1]
     return(data.frame(t(c(Xend, Yend))))
}

