plot2DVelocityFieldFromL = function(a, b, c, d, gridSpacing = 2) {
     
     #Create a grid of locations
     locations = createPointGrid(gridSpacing = gridSpacing)
     
     L = rbind(c(a, b),
               c(c, d))
     
     velocityDataFrame = velocitiesAndLocationsFromLocationsAndL(L, locations)
     
     plotField(velocityDataFrame, end = FALSE)

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


displacementFieldIncrement <- function(posGradTens, gridSpacing) {
     

     locs = createPointGrid(gridSpacing = gridSpacing)
     
     positionGradientTensor = posGradTens
     
     displacementsDataFrame = displacementsFromLocationsAndF(positionGradientTensor, locs)
     
     return(displacementsDataFrame)
}

strainEllipseIncrement = function(positionGradientTensor, radius = 8, numPoints = 32) {
     locs = createPointCircle(radius = radius, numPoints = numPoints)
     
     finiteStrain = computeFiniteStrainSemiAxes(positionGradientTensor, radius)
     
     maxStrain = data.frame(t(finiteStrain[[1]]))
     names(maxStrain) = c("xStart", "yStart")
     maxStrain$xEnd = -maxStrain$xStart
     maxStrain$yEnd = -maxStrain$yStart
     
     
     minStrain = data.frame(t(finiteStrain[[2]]))
     names(minStrain) = c("xStart", "yStart")
     minStrain$xEnd = -minStrain$xStart
     minStrain$yEnd = -minStrain$yStart
     
     maxStretch = finiteStrain[[3]]
     minStretch = finiteStrain[[4]]
     
     displacementsDataFrame = displacementsFromLocationsAndF(positionGradientTensor, locs)
     
     
     return(list(displacementsDataFrame, maxStrain, minStrain, maxStretch, minStretch))
     
}

materialLineIncrement <- function(positionGradientTensor, thetaStart, radius = 15) {
     startPos = list(c(cos(thetaStart) * radius ,sin(thetaStart) * radius))
     endPos = displacementsFromLocationsAndF(positionGradientTensor, startPos)
     stretch = (endPos$xEnd^2 + endPos$yEnd^2)^(1/2) / (endPos$xStart^2 + endPos$yStart^2)^(1/2)
     endPos$stretch = stretch
     return(endPos)
}




plotParticlePaths <- function(velocityGradientTensor, totalTime, numberOfIncrements = 100, gridSpacing = 2.5, superimposeStrain = FALSE, superimposeMaterialLine = FALSE, materialLineAngle = 60) {
     
     displacementsDataFrame = displacementFieldIncrement(rbind(c(1,0),c(0,1)), gridSpacing = gridSpacing)
     
     if (superimposeStrain == TRUE) {
         strain = strainEllipseIncrement(rbind(c(1,0),c(0,1)), radius = 10, numPoints = 100)
         strainIncrementDataFrame = strain[[1]]
         strainMaxDataFrame = strain[[2]]
         strainMinDataFrame = strain[[3]]
              
          strainIncrementDataFrame$frame = 0
          strainMaxDataFrame$frame = 0
          strainMaxDataFrame$maxStretch = strain[[4]]
          strainMinDataFrame$frame = 0
          strainMinDataFrame$minStretch = strain[[5]]
     
     }
     
     
     if (superimposeMaterialLine == TRUE) {
          
          theta = materialLineAngle*degree
          materialLineDataFrame = materialLineIncrement(rbind(c(1,0),c(0,1)), theta)
          materialLineDataFrame$frame = 0
     }
     
     displacementsDataFrame$frame = 0
     
     exDot = velocityGradientTensor[1,1]
     gammaDot = velocityGradientTensor[1,2]
     c = velocityGradientTensor[2,1]
     xyDot = velocityGradientTensor[2,2]
     
     
     for (i in 1:numberOfIncrements) {
                    t = totalTime/numberOfIncrements * i
                    if (exDot == 0) {
                         
                         posGradTens = rbind(c(exp(exDot*t), gammaDot * t          ),
                                             c(c          , exp(eyDot*t)               )) 
                         
                    } else {
                         
                         posGradTens = rbind(c(exp(exDot*t), gammaDot/exDot*sinh(exDot*t)),
                                             c(c           , exp(eyDot*t)               ))
                         
                    }
                    
                    increment = displacementFieldIncrement(posGradTens, gridSpacing = gridSpacing)
                    
                    increment$frame = i
                    
                    displacementsDataFrame = rbind(displacementsDataFrame, increment)
                    
                    if (superimposeMaterialLine == TRUE) {
                         materialIncrement = materialLineIncrement(posGradTens, theta)
                         materialIncrement$frame = i
                         
                         materialLineDataFrame = rbind(materialLineDataFrame, materialIncrement)
                    }
                    
                    if (superimposeStrain == TRUE) {
                         strain = strainEllipseIncrement(posGradTens, radius = 10, numPoints = 100)
                         
                         strainIncrement = strain[[1]]
                         strainMax = strain[[2]]
                         strainMin = strain[[3]]
                              
                         strainIncrement$frame = i
                         strainMax$frame = i
                         strainMax$maxStretch = strain[[4]]
                         
                         strainMin$frame = i
                         strainMin$minStretch = strain[[5]]
                         
                         strainIncrementDataFrame = rbind(strainIncrementDataFrame,strainIncrement)
                         strainMaxDataFrame = rbind(strainMaxDataFrame, strainMax)
                         strainMinDataFrame = rbind(strainMinDataFrame, strainMin)
                    }
                    
                    
     }
     plot <- plotAnimatedField(displacementsDataFrame, arrow = FALSE, start = FALSE)
     
     # if (superimposeMaterialLine == TRUE) {
     #      plot = plot +
     #           geom_segment(data = materialLineDataFrame, aes(x = -xEnd, y = -yEnd, xend = xEnd, yend = yEnd), color = "brown")
     #           
     # }
     
     if (superimposeStrain == TRUE && superimposeMaterialLine == TRUE) {
          plot = plot +
               geom_point(data = strainIncrementDataFrame, aes(x = xEnd, y = yEnd), color = "black", size = 1, alpha = 1) +
               geom_segment(data = strainMaxDataFrame, aes(x=xStart, y=yStart,xend=xEnd,yend=yEnd), color = "blue", size = 2) +
               geom_segment(data = strainMinDataFrame, aes(x=xStart, y=yStart,xend=xEnd,yend=yEnd), color = "orange", size = 2) +
               geom_segment(data = materialLineDataFrame, aes(x = -xEnd, y = -yEnd, xend = xEnd, yend = yEnd), color = "brown", size = 2) +
               shadow_mark(alpha = 0.3, size = 0.25, color = 'grey', exclude_layer = c(5,6,7))
          
           plot2 = ggplot() +
                geom_hline(yintercept = 1) +
                geom_point(data = strainMaxDataFrame, aes(x = frame, y = maxStretch), color = "blue", size = 3) +
                geom_point(data = strainMinDataFrame, aes(x = frame, y = minStretch), color = "orange", size = 3) +
                geom_point(data = materialLineDataFrame, aes(x = frame, y = stretch), color = "brown", size = 3) +
                theme_classic() +
                transition_time(frame) +
                shadow_mark(alpha = 0.3, size = 1) +
                ease_aes('linear')
         
           plot = list(plot, plot2) 
     }
     
     if (superimposeStrain == TRUE && superimposeMaterialLine == FALSE) {
          plot = plot +
               geom_point(data = strainIncrementDataFrame, aes(x = xEnd, y = yEnd), color = "black", size = 1, alpha = 1) +
               geom_segment(data = strainMaxDataFrame, aes(x=xStart, y=yStart,xend=xEnd,yend=yEnd), color = "blue", size = 2) +
               geom_segment(data = strainMinDataFrame, aes(x=xStart, y=yStart,xend=xEnd,yend=yEnd), color = "orange", size = 2) +
               shadow_mark(alpha = 0.3, size = 0.25, color = 'grey', exclude_layer = c(5,6))
          
          plot2 = ggplot() +
               geom_hline(yintercept = 1) +
               geom_point(data = strainMaxDataFrame, aes(x = frame, y = maxStretch), color = "blue", size = 3) +
               geom_point(data = strainMinDataFrame, aes(x = frame, y = minStretch), color = "orange", size = 3) +
               theme_classic() +
               transition_time(frame) +
               shadow_mark(alpha = 0.3, size = 1) +
               ease_aes('linear')
          
          plot = list(plot, plot2) 
     }
     
     if (superimposeStrain == FALSE && superimposeMaterialLine == TRUE) {
          plot = plot +
               geom_segment(data = materialLineDataFrame, aes(x = -xEnd, y = -yEnd, xend = xEnd, yend = yEnd), color = "brown", size = 2) +
               shadow_mark(alpha = 0.3, size = 0.25, color = 'grey', exclude_layer = c(5))
          
          plot2 = ggplot() +
               geom_hline(yintercept = 1) +
               geom_point(data = materialLineDataFrame, aes(x = frame, y = stretch), color = "brown", size = 3) +
               theme_classic() +
               transition_time(frame) +
               shadow_mark(alpha = 0.3, size = 1) +
               ease_aes('linear')
          
          plot = list(plot, plot2) 
     }
     
     return(plot)
     
     
     
}


plotAnimatedField = function(dataFrameWithStartsEnds, start = TRUE, arrow = TRUE, end = TRUE) {
     p =  ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          theme_void()
     
     if (start == TRUE && end == TRUE && arrow == TRUE) {p = p + 
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
     
     } 
     if (start == TRUE && end == TRUE && arrow == FALSE) {p = p + 
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75)
     } 
     
     if (start == TRUE && end == FALSE && arrow == FALSE) {p = p + 
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) 
     } 
     
     if (start == FALSE && end == TRUE && arrow == FALSE) {p = p + 
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "black", size = 1, alpha = 1) 
     
     } 
     if (start == FALSE && end == TRUE && arrow == TRUE) {p = p + 
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
     
     } 
     
     if (start == FALSE && end == FALSE && arrow == TRUE) {p = p + 
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
     
     } 
     
     if (start == TRUE && end == FALSE && arrow == TRUE) {p = p + 
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
     
     } 
     
     
     p = p + 
          transition_time(frame) +
          shadow_mark(alpha = 0.3, size = 0.25, color = 'grey')+
          ease_aes('linear')
     return(p)
}

# field headings in the data frame should have xStart, yStart, xEnd, yEnd whether the field is velocity or position gradients
plotField = function(dataFrameWithStartsEnds, start = TRUE, arrow = TRUE, end = TRUE) {
    p =  ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          theme_void()
     
     if (start == TRUE && end == TRUE && arrow == TRUE) {p = p + 
                                             geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
                                             geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
                                             geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
     
     } 
    if (start == TRUE && end == TRUE && arrow == FALSE) {p = p + 
         geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
         geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75)
    } 
    
    if (start == TRUE && end == FALSE && arrow == FALSE) {p = p + 
         geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) 
    } 
        
    if (start == FALSE && end == TRUE && arrow == FALSE) {p = p + 
         geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "black", size = 0.25, alpha = 1) 

    } 
    if (start == FALSE && end == TRUE && arrow == TRUE) {p = p + 
         geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
         geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
    
    } 
    
    if (start == FALSE && end == FALSE && arrow == TRUE) {p = p + 
         geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
    
    } 
    
    if (start == TRUE && end == FALSE && arrow == TRUE) {p = p + 
         geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
         geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue")
    
    } 
    return(p)
}

plotFieldWithMaxStrain = function(dataFrameWithStartsEnds, maxS, minS) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          geom_segment(data = maxS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 2) +
          geom_segment(data = minS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 2) +
          theme_void()
}


plotFieldWithOffsetMarker = function(dataFrameWithStartsEnds, offsetMarker) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          geom_segment(data = offsetMarker, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "black", size = 2) +
          theme_void()
}


plotFieldWithOffsetAndMaxStrain = function(dataFrameWithStartsEnds, maxS, minS, offsetMarker) {
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart), color = "black", size = 2) +
          geom_point(data = dataFrameWithStartsEnds, aes(x = xEnd, y = yEnd), color = "grey", size = 2, alpha = 0.75) +
          geom_segment(data = dataFrameWithStartsEnds, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          geom_segment(data = maxS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 2) +
          geom_segment(data = minS, aes(x = -X1, y = -X2, xend = X1, yend=X2), color = "grey", size = 2) +
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

createPointCircle = function(radius, numPoints = 32) {
     locations = list()
     counter = 1
     for(i in seq(0,2*pi, 2*pi/numPoints)){
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
     
     identityMat = rbind(c(1,0),
                         c(0,1))
     
     if(all(positionGradientTensor == identityMat)){
          maxStrain = c(0,0)
          minStrain = c(0,0)
     } else{ 

     
     maxStrain = c(maxValue * radius * maxVector[1], maxValue * radius * maxVector[2])
     minStrain = c(minValue * radius * minVector[1], minValue * radius * minVector[2])
}
     
     return(list(maxStrain, minStrain, maxValue, minValue))

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

computePrincipalLFromF = function(positionGradientTensor) {
     princL = logm(positionGradientTensor) 
     return(princL)
}

