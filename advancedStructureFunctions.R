plot2DVelocityFieldFromL = function(a, b, c, d, gridSpacing = 2) {
     
     #Create a grid of locations
     locations = createVelocityGrid(gridSpacing = gridSpacing)
     
     L = rbind(c(a, b),
               c(c, d))
     
     velocityDataFrame = velocitiesAndLocationsFromLocationsAndL(L, locations)
     
     ggplot() +
          coord_fixed(ratio = 1) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          geom_segment(data = velocityDataFrame, aes(x = xStart, y = yStart, xend = xEnd, yend=yEnd), arrow = arrow(length=unit(3,"pt"), ends="last", type = "closed"), color = "blue") +
          #geom_point(data = velocityDataFrame, aes(x = xStart, y = yStart), size = 1,color = "blue") +
          theme_void()

}


# Helper functions
velocitiesAndLocationsFromLocationsAndL = function(L, locations) {
     velocities = t(sapply(locations, function(s) L %*% s))
     locations = t(sapply(locations,function(s) s))
     locationsAndVelocities = data.frame(cbind(locations,velocities))
     names(locationsAndVelocities) = c("xStart", "yStart", "xVel", "yVel")
     locationsAndVelocities$xEnd = locationsAndVelocities$xStart + locationsAndVelocities$xVel
     locationsAndVelocities$yEnd = locationsAndVelocities$yStart + locationsAndVelocities$yVel
     return(locationsAndVelocities)
}

createVelocityGrid = function(xmin = -10, ymin = -10, xmax = 10, ymax = 10, gridSpacing = 0.5){
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

