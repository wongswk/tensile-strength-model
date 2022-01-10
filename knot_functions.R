# FUNCTION1
# function to decide if a point is inside an ellipse or not
is_inside_ellipse = function(point, ell_x, ell_y, var_x, var_y, cov){
  scale = 1
  
  indicator = FALSE
  x = point[1]
  y = point[2]
  
  # covariance matrix
  A = scale*matrix(c(var_x, cov, cov, var_y),byrow = T,2,2)
  eigenA = eigen(A)
  major = 2*sqrt(max(eigenA$values)) # major axis
  minor = 2*sqrt(min(eigenA$values)) # minor axis
  rotation = diag(eigenA$vectors)[1] # in radian 
  cos_angle = cos((180-rotation)*pi/180)
  sin_angle = sin((180-rotation)*pi/180)
  xc = x-ell_x
  yc = y-ell_y
  xct = xc*cos_angle - yc*sin_angle
  yct = xc*sin_angle + yc*cos_angle
  rad_cc = (xct^2/(minor/2)^2) + (yct^2/(major/2)^2)
  if(rad_cc<=1){
    indicator <- TRUE
  } else {indicator <- FALSE}
  return(indicator)
}


#### FUNCTION 2 ####
# sample points within given ellipse parameters
sample_points_in_ellipse = function(ell_x, ell_y, var_x, var_y, cov, num_points = 10){
  library(MASS)
  i = 1
  points = data.frame(matrix(ncol = 2, nrow = 0))
  while(i<= num_points){
    # randomly pick one point
    point = mvrnorm(n = 1, mu = c(ell_x, ell_y), Sigma = matrix(c(var_x, cov, cov, var_y),byrow = T,2,2))
    if(is_inside_ellipse(point = point, ell_x = ell_x, ell_y = ell_y, var_x = var_x, var_y = var_y, cov = cov)){
      i = i+1
      points = rbind(points, point)
    }
  }
  colnames(points) = c('x','y')
  return(points)
}


#### FUNCTION 3 ####
# get the matched knot_list 

get_matched_knot_list = function(board_dat,
                                 x_start,
                                 x_end){
  library(rlist)
  matching = unique(board_dat$matching)
  matched_knot_list = list()
  #browser()
  # check if the lumber has knots
  if(length(matching)!=0){
    for(i in matching){
      matched_knot = subset(board_dat, matching == i)
      # check if area_over>0 -- edge knots
      if(any(matched_knot$area_over>0)){
        matched_knot$edge_knot = 1
      } else {
        matched_knot$edge_knot = 0
        }
      # check if any one of the matched knots lie within the test span
      # if so, add to the matched knot list
      if(all(matched_knot$x >= x_start & matched_knot$x <= x_end))
        matched_knot_list = list.append(matched_knot_list,matched_knot)
    }
  }

  return(matched_knot_list)
}


#### FUNCTION 4 ####
# get the valid points that lies within ellipse and 3D lumber 
# NOTE the meaning of the coordinates
get_valid_points = function(matched_knots,x_start, x_end,
                            y_start, y_end,
                            z_start, z_end){
  
  library(tidyverse)
  scale = 1
  valid_points = data.frame(matrix(nrow = 0, ncol = 3))
  names(valid_points) = c('x','y','z')
  # browser()
  for (i in 1:nrow(matched_knots)){
    dat_row = matched_knots[i,]
    surface = dat_row$surface # surface
    x = dat_row$x
    y = dat_row$y
    var_x = dat_row$var_x
    var_y = dat_row$var_y
    cov = dat_row$cov
    Sigma = matrix(c(var_x, cov, cov, var_y), 2,2, byrow = T)
    RR = chol(Sigma)
    angle_sample <- seq(0, 2*pi, length.out=50)          
    ell  <- sqrt(scale) * cbind(cos(angle_sample), sin(angle_sample)) %*% RR 
    points_on <- as.data.frame(sweep(ell, 2, c(x,y), "+"))
    colnames(points_on) <- c('x','y')
    # points_on$z = dat_row$z
    points_in = sample_points_in_ellipse(ell_x = x, ell_y = y,
                                         var_x = var_x, var_y = var_y,
                                         cov = cov, num_points = 100)
    # points_in$z =  dat_row$z
    temp = rbind(points_on, points_in)
    # NOTE that the above x y are with respect to the current surface
    # Use the surface information to project the points onto 3D coordinates
    if(surface == 0){
      points_all= data.frame(x = temp$x, 
                             y = y_end - temp$y, 
                             z = z_start)
      colnames(points_all) = c('x', 'y', 'z')
    }
    

    if(surface==1){
      points_all = data.frame(x = temp$x, 
                              y = y_end, 
                              z = z_end - temp$y)
                          
      colnames(points_all) = c('x', 'y', 'z')
    }
    
    if(surface==2){
      points_all= data.frame(x = temp$x, 
                             y = temp$y, 
                             z = z_end)
      colnames(points_all) = c('x', 'y', 'z')
    }
    
    if(surface==3){
      points_all = data.frame(x = temp$x, 
                              y = y_start, 
                              z = temp$y)
      colnames(points_all) = c('x', 'y', 'z')
    }
    
    points_all %>% filter(x >= x_start & x <= x_end 
                          & y >= y_start & y <= y_end 
                          & z >= z_start & z <= z_end)
    
    valid_points = rbind(valid_points, points_all)
  }
  return(valid_points)
}


#### FUNCTION 5 ####
# get the matched 3D knots volumes
# matched_knot_list can be empty
get_knots_volume = function(matched_knot_list,
                            x_start, x_end,
                            y_start, y_end,
                            z_start, z_end
){
  library(geometry)
  knots_volume = data.frame(matrix(nrow = 0, ncol = 4))
  # browser()
  if(length(matched_knot_list)!=0){
  for (i in 1:length(matched_knot_list)) {
    matched_knots = matched_knot_list[[i]]
    valid_points = get_valid_points(matched_knots = matched_knots,  
                                    x_start, x_end,
                                    y_start, y_end,
                                    z_start, z_end)
    if(nrow(valid_points)!= 0){
      convexHull = convhulln(p = as.matrix(valid_points), output.options = 'FA') 
      knots_volume = rbind(knots_volume, c(i, convexHull$vol, mean(matched_knots$x), mean(matched_knots$y)))
      
    } 
  }
    }
  names(knots_volume) = c('matching', 'volume', 'x_coord', 'y_coord')
  return(knots_volume)
}