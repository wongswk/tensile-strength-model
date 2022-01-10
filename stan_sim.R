## Helper functions for simulating data from model

rm(list = ls())
library(mvtnorm)
library(spatstat)
library(raster)
library(geostatsp)

# simulate the number and location of knots.
get_PPlist = function(N = 1, lambda = 0.01, x_start, x_end, y_start, y_end){
  list_pp = list()
  K_vec = vector()
  for (i in 1:N) {
    pp = rpoispp(lambda, win=owin(c(x_start,x_end),c(y_start,y_end)))
    list_pp[[i]] = pp
    K_vec[i] = pp$n 
  }
  return( list("list_pp" = list_pp, "K_vec" = K_vec))
}

# function to set up a grid for each cell 
get_grid_nodes_list = function(grid_nrows = 5, grid_ncols = 4, J, x_start, x_end, y_start, y_end){
  x_nodes = seq(from = x_start, to = x_end, length.out = J+1)
  grid_nodes_list = list()
  for (j in 1:J) {
    left_x_node =  x_nodes[j]
    right_x_node = x_nodes[j+1]
    cell_partition = raster(nrows = grid_nrows, ncols = grid_ncols, 
                            xmn = left_x_node, xmx = right_x_node,
                            ymn = y_start, ymx = y_end)
    grid_nodes = xyFromCell(cell_partition, cell = 1:(grid_nrows*grid_ncols)) # nodes of the cell grid
    grid_nodes_list[[j]] = grid_nodes
  }
  return(grid_nodes_list)
}

# function to get the distance matrix based on the cell grid
get_Dlist = function(N = 1, J = 24, list_pp, x_start, x_end, y_start, y_end){
  list_grid_nodes = get_grid_nodes_list(J = J, x_start = x_start, x_end = x_end, 
                                        y_start = y_start, y_end = y_end)
  list_D = list() 
  for (i in 1:N) {
    pp = list_pp[[i]] # get the i^th pp object
    n_knot = pp$n
    D_mat = matrix(NA, nrow = J, ncol = n_knot)
    for (j in 1:J) {
      cell_grid_node = list_grid_nodes[[j]]
      dist_to_cell_j = vector()
      if (n_knot > 0) {
        for(k in 1:n_knot) {
          dist_to_nodes = crossdist(X = cell_grid_node[,'x'], Y = cell_grid_node[,'y'], 
                                    x2 = pp$x[k], y2 = pp$y[k])
          dist_k = mean(dist_to_nodes)
          dist_to_cell_j[k] = dist_k
        }
        D_mat[j,] = dist_to_cell_j
      }
    }
    list_D[[i]] = D_mat
  }
  return(list_D)
}



# function to calculated the weights 
get_Wlist = function(N = 1, beta , d_max, list_D){
  list_W = list() 
  for (i in 1:N) {
    D = list_D[[i]]
    #   W = D^(-beta)*as.numeric(D<= dmax)
    W = exp(-beta*D)*as.numeric(D<= d_max)
    list_W[[i]] = W
  }
  return(list_W)
}



# function to randomly generate the knot effects (percentage)
get_Zlist = function(K_vec){
  # random generate of the percentage
  # x_sections = as.vector(read.csv(file = "xsections.csv", sep = ",")[,2])
  list_Z = lapply(K_vec, 
                  function(K){ return(matrix(rgamma(n=K, shape=2, scale=6), nrow=K, ncol=1))})
                   # return(matrix(runif(n = K, min = 20, max = 120), 
                                  #nrow = K, ncol = 1))})

  return(list_Z)  
}

# function to generate the edge knot indicators
get_Elist = function(K_vec){
  list_E = lapply(K_vec, function(K){
    return(matrix(sample(c(0,1), size = K, replace = T),
                  nrow = K, ncol = 1))})
  return(list_E)
}
  
# function to scale the knot effects by edge indicator 
get_scaled_Zlist = function(gamma0, gamma1, list_Z, list_E){
  N = length(list_Z)
  list_scaled_Z = list()
  for (i in 1:N) {
    Z_i = list_Z[[i]]
    K_vec_i = nrow(Z_i)
    E_i = list_E[[i]]
    gamma_i = diag(gamma1, nrow = K_vec_i, ncol = K_vec_i)%*%E_i + 
              diag(gamma0, nrow = K_vec_i, ncol = K_vec_i)%*%(1-E_i)
    scaled_Z_i = gamma_i*Z_i
    list_scaled_Z[[i]] = scaled_Z_i
  }
  return(list_scaled_Z)
}


# simulate the data for using STAN 
simDataSTAN = function(eta0, eta1, rho, sigma, beta, gamma0, gamma1,
                       N, J, d_max, lambda = 0.015, x_start, x_end, y_start, y_end){
  K_output = vector()
  Z_output = list()
  D_output = list()
  E_output = list()
  Y = matrix(NA, ncol = N, nrow = J)
  MOE_vec = vector()
  scaled_weighted_effects =  matrix(NA, ncol = N, nrow = J) 
  # each column is the weighted effect for each cell
  # unscaled.
#  p_mat1 = matrix(NA, ncol = N, nrow = J-1)
#  p_mat2 = matrix(NA, ncol = N, nrow = J-1 )
#  p_mat3 = matrix(NA, ncol = N, nrow = J-1 )
#  p_mat4 = matrix(NA, ncol = N, nrow = J-1 )
#  p_mat5 = matrix(NA, ncol = N, nrow = J-1 )
  
  
  i = 1
  while(i <= N){
    show(i)
    PP_obj = get_PPlist(N = 1, lambda = lambda, x_start = x_start, x_end = x_end, 
                        y_start = y_start, y_end = y_end)
    list_pp = PP_obj[['list_pp']]
    K_vec = as.integer(PP_obj[['K_vec']])
    if (K_vec > 0) {
      list_D = get_Dlist(N = 1, J = J, list_pp = list_pp, 
                         x_start = x_start, x_end = x_end, 
                         y_start = y_start, y_end = y_end) 
      list_W = get_Wlist(N = 1, beta = beta, d_max = d_max, list_D = list_D)      
      list_Z = get_Zlist(K_vec = K_vec)
      list_E = get_Elist(K_vec = K_vec)
      list_scaled_Z = get_scaled_Zlist(gamma0, gamma1, list_Z, list_E)
    } else{ # there are no knots on this piece
      list_D <- list_Z <- list_E <- list()
      list_D[[1]] = matrix(10000, nrow = J, ncol = 1)
      list_W = get_Wlist(N = 1, beta = beta, d_max = d_max, list_D = list_D)      
      list_Z[[1]] = matrix(0, nrow=1, ncol=1)
      list_E[[1]] = matrix(0, nrow=1, ncol=1)
      list_scaled_Z = get_scaled_Zlist(gamma0, gamma1, list_Z, list_E)
    }
    
    MOE_i = round(rnorm(1, 1.9, 0.25),digits=2) #round(runif(n = 1, min = 1, max = 2.5), digits = 2)
    y_vec = vector(length = J)
    mean1 = (eta0 + eta1 * MOE_i) - as.vector(list_W[[1]]%*%list_scaled_Z[[1]])[1]
    margin_sd = sigma/sqrt(1-rho^2)
    y_vec[1] = rnorm(n = 1, mean = mean1, sd = margin_sd)
#    p1 = vector()
#    p2= vector()
#    p3 = vector()
#    p4 = vector()
#    p5 = vector()
    for (j in 2:J) {
      # some testing
      p1 = (1-rho)*(eta0+ eta1*MOE_i)
      p2 = rho* y_vec[j-1]
      p3 = - as.vector(list_W[[1]]%*%list_scaled_Z[[1]])[j]
      p4 =  rho * as.vector(list_W[[1]]%*%list_scaled_Z[[1]])[j-1]
      p5 =  rnorm(n = 1, mean = 0, sd = sigma)
      
      # y_vec[j] = (1-rho)*(eta0+ eta1*MOE_i) + rho* y_vec[j-1]
      # - as.vector(list_W[[1]]%*%list_scaled_Z[[1]])[j]
      # + rho * as.vector(list_W[[1]]%*%list_scaled_Z[[1]])[j-1]
      # + rnorm(n = 1, mean = 0, sd = sigma)
      y_vec[j] = p1 + p2 + p3 + p4 + p5
    }
    
    if(all(y_vec >0)){
#      p_mat1[,i] = p1
#      p_mat2[,i] = p2
#      p_mat3[,i] = p3
#      p_mat4[,i] = p4
#      p_mat5[,i] = p5
      Y[,i] = y_vec
      scaled_weighted_effects[,i] = list_W[[1]]%*%list_scaled_Z[[1]]
      MOE_vec[i] = MOE_i
      K_output[i] = K_vec[1]
      Z_output[[i]] = list_Z[[1]]
      D_output[[i]] = list_D[[1]]
      E_output[[i]] = list_E[[1]]
      i = i+1
    } else {
      show("neg")
    }
  }
  
  min_index = Rfast::colMins(Y, value = FALSE)
  Ym = Rfast::colMins(Y,value = TRUE)

return(list(
#    'p_mat1' = p_mat1,
#    'p_mat2' = p_mat2,
#    'p_mat3' = p_mat3,
#    'p_mat4' = p_mat4,
#    'p_mat5' = p_mat5,
    'Y_matrix' = Y,
    'min_index' = min_index,
    'Ym' = Ym,
    'K_vec' = K_output,
    'MOE_vec' = MOE_vec,
    'list_Z' = Z_output,
    'list_D' = D_output,
    'list_E' = E_output,
    'scaled_weighted_effects' = scaled_weighted_effects
  ))
}  
  
 





