# Pre-process Douglas fir dataset

rm(list = ls())
library(ggplot2)
library(tidyverse)
library(mvtnorm)
library(spatstat)
library(raster)
library(geostatsp)

source("knot_functions.R")

# get the board list 
common_path <- "matched_enhanced_matching/"

primary_dirs = list.files(common_path, pattern = "*.csv")
board_list= lapply(primary_dirs, function(file){
  read.csv(paste0(common_path,file))
})
names(board_list) = primary_dirs

# valid test span
grip = 24 # machine grip in inches 

x_start = grip
x_end = 144 - grip

y_start =  0 
y_end =  5.5

z_start = 0
z_end = 1.5

# partition
J = 24

# board_list =list(board_list$enhanced_matching063.csv)

get_knots_volume_list = function(board_list){
  knots_volume_list = list()
  for (i in 1:length(board_list)) {
    set.seed(0)
    board = board_list[[i]]
    print(names(board_list)[i])
    # matched_knot_list can be empty if the lumber does not have knots
    matched_knot_list = get_matched_knot_list(board_dat = board, 
                                              x_start = x_start,
                                              x_end = x_end)
    # edge knot can be list of 0
    edge_knot = sapply(matched_knot_list, function(knot_dat){unique(knot_dat$edge_knot)},simplify = T)
      
    # knots_volume can be empty data frame
    knots_volume = get_knots_volume(matched_knot_list,
                                      x_start = x_start, 
                                      x_end = x_end,
                                      y_start = y_start, 
                                      y_end = y_end,
                                      z_start = z_start, 
                                      z_end = z_end)
    knots_volume$edge_knot = edge_knot
    knots_volume_list[[i]] = knots_volume
    
  }
  return(knots_volume_list)
}


add_knot_effects = function(knots_volume_list, J){
  x_len = 96  # test span in inches
  y_len = 5.5
  z_len = 1.5
  total_volume = x_len*y_len*z_len
  cell_volume = total_volume/J
  
  for (i in 1:length(knots_volume_list)) {
    knots_volume = knots_volume_list[[i]]
    knots_volume$knot_effect = knots_volume$volume/cell_volume*100 # in percentage
    knots_volume_list[[i]] <- knots_volume
  }
  
  return(knots_volume_list)
}

# function to set up a grid for each cell 
get_grid_nodes_list = function(grid_nrows = 5, grid_ncols = 4, J){
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


# calculate the mean-distance-to-grid-node distance matrix list
get_Dlist = function(knot_effects_list, J){
  list_grid_nodes = get_grid_nodes_list(J = J)
  list_D = list()
  
    for(i in 1:length(knot_effects_list)){
      knots_dat = knot_effects_list[[i]]
      n_knot = nrow(knots_dat)
      
      D_mat = matrix(NA, nrow = J, ncol = n_knot)
      if(n_knot!=0){
        for (j in 1:J) {
          cell_grid_node = list_grid_nodes[[j]]
          dist_to_cell_j = vector()
          for (k in 1:n_knot) {
            knots_dat_k = knots_dat[k,]
            dist_to_nodes = crossdist(X = cell_grid_node[,'x'], Y = cell_grid_node[,'y'], 
                                      x2 = knots_dat_k$x_coord, y2 = knots_dat_k$y_coord)
            dist_k = mean(dist_to_nodes)
            dist_to_cell_j[k] = dist_k
          }
          D_mat[j,] <- dist_to_cell_j
        }
      }
      list_D[[i]] = D_mat
    }
  return(list_D)
}

get_Zlist = function(knot_effects_list){
  list_Z = list()
  for (i in 1:length(knot_effects_list)) {
    knots_dat = knot_effects_list[[i]]
    list_Z[[i]] = as.matrix(knots_dat$knot_effect)
  }
  return(list_Z)
}

get_Edgelist = function(knot_effects_list){
  list_Edge = list()
  for(i in 1:length(knot_effects_list)){
    knots_dat = knot_effects_list[[i]]
    list_Edge[[i]] = as.matrix(knots_dat$edge_knot)
  }
  return(list_Edge)
}


knots_volume_list = get_knots_volume_list(board_list = board_list)
knot_effects_list = add_knot_effects(knots_volume_list, J = J)

save(knots_volume_list,file = paste0("knots_volume_list.RData"))
save(knots_volume_list, file = paste0("knot_effects_list.RData"))
# temp = lapply(knot_effects_list, function(x){
#  temp = all(x$knot_effect <= 100)
# })


## Data for the Bayesian Model
# N
N = 113
#N = 10

# K_vec
K_vec = unlist(lapply(knot_effects_list, function(knot_dat){nrow(knot_dat)}))

# Kmax
Kmax = max(K_vec) 

# d_max
d_max = 96 # in unit of x coordinates

# D_list
D_list = get_Dlist(knot_effects_list = knot_effects_list, J = J)

# Z_list
Z_list = get_Zlist(knot_effects_list = knot_effects_list)

# Edge Knot indicators
Edge_list = get_Edgelist(knot_effects_list = knot_effects_list)

# y_obs
# Using the UTS in imperial units
strength_data = read.csv(paste0("data/sum_data.csv"))

# sanity check for the order
library(readr)
paste0(as.character(strength_data$Board)==parse_number(primary_dirs))

Ym = strength_data$UTS_psi
Ym_ksi = strength_data$UTS_ksi
log_Ym_ksi = log(Ym_ksi)

# MOE_vec
MOE_vec = strength_data$VibeE

# min_index
failure_location = strength_data$mean_Failure_x_coord
cell_edge_x_coord = seq(x_start, x_end, length.out = J+1)
min_index = findInterval(failure_location,cell_edge_x_coord, all.inside = TRUE)

list_Z_uniform = function(N, Z_list){
  Z_list_uniform = list()
  for (i in 1:N) {
    Z_i = Z_list[[i]]
    Z_new = matrix(0, nrow = Kmax, ncol = 1)
    if(length(Z_i)!=0){
      Z_new[1:length(Z_i),] = Z_i
      }
    Z_list_uniform[[i]] = Z_new
  }
  return(Z_list_uniform)
}

list_Edge_uniform = function(N, Edge_list){
  Edge_list_uniform = list()
  for (i in 1:N) {
    E_i = Edge_list[[i]]
    E_new = matrix(-1, nrow = Kmax, ncol = 1)
    if(length(E_i)!=0){
      E_new[1:length(E_i),] = E_i
    }
    Edge_list_uniform[[i]] = E_new
  }
  return(Edge_list_uniform)
}


list_D_uniform = function(N,D_list){
  D_list_uniform = list()
  for (i in 1:N) {
    D_i = D_list[[i]]
    D_new = matrix(1e4, nrow = J, ncol = Kmax)
    if(ncol(D_i)!=0){
      D_new[1:nrow(D_i), 1:ncol(D_i)] = D_i
    }
    D_list_uniform[[i]]= D_new
  }
  return(D_list_uniform)
}


Z_list_uniform = list_Z_uniform(N = N, Z_list = Z_list)
Z_array = array(as.numeric(unlist(Z_list_uniform)), dim=c(Kmax, 1,N))
D_list_uniform = list_D_uniform(N = N, D_list = D_list)
D_array = array(as.numeric(unlist(D_list_uniform)), dim=c(J, Kmax,N))
Edge_list_uniform = list_Edge_uniform(N = N, Edge_list = Edge_list)
Edge_array = array(as.numeric(unlist(Edge_list_uniform)), dim=c(Kmax, 1,N))
# View(Edge_array)








