## ---------------------------
##
## Script name: entomological_model
##
## Purpose of script: define functions to initialize le model and to run it
##
## Author: Alex DROUIN
##
## Date Created: 2023-07-05

# NB : theta is the array of temperature data (row and columns : latitude and longitude, slices : time steps)

################################################################################
# Discrete time model ##########################################################
################################################################################

model_initialization <- function(theta,
                                 nE = 0,
                                 nL = 0,
                                 nP = 0,
                                 nA_em = 0,
                                 nA_1h = 0,
                                 nA_1g = 0,
                                 nA_1o = 0,
                                 nA_2h = 0,
                                 nA_2g = 0,
                                 nA_2o = 0){
  nrows = dim(theta)[1]
  ncols = dim(theta)[2]
  
  init <- list(E = matrix(nE, nrow = nrows, ncol = ncols),
               L = matrix(nL, nrow = nrows, ncol = ncols),
               P = matrix(nP, nrow = nrows, ncol = ncols),
               A_em = matrix(nA_em, nrow = nrows, ncol = ncols),
               A_1h = matrix(nA_1h, nrow = nrows, ncol = ncols),
               A_1g = matrix(nA_1g, nrow = nrows, ncol = ncols),
               A_1o = matrix(nA_1o, nrow = nrows, ncol = ncols),
               A_2h = matrix(nA_2h, nrow = nrows, ncol = ncols),
               A_2g = matrix(nA_2g, nrow = nrows, ncol = ncols),
               A_2o = matrix(nA_2o, nrow = nrows, ncol = ncols))
  
  return(init)
  
}

################################################################################
# Function to initialize model by repeating loops on 1990 ######################
################################################################################
#we want to run the model on the first year (1990) until the difference between two following yearly runs is negligible

#on two consequent repetitions of 1990, we compute the difference, if > 0.01% we run the model again using the previous output as initialization
#the output is the 31 dec of the year, which is used as initialization for the following run



loop_on_first_year_initialization <- function(log_alphaa, link_kappa_p, threshold, species, initialization, theta, R,
                                              season_start, season_end, proba_presence, timeline, timestep, output_folder,
                                              output_file_note, africa_id_mask){
  
  start_time <- Sys.time()
  
  #we start with 2 repetitions of year 1 on which we run the model, and repeating while distance > 0.01%
  year1 <- unique(year(timeline))[1]
  timeline_year1 <- timeline[year(timeline) == year1]
  theta_year1 <- theta[,, timeline %in% timeline_year1, drop = FALSE]
  R_year1 <- R[,, timeline %in% timeline_year1, drop = FALSE]
  
  
  if(species %in% c("caspius", "detritus", "vexans")){
    first_initialization <- model_initialization(theta = theta_year1, nE = 10^4)
  }else{
    first_initialization <- model_initialization(theta = theta_year1, nA_em = 10^4)
  }
  
  #we run the model a first time to have an initialization
  entomological_model_discrete_time(species = species, 
                                    initialization = first_initialization, 
                                    theta = theta_year1, R = R_year1,
                                    season_start = season_start, season_end = season_end,
                                    proba_presence = proba_presence,
                                    timeline = timeline_year1,
                                    timestep = timestep,
                                    output_folder = output_folder,
                                    output_file_note = output_file_note,
                                    log_alphaa = log_alphaa,
                                    follow_up_comments = "none",
                                    link_kappa_p = link_kappa_p, 
                                    africa_id_mask = africa_id_mask)
  #opening outputs
  state_E <- read.csv(paste("outputs/", output_folder, species,"_output_E_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_L <- read.csv(paste("outputs/", output_folder, species,"_output_L_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_P <- read.csv(paste("outputs/", output_folder, species,"_output_P_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_em <- read.csv(paste("outputs/", output_folder, species,"_output_A_em_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_1h <- read.csv(paste("outputs/", output_folder, species,"_output_A_1h_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_1g <- read.csv(paste("outputs/", output_folder, species,"_output_A_1g_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_1o <- read.csv(paste("outputs/", output_folder, species,"_output_A_1o_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_2h <- read.csv(paste("outputs/", output_folder, species,"_output_A_2h_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_2g <- read.csv(paste("outputs/", output_folder, species,"_output_A_2g_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  state_A_2o <- read.csv(paste("outputs/", output_folder, species,"_output_A_2o_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
  
  state_E <- as.matrix(state_E)
  state_L <- as.matrix(state_L)
  state_P <- as.matrix(state_P)
  state_A_em <- as.matrix(state_A_em)
  state_A_1h <- as.matrix(state_A_1h)
  state_A_1g <- as.matrix(state_A_1g)
  state_A_1o <- as.matrix(state_A_1o)
  state_A_2h <- as.matrix(state_A_2h)
  state_A_2g <- as.matrix(state_A_2g)
  state_A_2o <- as.matrix(state_A_2o)
  
  #saving data of december 31st
  state_E_31_dec <- state_E[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_L_31_dec <- state_L[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_P_31_dec <- state_P[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_em_31_dec <- state_A_em[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_1h_31_dec <- state_A_1h[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_1g_31_dec <- state_A_1g[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_1o_31_dec <- state_A_1o[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_2h_31_dec <- state_A_2h[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_2g_31_dec <- state_A_2g[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  state_A_2o_31_dec <- state_A_2o[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
  
  #removing outputs
  rm(first_initialization, state_E, state_L, state_P, state_A_em, state_A_1h, state_A_1g, state_A_1o, state_A_2h, state_A_2g, state_A_2o)
  
  #preparing vectors to store the distances
  dist_E_31_dec <- c()
  dist_L_31_dec <- c()
  dist_P_31_dec <- c()
  dist_A_em_31_dec <- c()
  dist_A_1h_31_dec <- c()
  dist_A_1g_31_dec <- c()
  dist_A_1o_31_dec <- c()
  dist_A_2h_31_dec <- c()
  dist_A_2g_31_dec <- c()
  dist_A_2o_31_dec <- c()
  dist_E_31_dec_relative_per_pixel <- c()
  dist_L_31_dec_relative_per_pixel <- c()
  dist_P_31_dec_relative_per_pixel <- c()
  dist_A_em_31_dec_relative_per_pixel <- c()
  dist_A_1h_31_dec_relative_per_pixel <- c()
  dist_A_1g_31_dec_relative_per_pixel <- c()
  dist_A_1o_31_dec_relative_per_pixel <- c()
  dist_A_2h_31_dec_relative_per_pixel <- c()
  dist_A_2g_31_dec_relative_per_pixel <- c()
  dist_A_2o_31_dec_relative_per_pixel <- c()
  
  nb_pixel_data <- sum(!is.na(state_E_31_dec[1,]))
  
  #initialization of relative distance at 100%   
  relative_distance <- 100 
  number_loop <- 0
  initialization_loop <- list(E = matrix(state_E_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              L = matrix(state_L_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              P = matrix(state_P_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_em = matrix(state_A_em_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_1h = matrix(state_A_1h_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_1g = matrix(state_A_1g_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_1o = matrix(state_A_1o_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_2h = matrix(state_A_2h_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_2g = matrix(state_A_2g_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                              A_2o = matrix(state_A_2o_31_dec[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]))
  
  rm(state_E_31_dec,
     state_L_31_dec,
     state_P_31_dec,
     state_A_em_31_dec,
     state_A_1h_31_dec,
     state_A_1g_31_dec,
     state_A_1o_31_dec,
     state_A_2h_31_dec,
     state_A_2g_31_dec,
     state_A_2o_31_dec)
  
  while(relative_distance > threshold){
    
    print(paste("loop number : ", number_loop))
    
    entomological_model_discrete_time(species = species, 
                                      initialization = initialization_loop, 
                                      theta = theta_year1, R = R_year1,
                                      season_start = season_start, season_end = season_end,
                                      proba_presence = proba_presence,
                                      timeline = timeline_year1,
                                      timestep = timestep,
                                      output_folder = output_folder,
                                      output_file_note = output_file_note,
                                      log_alphaa = log_alphaa,
                                      follow_up_comments = "none",
                                      link_kappa_p = link_kappa_p, 
                                      africa_id_mask = africa_id_mask)
    
    
    #we compute the relative distance between year i+1 and year 1
    
    state_E_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_E_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_L_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_L_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_P_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_P_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_em_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_em_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_1h_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_1h_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_1g_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_1g_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_1o_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_1o_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_2h_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_2h_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_2g_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_2g_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    state_A_2o_year_iplus1 <- read.csv(paste("outputs/", output_folder, species,"_output_A_2o_discrete_time_timestep_0.04_days_", output_file_note, ".csv", sep = ""), header = FALSE)
    
    #saving data of december 31st
    state_E_31_dec_year_iplus1 <- state_E_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_L_31_dec_year_iplus1 <- state_L_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_P_31_dec_year_iplus1 <- state_P_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_em_31_dec_year_iplus1 <- state_A_em_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_1h_31_dec_year_iplus1 <- state_A_1h_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_1g_31_dec_year_iplus1 <- state_A_1g_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_1o_31_dec_year_iplus1 <- state_A_1o_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_2h_31_dec_year_iplus1 <- state_A_2h_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_2g_31_dec_year_iplus1 <- state_A_2g_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    state_A_2o_31_dec_year_iplus1 <- state_A_2o_year_iplus1[day(timeline_year1) == 31 & month(timeline_year1) == 12,, drop = FALSE]
    
    state_E_31_dec_year_iplus1 <- as.matrix(state_E_31_dec_year_iplus1)
    state_L_31_dec_year_iplus1 <- as.matrix(state_L_31_dec_year_iplus1)
    state_P_31_dec_year_iplus1 <- as.matrix(state_P_31_dec_year_iplus1)
    state_A_em_31_dec_year_iplus1 <- as.matrix(state_A_em_31_dec_year_iplus1)
    state_A_1h_31_dec_year_iplus1 <- as.matrix(state_A_1h_31_dec_year_iplus1)
    state_A_1g_31_dec_year_iplus1 <- as.matrix(state_A_1g_31_dec_year_iplus1)
    state_A_1o_31_dec_year_iplus1 <- as.matrix(state_A_1o_31_dec_year_iplus1)
    state_A_2h_31_dec_year_iplus1 <- as.matrix(state_A_2h_31_dec_year_iplus1)
    state_A_2g_31_dec_year_iplus1 <- as.matrix(state_A_2g_31_dec_year_iplus1)
    state_A_2o_31_dec_year_iplus1 <- as.matrix(state_A_2o_31_dec_year_iplus1)
    
    
    
    rm(state_E_year_iplus1,
       state_L_year_iplus1, 
       state_P_year_iplus1,
       state_A_em_year_iplus1,
       state_A_1h_year_iplus1,
       state_A_1g_year_iplus1,
       state_A_1o_year_iplus1,
       state_A_2h_year_iplus1,
       state_A_2g_year_iplus1,
       state_A_2o_year_iplus1)
    
    
    #reshapping in matrix lat*lon
    states_31_dec_year_iplus1 <- list(E = matrix(state_E_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      L = matrix(state_L_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      P = matrix(state_P_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_em = matrix(state_A_em_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_1h = matrix(state_A_1h_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_1g = matrix(state_A_1g_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_1o = matrix(state_A_1o_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_2h = matrix(state_A_2h_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_2g = matrix(state_A_2g_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]),
                                      A_2o = matrix(state_A_2o_31_dec_year_iplus1[1,], nrow = dim(theta)[1], ncol = dim(theta)[2]))
    rm(state_E_31_dec_year_iplus1,
       state_L_31_dec_year_iplus1,
       state_P_31_dec_year_iplus1,
       state_A_em_31_dec_year_iplus1,
       state_A_1h_31_dec_year_iplus1,
       state_A_1g_31_dec_year_iplus1,
       state_A_1o_31_dec_year_iplus1,
       state_A_2h_31_dec_year_iplus1,
       state_A_2g_31_dec_year_iplus1,
       state_A_2o_31_dec_year_iplus1)
    
    #E
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_E_31_dec_year_i <- initialization_loop[["E"]]
    
    dist_E_31_dec_i <- sum(abs(states_31_dec_year_iplus1$E - state_E_31_dec_year_i), na.rm = TRUE)
    dist_E_31_dec <- c(dist_E_31_dec, dist_E_31_dec_i)
    
    #mean value per pixel 
    mean_E_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$E), na.rm = TRUE)
    mean_dist_E_31_dec_i <- dist_E_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_E_31_dec_relative_per_pixel_i <- (mean_dist_E_31_dec_i / mean_E_31_dec_year_iplus1 )*100
    dist_E_31_dec_relative_per_pixel <- c(dist_E_31_dec_relative_per_pixel, dist_E_31_dec_relative_per_pixel_i)
    
    #L
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_L_31_dec_year_i <- initialization_loop[["L"]]
    
    dist_L_31_dec_i <- sum(abs(states_31_dec_year_iplus1$L - state_L_31_dec_year_i), na.rm = TRUE)
    dist_L_31_dec <- c(dist_L_31_dec, dist_L_31_dec_i)
    
    #mean value per pixel 
    mean_L_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$L), na.rm = TRUE)
    mean_dist_L_31_dec_i <- dist_L_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_L_31_dec_relative_per_pixel_i <- (mean_dist_L_31_dec_i / mean_L_31_dec_year_iplus1 )*100
    dist_L_31_dec_relative_per_pixel <- c(dist_L_31_dec_relative_per_pixel, dist_L_31_dec_relative_per_pixel_i)
    
    #P
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_P_31_dec_year_i <- initialization_loop[["P"]]
    
    dist_P_31_dec_i <- sum(abs(states_31_dec_year_iplus1$P - state_P_31_dec_year_i), na.rm = TRUE)
    dist_P_31_dec <- c(dist_P_31_dec, dist_P_31_dec_i)
    
    #mean value per pixel 
    mean_P_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$P), na.rm = TRUE)
    mean_dist_P_31_dec_i <- dist_P_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_P_31_dec_relative_per_pixel_i <- (mean_dist_P_31_dec_i / mean_P_31_dec_year_iplus1 )*100
    dist_P_31_dec_relative_per_pixel <- c(dist_P_31_dec_relative_per_pixel, dist_P_31_dec_relative_per_pixel_i)
    
    #A_em
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_em_31_dec_year_i <- initialization_loop[["A_em"]]
    
    dist_A_em_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_em - state_A_em_31_dec_year_i), na.rm = TRUE)
    dist_A_em_31_dec <- c(dist_A_em_31_dec, dist_A_em_31_dec_i)
    
    #mean value per pixel 
    mean_A_em_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_em), na.rm = TRUE)
    mean_dist_A_em_31_dec_i <- dist_A_em_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_em_31_dec_relative_per_pixel_i <- (mean_dist_A_em_31_dec_i / mean_A_em_31_dec_year_iplus1 )*100
    dist_A_em_31_dec_relative_per_pixel <- c(dist_A_em_31_dec_relative_per_pixel, dist_A_em_31_dec_relative_per_pixel_i)
    
    #A_1h
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_1h_31_dec_year_i <- initialization_loop[["A_1h"]]
    
    dist_A_1h_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_1h - state_A_1h_31_dec_year_i), na.rm = TRUE)
    dist_A_1h_31_dec <- c(dist_A_1h_31_dec, dist_A_1h_31_dec_i)
    
    #mean value per pixel 
    mean_A_1h_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_1h), na.rm = TRUE)
    mean_dist_A_1h_31_dec_i <- dist_A_1h_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_1h_31_dec_relative_per_pixel_i <- (mean_dist_A_1h_31_dec_i / mean_A_1h_31_dec_year_iplus1 )*100
    dist_A_1h_31_dec_relative_per_pixel <- c(dist_A_1h_31_dec_relative_per_pixel, dist_A_1h_31_dec_relative_per_pixel_i)
    
    #A_1g
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_1g_31_dec_year_i <- initialization_loop[["A_1g"]]
    
    dist_A_1g_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_1g - state_A_1g_31_dec_year_i), na.rm = TRUE)
    dist_A_1g_31_dec <- c(dist_A_1g_31_dec, dist_A_1g_31_dec_i)
    
    #mean value per pixel 
    mean_A_1g_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_1g), na.rm = TRUE)
    mean_dist_A_1g_31_dec_i <- dist_A_1g_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_1g_31_dec_relative_per_pixel_i <- (mean_dist_A_1g_31_dec_i / mean_A_1g_31_dec_year_iplus1 )*100
    dist_A_1g_31_dec_relative_per_pixel <- c(dist_A_1g_31_dec_relative_per_pixel, dist_A_1g_31_dec_relative_per_pixel_i)
    
    #A_1o
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_1o_31_dec_year_i <- initialization_loop[["A_1o"]]
    
    dist_A_1o_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_1o - state_A_1o_31_dec_year_i), na.rm = TRUE)
    dist_A_1o_31_dec <- c(dist_A_1o_31_dec, dist_A_1o_31_dec_i)
    
    #mean value per pixel 
    mean_A_1o_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_1o), na.rm = TRUE)
    mean_dist_A_1o_31_dec_i <- dist_A_1o_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_1o_31_dec_relative_per_pixel_i <- (mean_dist_A_1o_31_dec_i / mean_A_1o_31_dec_year_iplus1 )*100
    dist_A_1o_31_dec_relative_per_pixel <- c(dist_A_1o_31_dec_relative_per_pixel, dist_A_1o_31_dec_relative_per_pixel_i)
    
    #A_2h
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_2h_31_dec_year_i <- initialization_loop[["A_2h"]]
    
    dist_A_2h_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_2h - state_A_2h_31_dec_year_i), na.rm = TRUE)
    dist_A_2h_31_dec <- c(dist_A_2h_31_dec, dist_A_2h_31_dec_i)
    
    #mean value per pixel 
    mean_A_2h_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_2h), na.rm = TRUE)
    mean_dist_A_2h_31_dec_i <- dist_A_2h_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_2h_31_dec_relative_per_pixel_i <- (mean_dist_A_2h_31_dec_i / mean_A_2h_31_dec_year_iplus1 )*100
    dist_A_2h_31_dec_relative_per_pixel <- c(dist_A_2h_31_dec_relative_per_pixel, dist_A_2h_31_dec_relative_per_pixel_i)
    
    #A_2g
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_2g_31_dec_year_i <- initialization_loop[["A_2g"]]
    
    dist_A_2g_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_2g - state_A_2g_31_dec_year_i), na.rm = TRUE)
    dist_A_2g_31_dec <- c(dist_A_2g_31_dec, dist_A_2g_31_dec_i)
    
    #mean value per pixel 
    mean_A_2g_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_2g), na.rm = TRUE)
    mean_dist_A_2g_31_dec_i <- dist_A_2g_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_2g_31_dec_relative_per_pixel_i <- (mean_dist_A_2g_31_dec_i / mean_A_2g_31_dec_year_iplus1 )*100
    dist_A_2g_31_dec_relative_per_pixel <- c(dist_A_2g_31_dec_relative_per_pixel, dist_A_2g_31_dec_relative_per_pixel_i)
    
    #A_2o
    ###computing distance between 31st december of year  i (that we just computed) and year i-1 (ie the initialization_loop)
    state_A_2o_31_dec_year_i <- initialization_loop[["A_2o"]]
    
    dist_A_2o_31_dec_i <- sum(abs(states_31_dec_year_iplus1$A_2o - state_A_2o_31_dec_year_i), na.rm = TRUE)
    dist_A_2o_31_dec <- c(dist_A_2o_31_dec, dist_A_2o_31_dec_i)
    
    #mean value per pixel 
    mean_A_2o_31_dec_year_iplus1 <- mean(as.numeric(states_31_dec_year_iplus1$A_2o), na.rm = TRUE)
    mean_dist_A_2o_31_dec_i <- dist_A_2o_31_dec_i / nb_pixel_data
    
    #relative value per pixel
    dist_A_2o_31_dec_relative_per_pixel_i <- (mean_dist_A_2o_31_dec_i / mean_A_2o_31_dec_year_iplus1 )*100
    dist_A_2o_31_dec_relative_per_pixel <- c(dist_A_2o_31_dec_relative_per_pixel, dist_A_2o_31_dec_relative_per_pixel_i)
    
    relative_distance <- max(dist_E_31_dec_relative_per_pixel_i,
                             dist_L_31_dec_relative_per_pixel_i,
                             dist_P_31_dec_relative_per_pixel_i,
                             dist_A_em_31_dec_relative_per_pixel_i,
                             dist_A_1h_31_dec_relative_per_pixel_i,
                             dist_A_1g_31_dec_relative_per_pixel_i,
                             dist_A_1o_31_dec_relative_per_pixel_i,
                             dist_A_2h_31_dec_relative_per_pixel_i,
                             dist_A_2g_31_dec_relative_per_pixel_i,
                             dist_A_2o_31_dec_relative_per_pixel_i)
    
    initialization_loop <- states_31_dec_year_iplus1
    
    
    
    number_loop <- number_loop + 1
  }
  
  # #plotting
  # distance_df_31_dec_relative_mean <- data.frame(year = 1:number_loop,
  #                                                dist_E = dist_E_31_dec_relative_per_pixel, 
  #                                                dist_L = dist_L_31_dec_relative_per_pixel,
  #                                                dist_P = dist_P_31_dec_relative_per_pixel,
  #                                                dist_A_em = dist_A_em_31_dec_relative_per_pixel,
  #                                                dist_A_1h = dist_A_1h_31_dec_relative_per_pixel,
  #                                                dist_A_1g = dist_A_1g_31_dec_relative_per_pixel,
  #                                                dist_A_1o = dist_A_1o_31_dec_relative_per_pixel,
  #                                                dist_A_2h = dist_A_2h_31_dec_relative_per_pixel,
  #                                                dist_A_2g = dist_A_2g_31_dec_relative_per_pixel,
  #                                                dist_A_2o = dist_A_2o_31_dec_relative_per_pixel)
  # 
  # distance_df_31_dec_relative_mean <- gather(distance_df_31_dec_relative_mean, key = "state", value = "relative_mean_distance_per_pixel", -year)
  # 
  # #relative mean distance per pixel 31 dec
  # plot_distance <- ggplot(distance_df_31_dec_relative_mean, aes(x = ymd(date), y = relative_mean_distance_per_pixel, color = state)) +
  #   geom_line() +
  #   labs(x ="number of loops of year 1990", y = "difference as a % of the total value of year i+1 (mean per pixel)", fill = "State") +
  #   ggtitle("Total distance between predictions year i+1 and year i (on 31 dec)")
  # 
  # ggsave(paste("outputs/optim_alpha/plots/distance_between_years_", species,"_31_dec_relative_mean.jpeg", sep = ""),
  #        width = 29.7,
  #        height = 21)
  # 
  
  end_time <- Sys.time()
  time_to_compute <- end_time - start_time 
  
  cat(paste("Initialization done in ", round(time_to_compute, 2),  units(time_to_compute), ". Number of loops: ", number_loop, ".\n", sep = ""))
  
  
  return(initialization_loop)
}

################################################################################
# Discrete time model ##########################################################
################################################################################

entomological_model_discrete_time <- function(species, theta, R, initialization, proba_presence,
                                              season_start, season_end, africa_id_mask, timeline, 
                                              timestep, output_folder = "", output_file_note = "",
                                              log_alphaa = NA, follow_up_comments = "full",
                                              link_kappa_p){
  
  #time balises
  start_time_total <- Sys.time()
  
  #defining bisextile year, number of year in the data
  bisextile_years <- c(1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)
  
  #transforming timeline to Date if it is not the case, and extracting years
  timeline <- ymd(timeline)
  timeline_years <- unique(year(timeline))
  nyear <- length(timeline_years)
  
  
  #model initializating
  output_timestep_k_minus_1 <- initialization
  
  # creating files for the output (data will be appended each day)
  #creating output directory if not existing
  dir.create(file.path(paste("outputs/", output_folder, sep = "")), showWarnings = FALSE)
  
  #preparing names
  output_files_names <- c("output_E", 
                          "output_L", 
                          "output_P", 
                          "output_A_em", 
                          "output_A_1h", 
                          "output_A_1g", 
                          "output_A_1o", 
                          "output_A_2h", 
                          "output_A_2g",
                          "output_A_2o")
  
  for (i in 1:length(output_files_names)){
    #completing files names
    output_files_names[i] <- paste("outputs/", output_folder, species, "_", output_files_names[i], "_discrete_time_timestep_", round(timestep, 2), "_days_", output_file_note, ".csv", sep = "")
    #creating files
    file.create(output_files_names[i])
  }
  
  #opening connections for apending
  con_E <- file(output_files_names[1], "a")
  con_L <- file(output_files_names[2], "a")
  con_P <- file(output_files_names[3], "a")
  con_A_em <- file(output_files_names[4], "a")
  con_A_1h <- file(output_files_names[5], "a")
  con_A_1g <- file(output_files_names[6], "a")
  con_A_1o <- file(output_files_names[7], "a")
  con_A_2h <- file(output_files_names[8], "a")
  con_A_2g <- file(output_files_names[9], "a")
  con_A_2o <- file(output_files_names[10], "a")
  
  #storing into a list
  con_list <- list(con_E,
                   con_L,
                   con_P,
                   con_A_em,
                   con_A_1h,
                   con_A_1g,
                   con_A_1o,
                   con_A_2h,
                   con_A_2g,
                   con_A_2o)
  
  #expliciting parameters
  beta_1 <- parameters["beta_1", species]
  beta_2 <- parameters["beta_2", species]
  sigma <- parameters["sigma", species]
  gamma_Aem <- parameters["gamma_Aem", species]
  gamma_Ah <- parameters["gamma_Ah", species]
  gamma_Ao <- parameters["gamma_Ao", species]
  mu_E <- parameters["mu_E", species]
  mu_L <- parameters["mu_L", species]
  mu_P <- parameters["mu_P", species]
  mu_A <- parameters["mu_A", species]
  mu_em <- parameters["mu_em", species]
  mu_r <- parameters["mu_r", species]
  phi <- parameters["phi", species]
  
  if (link_kappa_p == "kappa_constant"){
    kappa_L <- parameters["kappa_L", species] 
    kappa_P <- parameters["kappa_P", species]
  } else if(link_kappa_p == "linear"){
    kappa_L <- parameters["kappa_L", species]*as.matrix(proba_presence) 
    kappa_P <- parameters["kappa_P", species]*as.matrix(proba_presence)
  } else if(link_kappa_p == "exponential"){
    kappa_L <- parameters["kappa_L", species]^(as.matrix(proba_presence)^exp(log_alphaa)) 
    kappa_P <- parameters["kappa_P", species]^(as.matrix(proba_presence)^exp(log_alphaa)) 
  }
  
  africa_id_mask <- as.matrix(africa_id_mask)
  
  
  if(!is.null(dim(kappa_L)) & any(dim(kappa_L) != dim(theta)[1:2])){
    stop("Error : wrong dim(kappa_L)")
  }
  
  if(!is.null(dim(kappa_P)) & any(dim(kappa_P) != dim(theta)[1:2])){
    stop("Error : wrong dim(kappa_P)")
  }
  
  if (any(c(beta_1, beta_2, sigma, gamma_Aem, gamma_Ah, gamma_Ao, mu_E, mu_L, mu_P, mu_A, mu_em, mu_r, kappa_L, kappa_P, phi) < 0, na.rm = TRUE)){
    stop("One fixed parameter is < 0")
  }
  
  # temporal loop on years i (useful for diapause periods)
  for (i in 1:nyear){
    
    #time balise
    start_time_year_i <- Sys.time()
    
    #defining year_i, vector of days of year i, and number of days of year i
    year_i <- timeline_years[i]
    days_year_i <- year(timeline) == year_i
    ndays_year_i <- ifelse(year_i %in% bisextile_years,
                           366, 
                           365)
    
    # subscript of theta and R for the year i
    theta_i <- theta[,,days_year_i, drop = FALSE]
    R_i <- R[,,days_year_i, drop = FALSE]
    
    #temporal loop on days j of year i
    for (j in 1:ndays_year_i){
      
      #time balise
      start_time_day_j <- Sys.time()
      
      # extracting temperature and rainfall for day j
      theta_i_j <- theta_i[,,j, drop = FALSE]
      theta_i_j <- matrix(theta_i_j, nrow = dim(theta_i_j)[1], ncol = dim(theta_i_j)[2])
      
      R_i_j <- R_i[,,j, drop = FALSE]
      R_i_j <- matrix(R_i_j, nrow = dim(R_i_j)[1], ncol = dim(R_i_j)[2])
      
      #updating functions 
      f_E_i_j <- f_E(species = species, theta = theta_i_j, R = R_i_j)   
      f_L_i_j <- f_L(species = species, theta = theta_i_j) 
      f_P_i_j <- f_P(species = species, theta = theta_i_j)
      f_Ag_i_j <- f_Ag(species = species, theta = theta_i_j)
      m_L_i_j <- m_L(species = species, theta = theta_i_j)
      m_P_i_j <- m_P(species = species, theta = theta_i_j)
      m_A_i_j <- m_A(species = species, theta = theta_i_j)
      m_A_dia_i_j <- m_A_dia(species = species, theta = theta_i_j)
      
      if(any(c(f_E_i_j,f_L_i_j, f_P_i_j, f_Ag_i_j, m_L_i_j, m_P_i_j, m_A_i_j, m_A_dia_i_j) < 0, na.rm = TRUE)){
        stop("One function < 0")
      }
      
      
      # temporal loop on each k timestep of the day
      
      for (k in 1:(1/timestep)){
        
        #print(paste("i=", i, "/j=", j, "/k=", k))
        
        #defining current values of states variables
        curE = output_timestep_k_minus_1$E
        curL = output_timestep_k_minus_1$L
        curP = output_timestep_k_minus_1$P
        curA_1o = output_timestep_k_minus_1$A_1o
        curA_2o = output_timestep_k_minus_1$A_2o
        curA_em = output_timestep_k_minus_1$A_em
        curA_1h = output_timestep_k_minus_1$A_1h
        curA_1g = output_timestep_k_minus_1$A_1g
        curA_1o = output_timestep_k_minus_1$A_1o
        curA_2h = output_timestep_k_minus_1$A_2h
        curA_2g = output_timestep_k_minus_1$A_2g
        curA_2o = output_timestep_k_minus_1$A_2o
        
        if(any(curE < 0, na.rm = TRUE)){value = curE[which(curE < 0)] ; stop(paste("One state variable < 0: curE. Value = ", value, "\n"))}
        if(any(curL < 0, na.rm = TRUE)){value = curL[which(curL < 0)] ; stop(paste("One state variable < 0: curL. Value = ", value, "\n"))}
        if(any(curP < 0, na.rm = TRUE)){value = curP[which(curP < 0)] ; stop(paste("One state variable < 0: curP. Value = ", value, "\n"))}
        if(any(curA_em < 0, na.rm = TRUE)){value = curA_em[which(curA_em < 0)] ; stop(paste("One state variable < 0: curA_em. Value = ", value, "\n"))}
        if(any(curA_1h < 0, na.rm = TRUE)){value = curA_1h[which(curA_1h < 0)] ; stop(paste("One state variable < 0: curA_1h. Value = ", value, "\n"))}
        if(any(curA_1g < 0, na.rm = TRUE)){value = curA_1g[which(curA_1g < 0)] ; stop(paste("One state variable < 0: curA_1g. Value = ", value, "\n"))}
        if(any(curA_1o < 0, na.rm = TRUE)){value = curA_1o[which(curA_1o < 0)] ; stop(paste("One state variable < 0: curA_1o. Value = ", value, "\n"))}
        if(any(curA_2h < 0, na.rm = TRUE)){value = curA_2h[which(curA_2h < 0)] ; stop(paste("One state variable < 0: curA_2h. Value = ", value, "\n"))}
        if(any(curA_2g < 0, na.rm = TRUE)){value = curA_2g[which(curA_2g < 0)] ; stop(paste("One state variable < 0: curA_2h. Value = ", value, "\n"))}
        if(any(curA_2o < 0, na.rm = TRUE)){value = curA_2o[which(curA_2o < 0)] ; stop(paste("One state variable < 0: curA_2o. Value = ", value, "\n"))}
        
        
        favourable_season = (j >= season_start) & (j <= season_end)
        
        # updating states changes: E
        if (favourable_season) {
          E.L <- (1 - exp(-f_E_i_j * timestep)) * curE
        } else {
          E.L <- matrix(NA, nrow = nrow(curE), ncol = ncol(curE))
          
          E.L <- ifelse(africa_id_mask == 0, 
                        phi * (1 - exp(-f_E_i_j * timestep)) * curE,
                        ifelse(africa_id_mask == 1, 
                               (1 - exp(-f_E_i_j * timestep)) * curE, 
                               NA))
        }
        
        .E <- (1 - exp(-gamma_Ao * timestep)) * 
          (beta_1 * curA_1o + beta_2 * curA_2o)
        E. <- (1 - exp(-mu_E * timestep)) * curE
        out <- E.L + E.
        
        if (any(out > curE, na.rm = TRUE)) {
          x <- which(out > curE)
          E.L[x] <- E.L[x] * curE[x] / out[x]
          E.[x] <- E.[x] * curE[x] / out[x]
        }
        
        # updating states changes: L
        L.P <- (1 - exp(-f_L_i_j * timestep)) * curL
        L. <- (1 - exp(-timestep * (m_L_i_j * (1 + curL/kappa_L)))) * curL
        out <- L.P + L.
        
        if (any(out > curL, na.rm = TRUE)) {
          x <- which(out > curL)
          L.P[x] <- L.P[x] * curL[x] / out[x]
          L.[x] <- L.[x] * curL[x] / out[x]
        }
        
        # updating states changes: P
        P.A_em <- (1 - exp(-timestep *
                             (f_P_i_j * sigma * 
                                exp(-mu_em * (1 + curP/kappa_P))))) * curP
        P. <- (1 - exp(-timestep *
                         (m_P_i_j + 
                            f_P_i_j * (
                              1 - sigma * exp(-mu_em * 
                                                (1 + curP/kappa_P)))))) * curP
        out <- P.A_em + P.
        
        if (any(out > curP, na.rm = TRUE)) {
          x <- which(out > curP)
          P.A_em[x] <- P.A_em[x] * curP[x] / out[x]
          P.[x] <- P.[x] * curP[x] / out[x]
        }
        
        
        # updating states changes: A_em
        if (favourable_season) {
          A_em.A_1h <- (1 - exp(-gamma_Aem * timestep)) * curA_em
          A_em. <- (1 - exp(-m_A_i_j * timestep)) * curA_em
        } else {
          A_em.A_1h <- matrix(NA, ncol = ncol(curA_em), nrow = nrow(curA_em))
          A_em. <- matrix(NA, ncol = ncol(curA_em), nrow = nrow(curA_em))
          
          A_em.A_1h <-ifelse(africa_id_mask == 0, 
                             (1 - phi) * (1 - exp(-gamma_Aem * timestep)) * curA_em, 
                             ifelse(africa_id_mask == 1, 
                                    (1 - exp(-gamma_Aem * timestep)) * curA_em, 
                                    NA))
          
          A_em. <- ifelse(africa_id_mask == 0, 
                          (1 - exp(-m_A_dia_i_j * timestep)) * curA_em, 
                          ifelse(africa_id_mask == 1, 
                                 (1 - exp(-m_A_i_j * timestep)) * curA_em, 
                                 NA))
        }
        out <- A_em.A_1h + A_em.
        
        if (any(out > curA_em, na.rm = TRUE)) {
          x <- which(out > curA_em)
          A_em.A_1h[x] <- A_em.A_1h[x] * curA_em[x] / out[x]
          A_em.[x] <- A_em.[x] * curA_em[x] / out[x]
        }
        
        # updating states changes: A_1h
        A_1h.A_1g <- (1 - exp(-gamma_Ah * timestep)) * curA_1h
        A_1h. <- (1 - exp(-timestep *(m_A_i_j + mu_r))) * curA_1h
        out <- A_1h.A_1g + A_1h.
        
        if (any(out > curA_1h, na.rm = TRUE)) {
          x <- which(out > curA_1h)
          A_1h.A_1g[x] <- A_1h.A_1g[x] * curA_1h[x] / out[x]
          A_1h.[x] <- A_1h.[x] * curA_1h[x] / out[x]
        }
        
        # updating states changes: A_1g
        A_1g.A_1o <- (1 - exp(-f_Ag_i_j * timestep)) * curA_1g
        A_1g. <- (1 - exp(-m_A_i_j * timestep)) * curA_1g
        out <- A_1g.A_1o + A_1g.
        
        if (any(out > curA_1g, na.rm = TRUE)) {
          x <- which(out > curA_1g)
          A_1g.A_1o[x] <- A_1g.A_1o[x] * curA_1g[x] / out[x]
          A_1g.[x] <- A_1g.[x] * curA_1g[x] / out[x]
        }
        
        # updating states changes: A_1o
        A_1o.A_2h <- (1 - exp(-gamma_Ao * timestep)) * curA_1o
        A_1o. <- (1 - exp(-timestep * (m_A_i_j + mu_r))) * curA_1o
        out <- A_1o.A_2h + A_1o.
        
        if (any(out > curA_1o, na.rm = TRUE)) {
          x <- which(out > curA_1o)
          A_1o.A_2h[x] <- A_1o.A_2h[x] * curA_1o[x] / out[x]
          A_1o.[x] <- A_1o.[x] * curA_1o[x] / out[x]
        }
        
        # updating states changes: A_2h
        A_2h.A_2g <- (1 - exp(-gamma_Ah * timestep)) * curA_2h
        A_2h. <- (1 - exp(-timestep * (m_A_i_j + mu_r))) * curA_2h
        out <- A_2h.A_2g + A_2h.
        
        if (any(out > curA_2h, na.rm = TRUE)) {
          x <- which((out > curA_2h))
          A_2h.A_2g[x] <- A_2h.A_2g[x] * curA_2h[x] / out[x]
          A_2h.[x] <- A_2h.[x] * curA_2h[x] / out[x]
        }
        
        # updating states changes: A_2g
        A_2g.A_2o <- (1 - exp(-f_Ag_i_j * timestep)) * curA_2g
        A_2g. <- (1 - exp(-m_A_i_j * timestep)) * curA_2g
        out <- A_2g.A_2o + A_2g.
        
        if (any(out > curA_2g, na.rm = TRUE)) {
          x <- which(out > curA_2g)
          A_2g.A_2o[x] <- A_2g.A_2o[x] * curA_2g[x] / out[x]
          A_2g.[x] <- A_2g.[x] * curA_2g[x] / out[x]
        }
        
        # updating states changes: A_2o
        A_2o.A_2h <- (1 - exp(-gamma_Ao * timestep)) * curA_2o
        A_2o. <- (1 - exp(-timestep * (m_A_i_j + mu_r))) * curA_2o
        out <- A_2o.A_2h + A_2o.
        
        if (any(out > curA_2o, na.rm = TRUE)) {
          x <- which(out > curA_2o)
          A_2o.A_2h[x] <- A_2o.A_2h[x] * curA_2o[x] / out[x]
          A_2o.[x] <- A_2o.[x] * curA_2o[x] / out[x]
        }
        
        # updating states variables and assessing if values < 0
        E_i_j_k <- curE + .E - E.L - E.
        if (any(E_i_j_k < 0, na.rm = TRUE)){
          
          x2 = which(E_i_j_k < 0)  
          
          if(all(abs(E_i_j_k[x2]) < 10^-6)) {
            
            E_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curE, "curE.RDS")
            saveRDS(.E, ".E.RDS")
            saveRDS(E.L, "E.L.RDS")
            saveRDS(E., "E..RDS")
            saveRDS(E_i_j_k, "E_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value = E_i_j_k[which(abs(E_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: E_i_j_k. Value = ", value, "\n"))
          }
        }
        
        
        L_i_j_k <- curL + E.L - L.P - L.
        if (any(L_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(L_i_j_k < 0)  
          
          if (all(abs(L_i_j_k[x2]) < 10^-6)) {
            
            L_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curL, "curL.RDS")
            saveRDS(E.L, "E.L.RDS")
            saveRDS(L.P, "L.P.RDS")
            saveRDS(L., "L..RDS")
            saveRDS(L_i_j_k, "L_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value = L_i_j_k[which(abs(L_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: L_i_j_k. Value = ", value, "\n"))
          }
        }
        
        
        P_i_j_k <- curP + L.P - P.A_em - P.
        if (any(P_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(P_i_j_k < 0)  
          
          if(all(abs(P_i_j_k[x2]) < 10^-6)) {
            
            P_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curP, "curP.RDS")
            saveRDS(L.P, "L.P.RDS")
            saveRDS(P.A_em, "P.A_em.RDS")
            saveRDS(P., "P..RDS")
            saveRDS(P_i_j_k, "P_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- P_i_j_k[which(abs(P_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: P_i_j_k. Value = ", value, "\n"))
          }
        }
        
        
        A_em_i_j_k <- curA_em + P.A_em - A_em.A_1h - A_em.
        if (any(A_em_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_em_i_j_k < 0)  
          
          if(all(abs(A_em_i_j_k[x2]) < 10^-6)) {
            
            A_em_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_em, "curA_em.RDS")
            saveRDS(P.A_em, "P.A_em.RDS")
            saveRDS(A_em.A_1h, "A_em.A_1h.RDS")
            saveRDS(A_em., "A_em..RDS")
            saveRDS(A_em_i_j_k, "A_em_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_em_i_j_k[which(abs(A_em_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: A_em_i_j_k. Value = ", value, "\n"))
          }
        }
        
        A_1h_i_j_k <- curA_1h + A_em.A_1h - A_1h.A_1g - A_1h.
        if (any(A_1h_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_1h_i_j_k < 0)  
          
          if(all(abs(A_1h_i_j_k[x2]) < 10^-6)) {
            
            A_1h_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_1h, "curA_1h.RDS")
            saveRDS(A_em.A_1h, "A_em.A_1h.RDS")
            saveRDS(A_1h.A_1g, "A_1h.A_1g.RDS")
            saveRDS(A_1h., "A_1h..RDS")
            saveRDS(A_1h_i_j_k, "A_1h_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_1h_i_j_k[which(abs(A_1h_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: A_1h_i_j_k. Value = ", value, "\n"))
          }
        }
        
        A_1g_i_j_k <- curA_1g + A_1h.A_1g - A_1g.A_1o - A_1g.
        if (any(A_1g_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_1g_i_j_k < 0)  
          
          if(all(abs(A_1g_i_j_k[x2]) < 10^-6)) {
            
            A_1g_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_1g, "curA_1g.RDS")
            saveRDS(A_1h.A_1g, "A_1h.A_1g.RDS")
            saveRDS(A_1g.A_1o, "A_1g.A_1o.RDS")
            saveRDS(A_1g., "A_1g..RDS")
            saveRDS(A_1g_i_j_k, "A_1g_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_1g_i_j_k[which(abs(A_1g_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: A_1g_i_j_k Value = ", value, "\n"))
          }
        }
        
        A_1o_i_j_k <- curA_1o + A_1g.A_1o - A_1o.A_2h - A_1o.
        if (any(A_1o_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_1o_i_j_k < 0)  
          
          if (all(abs(A_1o_i_j_k[x2]) < 10^-6)) {
            
            A_1o_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_1o, "curA_1o.RDS")
            saveRDS(A_1g.A_1o, "A_1g.A_1o.RDS")
            saveRDS(A_1o.A_2h, "A_1o.A_2h.RDS")
            saveRDS(A_1o., "A_1o..RDS")
            saveRDS(A_1o_i_j_k, "A_1o_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_1o_i_j_k[which(abs(A_1o_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: A_1o_i_j_k. Value = ", value, "\n"))
          }
        }
        
        A_2h_i_j_k <- curA_2h + A_1o.A_2h + A_2o.A_2h - A_2h.A_2g - A_2h.
        if (any(A_2h_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_2h_i_j_k < 0)  
          
          if (all(abs(A_2h_i_j_k[x2]) < 10^-6)) {
            
            A_2h_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_2h, "curA_2h.RDS")
            saveRDS(A_1o.A_2h, "A_1o.A_2h.RDS")
            saveRDS(A_2o.A_2h, "A_2o.A_2h.RDS")
            saveRDS(A_2h.A_2g, "A_2h.A_2g.RDS")
            saveRDS(A_2h., "A_2h..RDS")
            saveRDS(A_2h_i_j_k, "A_2h_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_2h_i_j_k[which(abs(A_2h_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: A_2h_i_j_k. Value = ", value, "\n"))
          }
        }
        
        A_2g_i_j_k <- curA_2g + A_2h.A_2g - A_2g.A_2o - A_2g.
        if (any(A_2g_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_2g_i_j_k < 0)  
          
          if (all(abs(A_2g_i_j_k[x2])) < 10^-6) {
            
            A_2g_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_2g, "curA_2g.RDS")
            saveRDS(A_2h.A_2g, "A_2h.A_2g.RDS")
            saveRDS(A_2g.A_2o, "A_2g.A_2o.RDS")
            saveRDS(A_2g., "A_2g..RDS")
            saveRDS(A_2g_i_j_k, "A_2g_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_2g_i_j_k[which(abs(A_2g_i_j_k[x2]) < 10^-6)]
            
            stop(paste("One state variable < 0: A_2g_i_j_k Value = ", value, "\n"))
            
          }
        }
        
        A_2o_i_j_k <- curA_2o + A_2g.A_2o - A_2o.A_2h - A_2o.
        if (any(A_2o_i_j_k < 0, na.rm = TRUE)){
          
          x2 <- which(A_2o_i_j_k < 0)  
          
          if (all(abs(A_2o_i_j_k[x2])) < 10^-6) {
            
            A_2o_i_j_k[x2] <- 0 
            
          } else {
            
            saveRDS(curA_2o, "curA_2o.RDS")
            saveRDS(A_2g.A_2o, "A_2g.A_2o.RDS")
            saveRDS(A_2o.A_2h, "A_2o.A_2h.RDS")
            saveRDS(A_2o., "A_2o..RDS")
            saveRDS(A_2o_i_j_k, "A_2o_i_j_k.RDS")
            saveRDS(theta_i_j, "theta_i_j.RDS")
            saveRDS(R_i_j, "R_i_j.RDS")
            
            value <- A_2o_i_j_k[which(abs(A_2o_i_j_k[x2]) < 10^-6), drop = FALSE]
            
            stop(paste("One state variable < 0: A_2o_i_j_k Value = ", value, "\n"))
            
          }
        }
        
        #updating timestep k minus 1
        output_timestep_k_minus_1 <- list(E = E_i_j_k,
                                          L = L_i_j_k,
                                          P = P_i_j_k,
                                          A_em = A_em_i_j_k,
                                          A_1h = A_1h_i_j_k,
                                          A_1g = A_1g_i_j_k,
                                          A_1o = A_1o_i_j_k,
                                          A_2h = A_2h_i_j_k,
                                          A_2g = A_2g_i_j_k,
                                          A_2o = A_2o_i_j_k)
        
        
        # store results of last timestep of the day : adding the results to the output files (one line per day i, one column per pixel)
        #first transforming output of day j from matrix to vector (reading column per colum)
        if (k == 1/timestep){
          
          output_list <- list(output_E = t(c(E_i_j_k)), 
                              output_L = t(c(L_i_j_k)), 
                              output_P = t(c(P_i_j_k)), 
                              output_A_em = t(c(A_em_i_j_k)), 
                              output_A_1h = t(c(A_1h_i_j_k)), 
                              output_A_1g = t(c(A_1g_i_j_k)), 
                              output_A_1o = t(c(A_1o_i_j_k)), 
                              output_A_2h = t(c(A_2h_i_j_k)), 
                              output_A_2g = t(c(A_2g_i_j_k)),
                              output_A_2o = t(c(A_2o_i_j_k))) 
          
          
          #then adding this line to the output file
          for (l in 1:length(output_list)){
            write.table(output_list[[l]], file = con_list[[l]], append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
          }
          
        }
      } # end loop k
      
      #cleaning memory
      rm(#erasing weather variables
        theta_i_j, R_i_j,
        #erasing output
        output_list
      )
      
      #time balise
      end_time_day_j <- Sys.time()
      time_to_compute_day_j <- end_time_day_j - start_time_day_j 
      
      if (follow_up_comments == "full"){
        cat("Done : day ", j, " of year ", year_i, " (", round(time_to_compute_day_j, 2), units(time_to_compute_day_j), ") \n", sep ="")
      }
      
    } #end loop day j
    
    #time balise
    end_time_year_i <- Sys.time()
    time_to_compute_year_i <- end_time_year_i - start_time_year_i 
    
    if (follow_up_comments == "full"){
      cat("Done : year ", year_i, " (", round(time_to_compute_year_i, 2), units(time_to_compute_year_i), ") \n", sep ="")
    }
  }
  
  #closing connections
  for (i in 1:length(con_list)){
    close(con_list[[i]])
  }
  
  #time balise
  end_time_total <- Sys.time()
  time_to_compute_total <- end_time_total - start_time_total
  
  if (follow_up_comments == "minimal"){
    cat("> Discrete model for ", species, " (timestep: ", round(timestep, 2), " days) done. (Total time = ", round(time_to_compute_total, 2), units(time_to_compute_total), ") \n", sep ="")
  }
}

