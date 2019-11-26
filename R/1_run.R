#-----------------#
# 1. Load data #### 
#-----------------#

# 1.a Packages, paths and source function script
rm(list = ls())
require(raster)
source(file.path(".", "R", "0_functions.R"))

data_path <- file.path(".", "data")

# 1.b Load covariates and land use
# covariates
dat <- readRDS(file.path(data_path, "cov.rds")) #dynamic bioclimatic variables

# observed land use
lu_all <- readRDS(file.path(data_path, "lu.rds"))
lu <- lu_all[,1:5] #There are 5 classes: urban, crop, grassland, open forest, closed forest. Other columns are time steps.

# mask for Sweden and Luxemburg
mask <- readRDS(file.path(data_path, "mask.rds"))

# 1.c Calculate neighbourhood statistics
weights <- list(matrix(1, 3, 3, byrow= TRUE)) #size of window
weights <- rep(weights, length.out = 5)
ln <- neighbourhood(lu, c(1 :5), weights, mask) #see functions script to see what's happening here
data <- cbind(dat, ln)

# 1.d Correlation analysis and subset uncorr preds
preds <- colnames(correlations(data)) #this does correlation analysis and removes predictors from correlated pairs.
data <- data[,colnames(data)%in%preds]

#-------------------------#
# 2. Suitability model ####
#-------------------------#

# 2.a Model formula (simple linear terms)
form <- paste(preds, collapse="+")

# 2.b Build model
subs <- sample(1:nrow(lu), 100000)
suitmod <- suitmodel(form = form, lu = lu[subs,], data = data[subs,], resolution = 1000, model = FALSE, maxit = 1000)

# 2.c Extract demand change from observed time steps (by country)
ts <- c(2000, 2006, 2012, 2018)
demand_all <- demand(landuse = lu_all, ts = ts ,inds = NULL, k = 5)[,1:(5+1)] #demand function pulls demand change from observed time steps and interpolates missing years

#----------------------------#
# 3. Allocation algorithm ####
#----------------------------#

# 3.a Parameters for allocation algorithm: these can and probably need to be adjusted to achieve convergence and improve realism
allocpars <- list(
  stepsi = 1, # Stepsize for iterator. smaller values, slower convergence, higher chance for solution
  max_dev = 1, # Maximumn allowed deviation from prescribed demand.
  max_change = c(1,1,1,1,1), #by how much can land use grow in a cell? 1 means a land use can grow to occupy entire cell in one time step.
  min_change = - c(0.01, 0.05, 1, 1, 1), # opposite of above, i.e. urban and crop can only decrease by tiny fractions
  no_change = c(0.8, 0, 0, 0.05, 0), # growth threshold: in this case, when i.e. only when urban neighbourhood is > 0.8, urban can be allocated in a cell that currently has 0 urban. Prevents i.e. urban growth in the middle of a forest.
  suit_adj = 1e-14, # Amount of landcape-wide adjustment of predicted suitability when demand allocations increase the difference between allocated and prescribed supply.
  resolution = 10000
)

ts <- c(2000, 2006, 2012, 2018)
lu_out <- lu
lu_ts <- list()

# 3.b Iterate through time steps
i <- 2
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  # 3.c Caluclate neighbourhood maps for current time step
  ln <- neighbourhood(lu_out, 1:5, weights, mask)
  
  #3.d Predict the suitability model, adding the neigbourhood maps as covariates
  sm <- predict(suitmod, newdata = cbind(dat, ln), type = "probs")
  
  lu_pred <- matrix(NA, nrow(lu), ncol(lu))

  # 3.d get demand change for time step
  dmd <- demand_all #j instead of 12
  dmd <- dmd[which(dmd[,1]%in%ts),-1]
  dmd_t0 <- dmd[i,]
  
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  # 3.e Allocate, see functions script for what this does
  lu_pred <- allocation(lu = lu_out,
                        ln = ln, 
                        sm = sm, 
                        params = allocpars, 
                        dmd = dmd_ts)
  
  
  cat('\n')
  # 3.f Store predicted maps
  lu_ts[[i]] <- lu_out <- lu_pred
}

test <- mask
test[!is.na(test[])] <- lu_ts[[3]][,1]
plot(test)
