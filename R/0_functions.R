#Four functions

# calculate neihgourhood rasters
# <- lu data
# <- weights list
# <- which classes
# <- mask
# neighbourhood(lu, weights, classes, mask)
# -> neighbourhood data

neighbourhood <- function(lu, cols, weights, mask){
  if(length(which(!is.na(mask[]))) != nrow(lu)) stop("mask and land use not same length")
  out <- lu[,cols]
  for (i in 1:length(cols)){
    w <- weights[[i]]
    l <- mask
    inds <- which(!is.na(getValues(l)))
    l[inds] <- lu[,i]
    l <- raster::focal(l, w, fun = mean, na.rm = TRUE, pad = TRUE)
    out[,i] <- getValues(l)[inds]
  }
  colnames(out) <- paste0(colnames(lu[,cols]), "_neigh")
  out
}

#Correlation analysis
# <- predictor matrix
# <- subset
# <- treshold
# -> predictor names

correlations <- function(covs, thresh = 0.7) {
  subs_cor <- sample(1:nrow(covs), size = 10000)
  cors <- cor(covs[subs_cor, ], method = "spearman")
  while (min(abs(cors[abs(cors) >= thresh])) != 1){
    values <- cors[which(abs(cors) > thresh)]
    corellated <- which(abs(cors) > thresh)
    values[values ==1] <- NA
    corellated[which(values== max(values, na.rm = T))]
    rows_highest_cor <- which(cors == max(values, na.rm = T), arr.ind = T)[,1]
    cors_cur <- abs(cors[rows_highest_cor,])
    '%ni%' <- Negate('%in%')
    m1 <- max(cors_cur[1,][cors_cur[1,]%ni%c(max(values, na.rm = T),1)])
    m2 <- max(cors_cur[2,][cors_cur[2,]%ni%c(max(values, na.rm = T),1)])
    out <- ifelse(m1 > m2, 1, 2)
    cors <- cors[-which(colnames(cors) == names(rows_highest_cor)[out]), -which(colnames(cors) == names(rows_highest_cor)[out])]
    nrow(cors)
  }
  return(cors)
}

# build and predict suitaboility models
# <- subset length
# <- lu
# <- cov
# <- neigh rasters
# <- resolution for "integerification"
# <- formula object
# suitability(lu, cov, ind, neigh)
# -> predicted suitability maps

suitmodel <- function(form, lu, data, sub = NULL, resolution,...){
  
  if(!is.null(sub)){
    subs <- sample(1:nrow(lu), sub)
  }else{
    subs <- 1:nrow(lu)
  }
  
  counts <- integerify(x = lu[subs,], resolution = resolution)
  data_sub <- as.data.frame(data[subs,])
  f <- as.formula(paste("counts", "~", form))
  suit_model <- nnet::multinom(f, data = data_sub, ...)
  suit_model
}

# estimate demand from existing land-use data
# <- indices for subsetting
# <- land use maps for different time steps
# demand(lu, ind)
# -> demand trajectory 
#landuse <- lu_all

demand <- function(inds = NULL, landuse, ts, path = NULL, k){
  #Calculate demand
  lu <- landuse
  if(!is.null(inds)){
    lu <- lu[inds, ]
  }
  class_demand <- colMeans(lu)
  
  test <- numeric()
  years <- ts[1]:tail(ts,1)
  demand <- matrix(NA, nrow = length(years), ncol = k + 2)
  demand[,1] <- years
  obs_ind <- which(demand[,1]%in%ts)
  i <- 1
  for(i in 1:(ncol(lu)/k)){
    demand[obs_ind[i], 2:6]<- class_demand[seq((i * k - 4), i * k, by = 1)]
  }
  
  for(i in 1:k){
    demand[, i + 1] <- approx(demand[,1], demand[,i + 1], xout = years)$y
  }
  
  demand[, k + 2] <- rowSums(demand[,-1], na.rm = TRUE)
  if(!is.null(path)){
    saveRDS(demand, path)
  }
  return(demand)
}

# allocate in landscape
# <- land-use data first time step
# <- suitability maps
# <- constraints (maximum/minimum change, threshold and which class is tresholded)
# <- demand
# allocate(land-use data, suitability maps, constraints, demand)
# -> predicted land-use
# lu <- lu[subs,]
# ln <- ln[subs,]
# sm <- sm[subs,]
# params <- allocpars
# dmd <- dmd_ts

allocation <- function(lu, ln, sm, params, dmd){
  dev_diff = rep(10,5)
  resolution <- params$resolution
  p_t <- p_t0 <- lu
  p_t1_candidate <- p_t0 <- integerify(p_t0, resolution = resolution) #function in functions script
  demand_t0 <- apply(p_t0, 2, sum)
  demand_traj <- integerify(x = dmd, resolution = sum(demand_t0))
  demand_t1 <- demand_traj[2,]
  
  candidate_supply <- supply_t0 <- demand_t0
  k <- ncol(lu)
  n <- nrow(lu)
  
  max_dev <- params$max_dev
  stepsize <- params$stepsi
  min_change <- params$min_change * params$resolution
  max_change <- params$max_change * params$resolution
  no_change <- params$no_change
  
  count <- 0
  proposal_valid <- FALSE
  
  demand_change <- demand_t1 - candidate_supply
  mindiff <- pmin(demand_change, rep(0,k))
  maxdiff <- pmax(demand_change, rep(0,k))
  cat("\n")
  while ((sum(dev_diff > rep(max_dev,k)) !=0 | !proposal_valid) & count < 20000) {
    #6.a Counter increment
    count <- count + 1
    
    #6.b Supply that still needs to be allocated in each class (equal to demand_t-1 - demand_t in first iteration)
    demand_change <- demand_t1 - candidate_supply
    
    #6.c Estimate ideal change
    ideal_amounts <- sm
    if(any(pmax(demand_change, maxdiff) > maxdiff)){
      sa <- which(maxdiff < demand_change)
      ideal_amounts[,sa] <- ideal_amounts[,sa] + params$suit_adj
      maxdiff <- pmin(demand_change, maxdiff)
    }
    
    if(any(pmin(demand_change, mindiff) < mindiff)){
      sa <- which(mindiff > demand_change)
      ideal_amounts[,sa] <- ideal_amounts[,sa] - params$suit_adj
      mindiff <- pmax(demand_change, mindiff)
    }
    
    #ii Calculate change factor, i.e. by how much do we have to multiply the observed land use proportions to satisfy the modelled suitability.
    ideal_change <- (ideal_amounts * resolution) / p_t1_candidate
    
    #iii Handle divide by 0 errors. We get these when original landuse was 0. We just add 1 + predicted suitability
    bad <- !is.finite(ideal_change)
    ideal_change[bad] <- ideal_amounts[bad]
    
    #6.d Estimate how much to add to each cell
    #i Relative suitability (like maxent, sums to 1 in each col)
    rel_suitability <- ideal_change %*% diag(1/colSums(ideal_change))
    
    #ii Allocate demand change between pixels, wieghted by rel_suitability and adjusted by stepsize
    target_lu_change_pixel <-  rel_suitability %*% diag(demand_change)  * stepsize
    
    #6.e Apply constraints
    #i Maximum increase/decrease for each cell and type
    for (class in seq_len(k)) {
      target_lu_change_pixel[, class] <- pmin(target_lu_change_pixel[, class], max_change[class])
      target_lu_change_pixel[, class] <- pmax(target_lu_change_pixel[, class], min_change[class])
    }
    
    #ii Add proposed changes and make everyting positive
    p_t1_proposal <- p_t1_candidate + target_lu_change_pixel
    p_t1_proposal <- pmax(p_t1_proposal, 0)
    
    
    #6.f Turn proposal into integers by sampling: this adds stochasticity and makes sure values across land-use classes sum to defined resolution
    # - integrify function samples integer vlaues from mnomial distr, with relative probabilities given by p_t1_proposal and also takes care of land-uses that can only change in one direction, or not at all (they are excluded from sampling).
    
    p_t1_proposal <- integerify(x = p_t1_proposal, z = p_t1_candidate, resolution = resolution, no_decrease = min_change >= 0)
    
    #6.g Allow growth only into empty cells above certain threshold in neighrbourhood raster
    if(!is.null(no_change)){
      ch_classes <- which(no_change != 0)
      chl <- length(ch_classes)
      for(ch in 1:chl){
        p_t1_proposal[which(ln[,ch] < no_change[ch_classes[ch]] & p_t1_candidate[,ch_classes[ch]] == 0), ch] <- 0
      }
    }
    
    #6.h Check min/max change condition (we want to make sure, changes don't exceed what we said is possible)
    #i Calculate currently proposed cell-level changes
    proposed_change <- p_t1_proposal - p_t
    
    #ii Check that each proposed change added to a cell is withhin class min/max constraints
    proposal_valid <- TRUE
    for (class in seq_len(k)) {
      class_change_proposal <- proposed_change[, class]
      valid_class <- all((class_change_proposal <= max_change[class]) & (class_change_proposal >= min_change[class]))
      if (!valid_class) {
        proposal_valid <- FALSE
      }
    }
    
    #6.h Propose new candidate land-use distribution
    #i Current proposal becomes new candidate only if min/max constraints are met, otherwise the iteration is redone exactly as before, with a new multinomial draw that hopefully meets conditions.
    #print(proposal_valid)
    if (proposal_valid) {
      p_t1_candidate <- p_t1_proposal
    }
    
    #ii Calcuate new candidate supply, i.e. the supply of the currently proposed candidate
    candidate_supply <- colSums(p_t1_candidate)
    
    #iii Work out how close candidate supply is to demand
    diff <- abs(demand_t1 - candidate_supply)
    dev_diff <- diff/demand_t1 * 100
    dev_diff[which(is.na(dev_diff))] <- 0
    #iv change the effect of stepsize iterator (make finer), if we keep getting the exact same multinbomial draws that lead to the same, not-quite-good-enough land-use maps over and over again.
    
    if(count == 1){
      dev_diff_cur <- dev_diff
      mean_counter <- 1
    }
    
    if(count > 1){
      if(all(dev_diff_cur == dev_diff)){
        mean_counter <- mean_counter + 1
      }
      if(all(dev_diff_cur != dev_diff)){
        mean_counter <- 0 
        dev_diff_cur <- dev_diff
      }
      
      if(mean_counter == 3){
        stepsize <- stepsize/10
        mean_counter <- 0
        dev_diff_cur <- dev_diff
        print(paste0("changed stepsize to ", stepsize))
      }
    }
    if(stepsize < 1e-10){
      stepsize <- params$stepsi
      print(paste0("reset stepsize"))
    }
    cat("\r", paste0("Iteration: ", count, "    "), "Deviation from target per class [%]: ", paste(round(dev_diff, 3), sep = " "))
    
    if (count == 20000) {
      stop("algorithm failed to converge")
    }
  }
  pred_out <- p_t1_candidate/resolution
  colnames(pred_out) <- colnames(lu)
  pred_out
}

# autocovariates
# weights <- list of weights for neighbourhood rasters
# neigh_classes <- which of lu_cols to estimate for?
# neigh_data1 <- first ts neighbourhood raster
# neigh_data <- neighbourhood raster
# data <- data stack for estimations


integerify <- function (x, z = NULL, resolution = resolution, no_decrease = NULL) {
  
  n <- nrow(x)
  k <- ncol(x)
  
  if (any(no_decrease)) {
    min_allocation <- matrix(0, n, k)
    for (class in which(no_decrease)) {
      min_allocation[, class] <- z[, class]
    }
  } else {
    min_allocation <- NULL
  }
  
  # if a minimum number must be placed in a certain column, withold them and just allocate the others
  if (!is.null(min_allocation)) {
    to_allocate <- resolution - rowSums(min_allocation)
    
    # need to reduce the probability of allocating the remainder to these cells too
    # so make x sum to resolution, subtract the minima, and make sure the result is positive
    x <- resolution * x / rowSums(x) #* resolution because gets scaled to 0-1
    x <- x - min_allocation
    x <- pmax(x, 0)
    
  } else {
    to_allocate <- rep(resolution, n)
  }
  
  # random draws of the non-reserved ones
  x <- extraDistr::rmnom(n, to_allocate, prob = x)
  
  # add on the reserved ones again
  if(!is.null(min_allocation)) {
    x <- x + min_allocation
  }
  x
}