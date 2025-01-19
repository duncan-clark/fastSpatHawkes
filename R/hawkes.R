#' @title Log likelihood for spatio temporal log likelihood
#' @description Calculates the log likelihood of a hawkes process realization
#' @param params list with elements mu,alpha,beta,K
#' @param realiz list with elements, n,lon,lat,t
#' @param windowT vector with elements, start,end
#' @param windowS spatstat.geom owin object or an object that can be coerced as such
#' @param dists list with elements space_dist,time_dist
#' @param zero_background_region spatstat.geom owin object
#' @param density_approx boolean to determine if we are using the density approximation
#' @param numeric_integral boolean to determine if we are using the numeric integral
#' @return numerical loglikelihood value
#' @export
loglik_hawk = function(params,
                       realiz,
                       windowT,
                       windowS,
                       optimized = F,
                       exponential_kernel = TRUE,
                       kernel = NULL,
                       dists = NULL,
                       zero_background_region = NULL,
                       density_approx = TRUE,
                       numeric_integral= FALSE
                       ){
  tval <- windowT[2]-windowT[1]
  max_t <- max(realiz$t)
  if(tval < max(realiz$t) - min(realiz$t)){
    stop("realization has points outside time window")
  }
  mu<-params[1]
  alpha<-params[2]
  beta<-params[3]
  K<-params[4]

  # don't allow negative parameters
  if(min(mu,K,alpha,beta)<0) return(-99999)
  # don't allow explosive growth
  if(K>.99999) return(-99999)
  # don't allow k=0
  if(K==0){
    return(-99999)
  }

  # Throw warnings if not optimzed bu the data is not too big:
  if(optimized == FALSE & dim(realiz)[1] <= 10000){
    warning("loglik_hawk is not optimized and the data has less that 10000 rows
         it should be possible to use the optimized version without running out of memory.")
  }

  # define the kernal if exponential
  if(exponential_kernel){
    kernel <- function(time_dist,space_dist,beta = beta,alpha = alpha){
      exp(-beta * time_dist - alpha * space_dist)
    }
  }

  # ==============================
  # Noramlizing Integral Portion
  # ==============================

  # Adjust for areas with no background intensity
  if(!is.null(zero_background_region)){
    area_zero_background <- spatstat.geom::area(owin(zero_background_region))
    adjust_factor <- (spatstat.geom::area.owin(windowS) - area_zero_background) / spatstat.geom::area.owin(windowS)
    in_zero_background <- 1*inside.owin(realiz[,c("x","y")],w=zero_background_region)
  }else{
    area_zero_background <- 0
    adjust_factor <- 1
    in_zero_background <- rep(0,dim(realiz)[1])
  }
  mu_star <- mu/(spatstat.geom::area(windowS))
  const  <- K*(alpha/pi)*beta


  if(density_approx){
    # For intlam to work according to
    # Facilitated estimation of ETAS. Bulletin of the Seismological Society of America, 103(1), 601-605
    # We need \int_{xy} \mu(x,y,t) dx,dy = \mu T
    # This means we need some rescaling of mu in later in the code
    # that is if the window is not of area 1 the "mu" in the likelihood is rescaled below
    # K is E # points, alpha and beta from exponential decay
    # note that we need alpha/pi to make sure the triggering function is a density
    intlam <- adjust_factor*mu*tval + K*dim(realiz)[1]
  }else{
    intlam <- adjust_factor*mu*tval
    if(exponential_kernel){
      # neat trick to calculate the integral when we have an exponential kernel
      # do analytic method:
      # each piece consider the integral only caring about the k-th time point
      func_3 <- function(x,y,t){
        t_comp <- (1-exp(-beta*(max_t-t)))
        s_1 <- (pnorm(windowS$xrange[2],mean = x,sd = 1/sqrt(2*alpha)) -
                  pnorm(0,mean = x,sd = 1/sqrt(2*alpha)))
        s_2 <- (pnorm(windowS$yrange[2],mean = y,sd = 1/sqrt(2*alpha)) -
                  pnorm(0,mean = y,sd = 1/sqrt(2*alpha)))
        return(t_comp*s_1*s_2)
      }
      pieces <- func_3(realiz$x,realiz$y,realiz$t)


    }else{
      warning("loglik_hawk is doing numerical integration this will probably take a long time")
      int_func <- function(x,realiz){
        t<- x[1]
        y<- x[3]
        x<- x[2]
        x_diff <- x - realiz$x
        y_diff <- y - realiz$y
        t_diff <- t - realiz$t
        space_dist <- x_diff^2 + y_diff^2
        time_dist <- t_diff
        # cut down the dist matrices based on integral
        result <- sum(kernel(t_diff,space_dist,beta,alpha))
        if((result==0)){
          return(1e-10)
        }
        return(result)
      }
      func <- function(i){
        max_t <- realiz$t[i]
        min_t <- realiz$t[i-1]
        hcubature(int_func,
                  realiz = realiz[realiz$t==min_t,],
                  lowerLimit = c(min_t, windowS$xrange[1], windowS$yrange[1]),
                  upperLimit = c(max_t, windowS$xrange[2], windowS$yrange[2])
        )$integral
      }
      pieces <- sapply(2:dim(realiz)[1],func)
    }
    if(K!=0 & (dim(realiz)[1] !=1)){
      intlam <- intlam + K*sum(pieces)
    }
  }

  # ==============================
  # Kernel Weighted Sum Portion
  # ==============================
  # initialize log sum:
  sum_log <- log(mu_star)
  if(optimized){
    # Precompute distances if not supplied:
    if(is.null(dists)){
      realiz <- realiz[order(realiz$t),]
      x_diff <- outer(realiz$x, realiz$x, "-")
      y_diff <- outer(realiz$y, realiz$y, "-")
      space_dist <- x_diff^2 + y_diff^2
      time_dist <- outer(realiz$t, realiz$t, "-")
    }else{
      space_dist <- dists$space_dist
      time_dist <- dists$time_dist
    }
    # check that the kernel takes in both time and space distances
    if(length(formals(kernel)) != 4){
      stop("kernel function must take in both time and space distances as wells as beta and alpha params")
    }
    kernel_mat <- kernel(time_dist,space_dist,beta,alpha)
    # for numerical stability
    kernel_mat[upper.tri(kernel_mat, diag = FALSE)] <- 0
    nzbg <- 1 - in_zero_background
    # for row column j sum across columns up to j - 1
    gij_vec <- rowSums(
      lower.tri(kernel_mat, diag = FALSE) * kernel_mat
    )
    # only credit points not in zero background region with mu
    # don't credit first row of gij vec

    lamjs <- (1-in_zero_background)*mu_star + (const * gij_vec)
    lamjs <- lamjs[2:length(lamjs)]

    # add the log lamjs in:
    if(any(is.na(lamjs)) || any(lamjs < 0)){
      return(-99999)
    }else{
      sum_log <- sum_log + sum(log(lamjs))
    }

    }else{
      if(nrow(realiz) >= 2){
        for (j in 2:nrow(realiz)){
          gij <- 0
          for (i in 1:(j-1)) {
            space_dist <- (realiz$x[j] - realiz$x[i])^2 + (realiz$y[j] - realiz$y[i])^2
            time_dist <- realiz$t[j] - realiz$t[i]
            gij <- gij + kernel(time_dist,space_dist,beta = beta,alpha = alpha)
          }
          lamj <- mu_star + const * gij
          if (is.na(lamj) || lamj < 0) {
            return(-99999)
          }
          sum_log <- sum_log + log(lamj)
        }
      }else{
        sum_log <- 0
      }
    }

  loglik <- sum_log - intlam
  if(loglik == -Inf){
    return(-99999)
  }
  return(loglik)
}

# roxgen documentation
#' @title Using nelder mead to estimate Hawkes process MLEs
#' @description Fit a hawkes process by maximizing th log likelihood with the nelder mead algorithm
#' @param params_init list with elements,
#' @param realiz list with elements, n,lon,lat,t
#' @param windowT vector with elements, start,end
#' @param windowS os.win
#' @param trace trace for optim
#' @param maxit maxit for optim
#' @param poisson_flag flag for if this process is really a pure poisson process
#' @param ... other arguments to pass to optim
#' @return optim fit object
#' @export
fit_hawkes <- function(params_init,
                       realiz,
                       windowT,
                       windowS,
                       maxit = 1000,
                       exponential_kernel = TRUE,
                       kernel = NULL,
                       optimized = TRUE,
                       trace = 0,
                       poisson_flag = FALSE,
                       zero_background_region = NULL,
                       ...){
  if(class(params_init) == "list"){params_init <- unlist(params_init)}
  if(poisson_flag | params_init[4] == 0){
    if(is.null(zero_background_region)){
      return(list(par = list(mu = dim(realiz)[1]/windowT[2] ,alpha = 1,beta = 1,K = 0), converged = T))
    }else{
      win_area <- spatstat.geom::area(windowS)
      non_zero_area <- spatstat.geom::area.owin(zero_background_region)
      adjust_factor <- (win_area - non_zero_area) / win_area
      return(list(par = list(mu = dim(realiz)[1]/(adjust_factor*windowT[2]) ,alpha = 1,beta = 1,K = 0), converged = T))
    }

  }else{
    if(optimized == FALSE & dim(realiz)[1] <= 10000){
      warning("loglik_hawk is not optimized and the data has less that 10000 rows
         it should be possible to use the optimized version without running out of memory.")
    }
    if(optimized == TRUE & dim(realiz)[1] > 20000){
      warning(paste0("Fitting with optmized=TRUE results in calculating a n x n matrix, n= ",dim(relaiz)[1]))
    }
    # order the realization
    realiz <- realiz[order(realiz$t),]
    # calculate distances if optimized
    if(optimized){
      x_diff <- outer(realiz$x, realiz$x, "-")
      y_diff <- outer(realiz$y, realiz$y, "-")
      space_dist <- x_diff^2 + y_diff^2
      time_dist <- outer(realiz$t, realiz$t, "-")
      dists <- list(space_dist = space_dist,time_dist = time_dist)
    }else{
      dists <- NULL
    }
    fit <- optim(par = params_init,
                 fn = loglik_hawk,
                 method = "Nelder-Mead", # since no hessian is available
                 control = list(fnscale = -1,trace=trace,maxit=maxit),
                 realiz = realiz,
                 windowT = windowT,
                 windowS = windowS,
                 optimized = optimized,
                 exponential_kernel = exponential_kernel,
                 kernel = kernel,
                 dists = dists,
                 zero_background_region = zero_background_region,
                 ...)
  }
  return(fit)
}

# roxgen documentation
#' @title Simulate Spatio Temporal Hawkes Processes
#' @description
#' simulating a recursive Hawkes process with lambda(t,x,y) = mu + SUM g(t_i,x_i,y_i) for g some kernel function
#' where:
#' - g(t) is kenel function with parameters alpha and beta
#' - mu is the background rate,
#' - K is the average number of aftershocks generated by an event.
#' @param params list with elements,
#' @param windowT vector with elements, start,end
#' @param windowS owin object
#' @param background_realization list with elements, n,lon,lat,t which is used as the background if supplied
#' @param filtration - points outside the time window that can trigger new points still (ie. a past history)
#' @param optimized boolean to determine if we are using the optimized version (i.e. slight faster due to vectorization)
#' @return list with t,x,y for point process
#' @export
sim_hawkes = function(params,
                      windowT,
                      windowS,
                      background_realization = NULL,
                      filtration = NULL,
                      optimized = TRUE){
  mu<-params$mu
  alpha<-params$alpha
  beta<-params$beta
  K<-params$K

  # Initialize the output list of events
  events = list()
  events$n = 0
  events$t = numeric()

  # Generate background events
  if(is.null(background_realization)){
    n_bg = rpois(1, mu * (windowT[2] - windowT[1]))
    events$n = n_bg
    events$x = runif(n_bg, min=windowS[1], max=windowS[2])
    events$y = runif(n_bg, min=windowS[3], max=windowS[4])
    events$t = sort(runif(n_bg, min=0, max=windowT[2]))
    events$background = rep(TRUE,n_bg)
  }else{
    events <- background_realization
    events$n <- length(events$t)
    events$t <- sort(events$t)
    events$background = rep(TRUE,length(events$t))
    # remove points outside window:
    valid = events$x >= windowS[1] &
      events$x <= windowS[2] &
      events$y >= windowS[3] &
      events$y <= windowS[4] &
      events$t <= windowT[2]
    events = lapply(events,function(x){x[valid]})
  }
  if(!is.null(filtration)){
    # these events can keep triggering - but new events are only added if they are affter the start time
    events <- bind_rows(as.data.frame(events),as.data.frame(filtration))
    events <- as.list(events)
  }

  # Initialize the queue with background events
  event_queue <- data.table(time = events$t, x = events$x, y = events$y,background = events$background)
  # Initialize the list to store new events
  new_events_list <- vector("list", nrow(event_queue)*(1/(1-K))*2)  # Pre-allocate space for new events, expected number of points *2 for safety
  list_index <- 1
  tot <- 0
  tot_attempt <- 0
  if(!optimized){
    # Non - optimized version - one aftershock at a time
    while (nrow(event_queue) > 0) {
      current_event <- event_queue[1, ]
      event_queue <- event_queue[-1, ]  # Remove the processed event

      num_aftershocks <- rpois(1, K)
      if (num_aftershocks > 0) {
        # NOTE the 1/beta and just alpha
        # This is due to coordinate change - DOUBLE CHECK!
        delay <- rexp(num_aftershocks, rate = beta)
        dist <- sqrt(rexp(num_aftershocks, rate = alpha))
        angle <- runif(num_aftershocks, min = 0, max = 2 * pi)
        new_lon <- current_event$x + dist * cos(angle)
        new_lat <- current_event$y + dist * sin(angle)
        new_time <- current_event$time + delay

        valid <- new_lon >= windowS[1] &
          new_lon <= windowS[2] &
          new_lat >= windowS[3] &
          new_lat <= windowS[4] &
          new_time <= windowT[2] &
          new_time >= windowT[1]

        valid = which(valid)
        tot_attempt <- tot_attempt + num_aftershocks
        tot <- tot + length(valid)
        if(tot > tot_attempt){browser()}

        if(length(valid) !=0){
          events$x <- c(events$x, new_lon[valid])
          events$y <- c(events$y, new_lat[valid])
          events$t <- c(events$t, new_time[valid])
          events$background <- c(events$background, rep(FALSE,length(valid)))
          events$n <- events$n + length(valid)

          # Add the valid new events to the queue
          new_events_list[[list_index]] <- data.table(time = new_time[valid],
                                                      x = new_lon[valid],
                                                      y = new_lat[valid],
                                                      background = rep(FALSE,length(valid)))
          list_index <- list_index + 1

        }
      }
      # Concatenate new events to event_queue only if we have only one event left to go
      # note hawkes process has no inhibition, i.e. points only trigger other points, so can do in any time order
      if (nrow(event_queue)[1] <= 1 & length(new_events_list) > 1) {
        old_n <- dim(event_queue)[1]
        event_queue <- rbindlist(list(event_queue, rbindlist(new_events_list, use.names = TRUE)))
        new_events_list <- vector("list", length(new_events_list))  # Reset the list
        list_index <- 1
      }
    }
  }else{
    # Optmized - vectorize the aftershock generation:
    while (nrow(event_queue) > 0) {
      num_events <- nrow(event_queue)
      num_aftershocks <- rpois(num_events, K)

      # Filter events with non-zero aftershocks
      idx_nonzero <- which(num_aftershocks > 0)
      if (length(idx_nonzero) == 0) {
        break
      }

      # Repeat parent events according to the number of aftershocks
      parent_x <- rep(event_queue$x[idx_nonzero], num_aftershocks[idx_nonzero])
      parent_y <- rep(event_queue$y[idx_nonzero], num_aftershocks[idx_nonzero])
      parent_time <- rep(event_queue$time[idx_nonzero], num_aftershocks[idx_nonzero])

      total_aftershocks <- length(parent_x)

      # Simulate aftershock properties
      delay <- rexp(total_aftershocks, rate = beta)
      dist <- sqrt(rexp(total_aftershocks, rate = alpha))
      angle <- runif(total_aftershocks, min = 0, max = 2 * pi)
      new_lon <- parent_x + dist * cos(angle)
      new_lat <- parent_y + dist * sin(angle)
      new_time <- parent_time + delay

      # Validate new events
      valid <- new_lon >= windowS[1] &
        new_lon <= windowS[2] &
        new_lat >= windowS[3] &
        new_lat <= windowS[4] &
        new_time <= windowT[2] &
        new_time >= windowT[1]

      if (any(valid)) {
        valid_new_lon <- new_lon[valid]
        valid_new_lat <- new_lat[valid]
        valid_new_time <- new_time[valid]
        num_valid_new_events <- length(valid_new_lon)

        # Add valid new events to events
        events$x <- c(events$x, valid_new_lon)
        events$y <- c(events$y, valid_new_lat)
        events$t <- c(events$t, valid_new_time)
        events$background <- c(events$background, rep(FALSE, num_valid_new_events))
        events$n <- events$n + num_valid_new_events

        # Update event_queue with valid new events
        event_queue <- data.table(time = valid_new_time,
                                  x = valid_new_lon,
                                  y = valid_new_lat,
                                  background = rep(FALSE, num_valid_new_events))
      } else {
        # No valid new events
        event_queue <- data.table()  # Empty
      }
    }
  }
  # Final events list
  events$n = rep(length(events$t),length(events$t))
  events
}
