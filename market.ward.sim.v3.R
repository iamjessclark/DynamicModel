# function to simulate movements
# inputs:
#   - inter-node movement probabilities matrix, mu.z
#   - inter-node movement rate matrix (assuming zero-truncated 
#     negative binomial), mu.c
#   - alpha, the standard dispersion parameter (aka "size") for the
#     negative binomial distribution
#   - a logical vector indicating which nodes (rows/columns) are markets 

# the major change in v3 (2017-12-01) is to allow the mu matrices not to
# be square, i.e. not all wards can generate movements, but all wards 
# can receive movements, so nrow(mu) < ncol(mu)

sim.hurdle.move <-
  function(mu.z, mu.c, alpha, is.mkt) {
    
    # if the probability of movement is zero for all pairs of nodes, return 
    # a matrix of zero moves
    if(all(mu.z == 0)) return(mu.z)
    
    # simulate numbers being moved backwards because movements
    # out of markets are more likely to be recorded than movements into markets
    
    # create empty matrix to store one round of moves
    move.c <- move.z <- rep(NA, length = prod(dim(mu.z)))
    attributes(move.c) <- attributes(move.z) <- attributes(mu.z)
    
    # expected number moved
    (mu <- mu.z * mu.c)
    
    # which pairs had any movement at all
    move.z[, ] <- rbinom(prod(dim(mu.z)), 1, c(mu.z))
    # mu.z
    # move.z
    
    # move 4: ward-ward movements
    print("simulating ward-ward movements")
    # first, which wards have any movement
    # move.z[!is.mkt, !is.mkt] 
    
    # and how many movements are there given that movement has happened
    move.c[!is.mkt, !is.mkt] <- 
      rztnbinom.mu(n = sum(!is.mkt)^2, 
                   mu = c(mu.c[!is.mkt, !is.mkt]), 
                   size = alpha) #* c(mu.z[!is.mkt, !is.mkt])
    
    # move 3: movements from markets to wards
    print("simulating market-ward movements")
    # rates:
    # mu.z[is.mkt, !is.mkt]
    # mu.c[is.mkt, !is.mkt]
    # move.z[is.mkt, !is.mkt]
    move.c[is.mkt, !is.mkt] <- 
      rztnbinom.mu(n = sum(is.mkt) * sum(!is.mkt), mu = c(mu.c[is.mkt, !is.mkt]), 
                   size = alpha)
    (move <- move.c * move.z)
    
    # the number that moved out of each market
    (n.move3 <- sum(move[is.mkt, !is.mkt]))
    # rowSums(move[is.mkt, !is.mkt])
    
    # move 2: moves among markets
    print("simulating market-market movements")
    # assume that the between market records are as good as the out of market records,
    # so a proportion of the moves out of markets (move 3), had moved between markets
    # e.g. ...  
    # sum(move["M1", !is.mkt]) # ...animals moved out of market 1. 
    # on average, this number...
    # mu[is.mkt, "M1"]
    # came into M1 from other markets
    
    # the expected number moved among markets:
    # sum(mu[is.mkt, is.mkt])
    # expected proportion of movements from markets to wards that had moved between markets
    # sum(mu[is.mkt, is.mkt]) / sum(mu[is.mkt, !is.mkt])  #### constraint! #### must be < 1 
    
    n.move2 <- 
      if(n.move3 < 0.5) 0 else
        rbinom(1, n.move3, min(1, sum(mu[is.mkt, is.mkt]) / sum(mu[is.mkt, !is.mkt])))
    
    if(n.move2 == 0) {
      move[is.mkt, is.mkt] <- 0
      n.move1.tot <-
        rowSums(move[is.mkt, !is.mkt]) + rowSums(move[is.mkt, is.mkt]) - colSums(move[is.mkt, is.mkt])
    } else {
      repeat {
        move[is.mkt, is.mkt] <-
          rztmultinom.dir(n = 1, size = n.move2, prob = prop.table(c(mu[is.mkt, is.mkt])), 
                          alpha = alpha * length(c(mu[is.mkt, is.mkt])))

        # to recap, this number of animals
        # n.move3
        # have moved from market to wards.
        # this number of them
        # n.move2
        # moved among markets.
        # a total of 
        # n.move3
        # animals must move from wards to markets in move 1.
        # the number moved from each market is 
        # rowSums(move[is.mkt, !is.mkt])
        # however inter-market transfers mean that the numbers moving into
        # the markets must be adjusted
        
        # inflow to each market in move 2
        # colSums(move[is.mkt, is.mkt])
        # outflow from each market in move 2
        # rowSums(move[is.mkt, is.mkt])
        
        print("simulating ward-market movements")
        # to have a net flow of zero through the markets, the number put in from wards must be
        n.move1.tot <-
          rowSums(move[is.mkt, !is.mkt]) + rowSums(move[is.mkt, is.mkt]) - colSums(move[is.mkt, is.mkt])
        
        # due to the patchy nature of the data, it's possible for some of these input totals to be negative
        # which is not possible. this might not happen with real data, but if it does, re-draw until all >= 0
        if(all(n.move1.tot >= 0)) break
      }
    }
    
    # these can be drawn from a multinomial distribution
    if(sum(mu.z[!is.mkt, is.mkt]) > 0) {
      move[!is.mkt, is.mkt] <-
        sapply(rownames(mu.z)[is.mkt], 
               function(mkt) {
                 print(mkt)
                 lambda <- c(mu[!is.mkt, mkt])
                 rztmultinom.dir(n = 1, size = n.move1.tot[mkt], prob = lambda / sum(lambda), alpha = alpha)
               })
    } else {
      move[!is.mkt, is.mkt] <- 0
    }
    
    stopifnot(sum(move * move.z) == sum(move * move.z))
    move
  }

# function to simulate from zero-truncated negative binomial
rztnbinom.mu <-
  function(n, size = Inf, prob, mu) {
    require(actuar)
    if(missing(prob)) prob <- size / (size + mu)
    rztnbinom(n = n, size = size, prob = prob)
  }



# function to approximate a zero-truncated negative binomial using multinomial-dirichlet
# so that the total can be fixed
rztmultinom.dir <- 
  function(n, size, prob, alpha) {
    require(dirmult)
    # if the total is zero, the output elements must all be zero
    if(size == 0) return(rep(0, length(prob)))
    which.zero <- prob == 0
    prob.nz <- prob[!which.zero]
    # because this is a zero-truncated count distribution, all the elements with non-zero mu must generate 1s 
    # so total must be at least the number of non-zero elements
    # if not, output 1s for the highest probability elements, up to a total of "size"
    if(size <= length(prob.nz)) {
      out.nz <- rep(0, length(prob.nz))
      out.nz[(1 + length(prob.nz) - rank(prob.nz)) <= size] <- 1
    } else {
      # only need to simulate for the nonzero elements, so select these
      p <- rdirichlet(n = n, alpha = alpha * prob.nz / sum(prob.nz))
      print(paste("size", size))
      out.nz <- apply(p, 1, function(pp) 1 + rmultinom(n = 1, size = size - length(pp), prob = pp))[, 1]
    }
    out <- rep(NA, length(prob))
    out[which.zero] <- 0
    # replace nonzero elements with out.nz. first check they are the same length
    stopifnot(length(out[!which.zero]) == length(out.nz))
    out[!which.zero] <- out.nz
    out
  }


