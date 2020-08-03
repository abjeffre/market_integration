
k = 1,                                      # Learning efficency, note that this should be a vector of length ngoods, but single values can be inputed for parameter sweeps.  These are then transformed within the model.
alpha = 1,                                  # Production elasticity, note that this should be a vector of length ngoods, but single values can be inputed for parameter sweeps. These are then transformed within the model.
skew = 2,                                   # Shape on dirichlet distribution that determines initial time allocations, note that this should be a vector of length ngoods, but single values can be inputed for parameter sweeps. These are then transformed within the model
yeild = 100,                                # The max group yeild for each good, , note that this should be a vector of length ngoods, but single values can be inputed for parameter sweeps. These are then transformed within the mode. 
rank = rep(1, ngoods),                      # Elasticies/rank importance of each good in utility functions.  NOTE A: The this parameter is contrained by the choice of utility function and number of goods. Note B. This should be a vector of length ngoods, but single values can be inputed for parameter sweeps, these are then transformed within the model
max_visits = 20,                            # Number of trading visits before individuals stop.
relatedness = 0,                            # Relatedness within the group,  1 is everyone is a clone.
trade_cost = 0,                             # The cost of each trading trip.
utility = "leontif",                        # Currently the simulation accepts Cobb douglas, Log-Linear, CES, Stone-Greary and Leontif, please be aware of constraints imposed by the utility function on "Rank"
gamma = 1,                                  # This is the subsitence level for Stone Greary utility functions: Note that this value is scaled 1/ngoods  
rho = .5,                                   # This controls the subsitution in CES functions. 
beta = 1,                                   # This controls how returns are a function of total group effort.
omega = 0.1,                                # Probability of environmental change.
trading = T,                                # This is a meta variable that toggles if the trading portion of the model runs.
sharing = T,                                # This is a meta variable that toggles if the trading portion of the model runs.
specialization = T,                         # This is a meta variable that toggles if the specialization portion of the model runs.
fluctuations = F                            # This is a meta variable that toggles if there is environmental variation.

){
  
  # Make useful functions
  softmax <- function(x) {exp(x)/sum(exp(x))}
  
  # Make output
  output <- list()
  extinction <- list()
  
  ##################################
  ##### Build full parameters ######
  # Note that if you want to have specific yeilds, skews or requried resources that vary in accord with a given good, then you must specifiy each manually.
  
  for(simn in 1:nsim){
    #check basic inputs
    if(length(rank)==1) rank <- rep(rank, ngoods)
    
    #check constraints on utility functions
    if(utility == "stone"){
      if(rank(sum(rank)!=1)){
        warning("Rank must sum to 1 for Stone Greary.  Applying softmax transformation " ,immediate = T)
        rank <- softmax(rank)
      }
    }
    
    if(utility == "ces"){
      sigma <- 1/(1-rho)  
      if(rank(sum(rank)!=1)){
        warning("Rank must sum to 1 for CES.  Applying softmax transformation " ,immediate = T)
        rank <- softmax(rank)
      }
    }
    
    
    if(length(yeild) == 1)yeild <- rep (yeild/ngoods, ngoods)  #NOTE the yeild is stable regardless of total goods.
    if(length(beta) == 1)beta <- rep (beta, ngoods)  #NOTE the yeild is stable regardless of total goods.
    if(length(skew) == 1)skew = rep(skew, ngoods)
    if(length(gamma) == 1)gamma = rep(gamma/ngoods ,ngoods)
    
    #Initalize Arrays
    p <- array(0, c(n, ngoods, ngen)) #proption of lifetime effort dedicated to producing each good. 
    prices <- array(0, c(n, ngoods, ngen)) # Price 
    learn <- array(0, c(n, ngoods, ngen)) # Time spent learning
    do <- matrix(0, nrow = n, ncol = ngoods) # Time spent producing
    produce <- matrix(0, nrow = n, ncol = ngoods) #goods produced - unscaled by yeild 
    inv <- array(0, c(n, ngoods, ngen)) #inventory of goods produced after scaled by total group production
    pre_sale_inv <- array(0, c(n, ngoods, ngen)) # Pre-sale inventory - mostly used for debugging
    optima <- array(0, c(n, ngoods, ngen)) # Optimal amount of good to be aquired as per marshallian demand
    need <- array(0, c(n, ngoods, ngen)) # Amount of each good needed to be purchased or sold -   n = (inventory- optima)
    sell <-  array(0, c(n, ngoods, ngen)) # A list of sales
    buy <-  array(0, c(n, ngoods, ngen)) # A list of purcahses
    visits <- matrix(0, nrow =n, ncol = ngen) # Number of visits in each time step
    
    
    #set up names in the agent array  
    columns <- c("gid","id", "p", "trade_type", "share_type", "utility")
    agents <- array(0, c(n, length(columns), ngen))
    colnames(agents) <- columns
    
    #give id numbers
    agents[,"id",] <- seq(1:n)
    
    #Set individuals into groups NOTE: Individuals stay in the same groups as their parents 
    agents[, "gid", 1] <- as.integer(cut(agents[,"id",1], groups))
    
    #Set up initial time allocations
    p[,, 1] <- rdirichlet(n, skew)
    
    #Set up share_types 1 type one does not share and type 2 does
    agents[, "share_type", 1] <- sample(1:2, n, replace = T)
    
    #Set up trade_types 1 type one does not trade and type 2 does
    agents[, "trade_type", 1] <- sample(1:2, n, replace = T)
    
    #Set up price vectors
    for (g in 1:ngoods){
      prices[, g, 1] <- runif(n, 0, 10)
    }
    
    
    #########################################
    ########### Begining of Time ############
    
    # loop over generations
    for ( t in 1:(ngen) ) {
      
      # If specialization is false the spec part of the model runs to determine production, but individuals simply get a new p vectors each round to ensure that 
      # evolution of specialization does not happen. 
      if(specialization == F)  p[,, t] <- rdirichlet(n, skew)
      
      # Environmental Change, not yet implemented fully
      # For initial testing purposes we shall begin with introducing production shocks on only one good
      
      if(fluctuations == T){
        if(runif(1, 0, 1) < omega){
          if(beta[1] ==1) ifelse(rbinom(1, 1, .5)==1, 2, .5) 
          if(beta[1] == 0.5) ifelse(rbinom(1, 1, .5)==1, 2, 1)
          if(beta[1] == 0.5) ifelse(rbinom(1, 1, .5)==1, 2, 1)
        }
      }
      #Determine the amount of time spent learning and doing, NOTE: ALL equations come from the feldman paper    
      
      for (i in 1:n){
        for (g in 1:ngoods){
          
          learn[i,g,t] <- ifelse(p[i,g,t] <= 1/(alpha*k),   0,    (alpha*k*p[i,g,t]-1)/((alpha +1)*k) )
          do[i,g] <- ifelse(p[i,g,t] <= 1/(alpha*k),    p[i,g,t],    (1/(alpha*k))*((alpha*(k*p[i,g,t]+1))/((alpha +1)))^(alpha+1) )
          produce[i,g] <- do[i,g]*(1+k*learn[i,g,t])^alpha
          
        }#end ngoods
      }#end time allocations  
      
      #Calculate Production
      for (i in 1:n) {
        for(g in 1:ngoods) {
          
          #scale production by number of producers and total yeild of resource
          # This is the old inv based on feldman model 1
          
          # inv[i,g,t] <- (produce[i,g,t]/sum(produce[,g,t][agents[,"gid",t] == agents[i,"gid",t]]))*yeild[g]
          
          
          
          #This new inv is defined as Yg|group =beta*E(produce|group)/(1+E(produce|group))
          # Based on Feldman model 2. 
          #This ensures that production is an increasing  concave downward function of production 
          # to see how the function behaves use code below
          #  seq <- seq(1, 10, length.out = 100)
          #  dat <- sapply(seq, function(x) x*beta/(1+x))
          #  plot(dat ~ c(1:100), type = "l")
          #  if we log 
          #  dat <- sapply(seq, function(x) log(x*beta)/(1+log(x)))
          #  plot(dat ~ c(1:100), type = "l")
          
          inv[i,g,t] <- (produce[i,g]/sum(produce[,g][agents[,"gid",t] == agents[i,"gid",t]]))*(mean(produce[,g][agents[,"gid",t] == agents[i,"gid",t]])*beta[g])/
            (1+mean(produce[,g][agents[,"gid",t] == agents[i,"gid",t]]))*yeild[g]
          
          #correct na's and nans for dividing by zero
          inv[i,g,t][is.nan(produce[i,g])] <- 0
          inv[i,g,t][is.na(produce[i,g])] <- 0
          inv[i,g,t][is.nan(inv[i,g,t])] <- 0
          inv[i,g,t][is.na(inv[i,g,t])] <- 0
        }#end goods
        if(all(inv[i,,t] == 0)) warning("Individuals Producing zero total goods, check inv")
      }#end n
      
      
      ################################################################
      ################## Marshallian Demand Functions ################
      
      # Individuals calculate their demand functions given their utility functions.
      # All demand functions are calculated using marshallian demand.
      # Marshallian demand means that individuals solve the constrained optimization problem given their utility function and their prices and derive the optimal
      ## bundle of goods that they would like to consume. 
      # This section contains the marshallian demand functions for 5 diffrent utility functions. 
      # Only one utility function is active in any given simulation
      
      #This toggles whether or not the trading portion of the model runs       
      if(trading == T){    
        
        
        
        # The leontif utility function is a function where the goods in the consumption bundle are perfect compliments. 
        
        if(utility == "leontif"){
          for (i in 1:n){
            #if(agents[i,"trade_type",t ] ==2){
            lambda <- sum(inv[i,,t]*prices[i,,t])/sum(rank*prices[i,,t])
            
            # Calculate thier optimum consmuption basket
            
            for(g in 1:ngoods){
              optima[i,g,t] <- lambda*rank[g]
            }
            
            #Calculate the "need" for goods to be bought and sold - positive is the amount of good x that the indiviudal must aquire, negative amount to sell
            
            for (g in 1:ngoods){
              need[i, g, t] <- optima[i,g,t] - inv[i,g,t]
              
              #Round off the need in order solve the market equations
              need[i, g, t]  <- round(need[i, g, t], 3)
              
            }#setup needs
          }#ends indv
        }#ends leontif
        
        
        # Cobb douglas and loglinear utility functions are maximized using the marshallian demand functions below. 
        
        temp <- c(0)
        if(utility == "cobb") temp <-1 #solves problem with | not working on char
        if(utility == "linear") temp <- 1 
        if(temp == 1){
          for(i in 1:n){ 
            b <-sum(inv[i,,t]*prices[i,,t])
            for(g in 1:ngoods){
              
              optima[i,g,t] <- (rank[g]/sum(rank))*(b/prices[i,g,t])
            }
            #Calculate the "need" for goods to be bought and sold - positive is the number of goods to aquire, negative is the number to sell
            
            for (g in 1:ngoods){
              need[i, g, t] <- optima[i,g,t] - inv[i,g,t]
              
              #Round off the need in order solve the market equations
              need[i, g, t]  <- round(need[i, g, t], 3)
            }#end goods
          }#end indv
        }#end cobb douglas/loglinear
        
        
        # The CES function allows for us to specify a constant rate subsitution across diffrent goods. 
        # Subsitution parameter is specified by Rho. 
        
        if (utility == "ces") {
          for (i in 1:n){
            b <-sum(inv[i,,t]*prices[i,,t])
            for(g in 1:ngoods){
              
              optima[i,g,t] <- (((rank[g]^sigma)*prices[i,g,t]^(1-sigma))/sum(rank^sigma*prices[i,,t]^(1-sigma))) *(b/prices[i,g,t])
              
            }#Calculate the "need" for goods to be bought and sold - positive is the number of goods to aquire, negative is the number to sell
            
            for (g in 1:ngoods){
              need[i, g, t] <- optima[i,g,t] - inv[i,g,t]
              
              #Round off the need in order solve the market equations
              need[i, g, t]  <- round(need[i, g, t], 3)
            } #ends goods
          }#ends invididuals
        }#ends ces
        
        
        # The Stone Greary Utility function is a derivate of the Cobb douglas but has a baseline subsistence component that specifies the minimum amount of each good that must be consumed
        # Note the gamma parameter is needed for stone greary
        if (utility == "stone") { 
          for(i in 1:n){
            b <-sum(inv[i,,t]*prices[i,,t])  
            
            for(g in 1:ngoods){
              
              
              optima[i,g,t] <- gamma[g]+(rank[g]/prices[i,g,t])*(b-sum( prices[i,,t]*gamma))
            }#end goods
            
            for (g in 1:ngoods){
              need[i, g, t] <- optima[i,g,t] - inv[i,g,t]
              
              #Round off the need in order solve the market equations
              need[i, g, t]  <- round(need[i, g, t], 3)
            }
          }#end indv
        }#end stone
        
        
        ###############################################
        ############# The Market Place ################
        
        # Randomize buyer order
        # For each buy they randomize visitation order
        # Individuals visit a number of people upto a maximum specified by max visits
        # For each person they visit they randomize the order of goods they are trying to buy
        # They calculate the amount of good x they need to purcahse
        # They then construct a basket of goods they are willing to sell that is equal in value (according to their own price vector) that is equal to the value of goods 
        #        they are trying to buy (again according to their own price vector).
        # The reciever then evaluates value of the goods according to their own price vector and either rejects or accepts the trade.
        # If the trade takes place goods are exchanged and inventories updated. 
        # The individual then moves on to the next good that they "Need"
        # Individuals stop trading if they do not need any more of any good in their consumption basket or if they reach max visits.
        
        
        
        pre_sale_inv[,,t]<-inv[,,t]
        
        #Set up the markets    
        #Randomize bidder order
        bid_order <- sample(1:n, n, replace = FALSE)
        
        
        
        
        for (i in bid_order) {
          if(agents[i,"trade_type",t ] ==2){
            visit <- 1   
            #randomize reciever order
            recieve_order <- sample(1:n, n, replace = FALSE)
            
            for(j in recieve_order){
              #Check to see if goods need exists and max visits not reached
              if(visit < max_visits & !all(need[i,,t]>=0)){
                if(agents[j, "trade_type",t]==2){
                  if(i != j){
                    #randomize buying order order
                    buy_order <- sample(1:ngoods, ngoods, replace = FALSE)
                    
                    for(g in buy_order){
                      proposal_sell_quant <- c(0)
                      proposal_buy_quant <- c(0)
                      # Bidder evaluates whether they need to purchase good X 
                      # NOTE: negative need is equvelent to needing to sell. 
                      if(need[i, g, t] > 0 & need[j, g, t] < 0){
                        #Reciver proposes max quantity for sale
                        maxAvailibe <- -1*need[j, g, t]
                        if(maxAvailibe >= need[i, g, t]) {  #the reciever evaluates their stock and sees if they can sell to the bidders total demand
                          proposal_buy_quant <- need[i, g, t] #the bidder tries to buy the max needed to satiate their desires. 
                          # The offer:
                          #the bidder chooses an amount of alternative good equivenlent to the amount they want to buy. They loop through their inventroy choosing each good they dont need and adding it to the basket until they reach what they think the price is.
                          proposal_sell_quant <- rep(0, ngoods)
                          
                          for (u in 1:ngoods){
                            
                            if(u !=g){
                              
                              if (need[i, u, t] < 0 & proposal_buy_quant*prices[i,g,t] >= sum(proposal_sell_quant*prices[i,,t])){
                                
                                mental_account <- need[i,u,t]  
                                
                                while(mental_account <=0 & proposal_buy_quant*prices[i,g,t]  >= sum(proposal_sell_quant*prices[i,,t])){
                                  
                                  proposal_sell_quant[u] <- proposal_sell_quant[u] +0.001  
                                  
                                  mental_account <- mental_account+0.001
                                  
                                }
                                
                              }#end while loop on good
                              
                            }#ran out of good
                            
                          }#check need on good
                          
                          
                          # The reciver insepects the value
                          
                          if(sum(proposal_sell_quant*prices[j, , t]) >= (proposal_buy_quant*prices[j, g, t])) {
                            
                            #Record Transaction
                            sell[j,g,t] <- proposal_buy_quant + sell[j,g,t]
                            buy[j, , t] <- proposal_sell_quant + buy[j,,t]
                            
                            buy[i,g,t] <- proposal_buy_quant + sell[i,g,t]
                            sell[i, , t] <- proposal_sell_quant + buy[i,,t]
                            
                            
                            #Transfer goods
                            
                            inv[j,,t] <- inv[j,,t] + proposal_sell_quant # seller gets goods
                            inv[j,g,t] <- inv[j,g,t] - proposal_buy_quant # seller looses goods
                            
                            inv[i,g,t] <- inv[i,g,t] + proposal_buy_quant # buyer gets goods
                            inv[i,,t] <- inv[i,,t] - proposal_sell_quant # buyer looses goods
                            
                            #update needs and desires
                            need[j,,t] <- need[j,,t] - proposal_sell_quant # seller desires purchased goods less
                            need[j,g,t] <- need[j,g,t] + proposal_buy_quant # seller may for whatever reason want these goods more
                            
                            need[i,,t] <- need[i,,t] + proposal_sell_quant # buyer desires sold goods more
                            need[i,g,t] <- need[i,g,t] - proposal_buy_quant #  buy desires purchased goods less
                            
                          }# trade eval
                          
                        }else{                 # try to find an amount the seller is willing to sell
                          
                          proposal_buy_quant <- need[j, g, t]*-1 # the amount traded is the max of the seller because their inventory is too small for i
                          
                          proposal_buy_value <- prices[i, g, t]*proposal_buy_quant # bidder evaluates the amount and proposes an amount of the gathered good. 
                          
                          proposal_sell_quant <- rep(0, ngoods)
                          for (u in 1:ngoods){
                            if(u !=g){
                              if (need[i, u, t] < 0 & proposal_buy_value >= sum(proposal_sell_quant*prices[i,,t])){
                                mental_account <- need[i,u,t]  
                                while(mental_account <=0 & proposal_buy_value >= sum(proposal_sell_quant*prices[i,,t])){
                                  proposal_sell_quant[u] <- proposal_sell_quant[u] +0.001  
                                  mental_account <- mental_account+0.001
                                }
                              }#end while loop on good
                            }#ran out of good
                          }#check need on good
                          
                          # The reciver insepects the value
                          if(sum(proposal_sell_quant*prices[j, , t]) >= (proposal_buy_quant*prices[j, g, t])) {
                            
                            #Record Transaction
                            sell[j,g,t] <- proposal_buy_quant + sell[j,g,t]
                            buy[j, , t] <- proposal_sell_quant + buy[j,,t]
                            
                            buy[i,g,t] <- proposal_buy_quant + sell[i,g,t]
                            sell[i, , t] <- proposal_sell_quant + buy[i,,t]
                            
                            
                            #Transfer goods
                            
                            inv[j,,t] <- inv[j,,t] + proposal_sell_quant # seller gets goods
                            inv[j,g,t] <- inv[j,g,t] - proposal_buy_quant # seller looses goods
                            
                            inv[i,g,t] <- inv[i,g,t] + proposal_buy_quant # buyer gets goods
                            inv[i,,t] <- inv[i,,t] - proposal_sell_quant # buyer looses goods
                            
                            #update needs and desires
                            need[j,,t] <- need[j,,t] - proposal_sell_quant # seller desires purchased goods less
                            need[j,g,t] <- need[j,g,t] + proposal_buy_quant # seller may for whatever reason want these goods more
                            
                            need[i,,t] <- need[i,,t] + proposal_sell_quant # buyer desires sold goods more
                            need[i,g,t] <- need[i,g,t] - proposal_buy_quant #  buy desires purchased goods less
                            
                            
                          } # accept or reject trade
                          
                        }#end lowstock transfer
                      }#end check need for good
                    } # end good
                  } #end if not self 
                  visit <- visit + 1
                }#end trader condition
                
              } #end j loop
            } #end while loop
            visits[i,t] <- visit 
          } #end agents type condition
        } #end i loop
        
      }# end trading
      
      
      ######################
      ###### Sharing #######
      
      # Individuals in each group can share with members of their own group
      # If they do share they put all of their goods into a local "pot"
      # The pot is then subdivided amongst all individuals in the group regardless of contribution
      
      
      if(sharing == TRUE){ # this toggles whether the sharing component of the model is active     
        group_pot <- matrix(0, nrow = groups, ncol = ngoods)
        
        #Indiviuals who are sharers will contribute to a group pot
        for (i in 1:n){
          if(agents[i,"share_type", t]==2){
            for(g in 1:ngoods){
              if (inv[i,g,t] >= 0) group_pot[agents[i,"gid",t], g] <- group_pot[agents[i,"gid",t],g] + inv[i,g,t]
              inv[i,g,t] <- 0
            }
          }    
        }  
        
        #divide the pot
        
        for(i in 1:n){
          inv[i,,t] <- inv[i,, t] + group_pot[agents[i, "gid", t],]/sum(agents[, "gid", t]==agents[i, "gid", t])
        }
        
      } # End sharing
      
      
      #########################   
      ####CALCULATE UTILITY####
      
      for (i in 1:n){
        
        if(utility == "stone") agents [i, "utility", t] <- ifelse(!is.nan(base::prod((inv[i,,t]-gamma)^rank) - trade_cost*visits[i,t]), base::prod((inv[i,,t]-gamma)^rank) - trade_cost*visits[i,t], 0)  #note the ifelse basically kills individuals with no utility 
        if(utility == "ces") agents[i, "utility", t] <-  sum(rank*inv[i,,t]^rho)^(1/rho) - trade_cost*visits[i,t]
        if(utility == "cobb") agents[i,"utility",t] <- base::prod(inv[i,,t]^rank) - trade_cost*visits[i,t]
        if(utility == "loglinear") agents[i,"utility",t] <- sum(log(inv[i,,t]+1)*rank) - trade_cost*visits[i,t]
        if(utility == "leontif")agents[i,"utility",t] <- min(inv[i,,t]/rank) - trade_cost*visits[i,t]
        
        #This is here incase we get negative utility - it allows for the sample function to opperate properly
        if(any(agents[,"utility", t] <0)) warning("Negative utility encountered. Converted to zero")
        agents[,"utility", t][agents[,"utility", t] < 0] <- 0
      }
      
      
      
      ##########################      
      ### CONSUME INVENTORY ####
      # IGNORE FOR NOW
      # Note: This section is not fully opperational, it is designed for once we add in multiple stages of trade per lifetime or if we assume intergenerational wealth 
      ## transmission mechanism 
      #Note: this is only necessary for once we add in some sort of transmission mechanism
      #Note: Must add in for linear utlity function
      
      if(utility == "leontif") inv[,, t] <- inv[,, t] - agents[,"utility",t]
      
      ###########################################   
      ### CALCULATE UTILITY WITH RELATEDNESS #### 
      
      # for(i in 1:n){
      # agents[i, "utility", t] <- agents[i, "utility", t] + relatedness*sum(agents[, "utility", t]==agents[i, "gid", t])    
      # }
      
      ############################################
      ##### Evolutionary Dynamics ################     
      
      # First check to see if we get extinction
      if(all(agents[,"utility",t] <= 0)) {
        extinction[simn] <- paste0("Extinction event at generation ", t)
        output[[simn]] <- list(agents = agents,
                               p = p,
                               prices = prices,
                               learn = learn,
                               need = need,
                               inv = inv,
                               sell = sell,
                               buy = buy,
                               visits = visits,
                               extinction = extinction)
        break
      }
      
      
      # Replication and Transmission  
      # Babies are born   
      # Babies get parents who the inheret price vectors, production values, trading types, and sharing types
      
      
      pid <- sample(agents[,"id",t] , size=n , replace=TRUE , prob=agents[,"utility",t])
      babies_prices <- prices[pid, , t] 
      babies_p <- p[pid,, t]
      babies_share_type <- agents[pid, "share_type", t]
      babies_trade_type <- agents[pid, "trade_type", t]
      
      
      #Mutation
      for ( i in 1:n) {
        #Note: The multiplier on the dirichlet controls the varainace of the mutation
        if ( runif(1) < mu ) temp <- babies_p[i, ] <- rdirichlet(1, babies_p[i, ]*100)
        
        if ( runif(1) < mu ) babies_prices[i,] <- ifelse(rbinom(1, 1, .5)==1, babies_prices[i,]*.95, babies_prices[i,]*1.05)
        
        if ( runif(1) < mu ) babies_share_type[i] <- ifelse(babies_share_type[i]==1, 2, 1) 
        
        if ( runif(1) < mu ) babies_trade_type[i] <- ifelse(babies_trade_type[i]==1, 2, 1) 
        
      }#end mutation
      
      #   Babies become adults
      
      if(t< ngen){
        p[,,t+1] <- babies_p
        prices[,,t+1] <- babies_prices
        agents[,"gid",t+1] <- agents[,"gid",t]
        agents[,"share_type",t+1] <- babies_share_type
        agents[,"trade_type",t+1] <- babies_trade_type
      }
    }#end time
    output[[simn]] <- list(agents = agents,
                           p = p,
                           prices = prices,
                           learn = learn,
                           need = need,
                           inv = inv,
                           sell = sell,
                           buy = buy,
                           visits = visits)
  }#end sims
  return(output)  
  
  
}#end function


###############################
########PARAMETER SWEEP #######


seq_trade <- expand.grid(
  nsim = 1,                                                                    # Simulations Per Parameter
  ngen = 1e3,                                                                  # Generations
  n = 200,                                                                     # Population Size
  groups = c(100, 10, 2),                                                      # Groups
  ngoods = c(2, 5, 10),                                                        # Number of Goods
  mu = 0.001,			                                                             # Mutation rate
  k = c(1),                                                                    # Learning efficency
  alpha = c(1),                                                                # Production elasticity
  skew = 2,                                                                    # Shape on dirichlet distribution
  yeild = 100,                                                                 # The max group yeild for each good 
  rank = 1,                                                                    # The relative importance of each good NOTE: This value depends on the choice of utility function and number of goods
  max_visits = 20,                                                             # Number of trading visits before they give up
  relatedness = c(0, .125, .25),                                               # Relatedness within the group,  1 is everyone is a clone 
  trade_cost = c(0, .1, 1),                                                     # The cost of each trading trip I suspect a value of .25 is a reasonable cost to begin with
  utility = c("leontif", "cobb", "loglinear", "stone", "ces"),                 # Currently the simulation accepts Cobb douglas, Log-Linear, CES, Stone-Greary and Leontiff
  gamma = 1,                                                                   #this is the subsitence level for stone greary utility functions  
  rho = .5,
  beta = 1,
  omega = 0,
  trading = TRUE,
  sharing = FALSE,
  specialization = FALSE,
  fluctuations = FALSE
)


############################################################################################
##################### PARALELL PARAMETER  SWEEP FOR TRADE ##################################

library(parallel)


for (i in 1:9){
if (i > 1) {
  result_dol <- mclapply(
    ((i-1)*45+1):(i*45),
    function(i) sim_trade_spec(seq_trade$nsim[i], seq_trade$ngen[i], seq_trade$n[i], seq_trade$groups[i], seq_trade$ngoods[i],
                               seq_trade$mu[i], seq_trade$k[i], seq_trade$alpha[i], seq_trade$skew[i],
                               seq_trade$yeild[i], seq_trade$rank[i], seq_trade$max_visits[i], seq_trade$relatedness[i],
                               seq_trade$trade_cost[i], seq_trade$utility[i], seq_trade$gamma[i], seq_trade$rho[i], 
                               seq_trade$beta[i],  seq_trade$omega[i], seq_trade$trading[i], seq_trade$sharing[i],
                               seq_trade$specialization[i], seq_trade$fluctuations[i]),
    mc.cores=45)
  
  saveRDS(result_dol, file = paste0("market_integration/ABM_spec_trade/sweeps/trade/trade_sweep_",((i-1)*45+1),"_", (i*45),".rds"))
}
  
if(i == 1){
  
result_dol <- mclapply(
  1:45 ,
  function(i) sim_trade_spec(seq_trade$nsim[i], seq_trade$ngen[i], seq_trade$n[i], seq_trade$groups[i], seq_trade$ngoods[i],
                             seq_trade$mu[i], seq_trade$k[i], seq_trade$alpha[i], seq_trade$skew[i],
                             seq_trade$yeild[i], seq_trade$rank[i], seq_trade$max_visits[i], seq_trade$relatedness[i],
                             seq_trade$trade_cost[i], seq_trade$utility[i], seq_trade$gamma[i], seq_trade$rho[i], 
                             seq_trade$beta[i],  seq_trade$omega[i], seq_trade$trading[i], seq_trade$sharing[i],
                             seq_trade$specialization[i], seq_trade$fluctuations[i]),
  mc.cores=45)

saveRDS(result_dol, file = "market_integration/ABM_spec_trade/sweeps/trade/trade_sweep_1_45.rds")
}
}


## For these models I should first extract the change in trader type across all diffrent parameter combos

means <- matrix(NA, nrow = 1000, ncol = 45)
temp <- readRDS("market_integration/ABM_spec_trade/sweeps/trade/trade_sweep_1_45.rds")
prices <- array(NA, c(1000, 10, 45)) 
for(j in 1:45){
  agents <- temp[[j]][[1]][[1]]
  price <- temp[[j]][[1]][["prices"]]
  for(i in 1:1000){
  means[i,j] <- mean(agents[,"trade_type",i])
  prices[i,,j] <- apply(price[,,i], 2, mean)
  }  
  
}

for (i in 1:45) plot(1:1000, means[,i], type = "l", main = paste0("ROW:", i,"trade cost: ", seq_trade$trade_cost[i], " ngoods: ", seq_trade$ngoods[i], 
                                                                  " relatedness: ", seq_trade$relatedness[i], "ngroups: ", seq_trade$groups [i]))
for(i in 1:45){ plot(1:1000, prices[,1,i], type = "l", ylim = c(0, 10), main = paste0("ROW:", i, "trade cost: ", seq_trade$trade_cost[i], " ngoods: ", seq_trade$ngoods[i], 
                                                                                      " relatedness: ", seq_trade$relatedness[i], " ngroups: ", seq_trade$groups [i]))
for(j in 2:10) {
  lines(1:1000, prices[,j,i])
}
}




means <- rep(NA, 45)
for(i in 1:length(sim_trade_sweep45)){
  agents <- sim_trade_sweep45[[i]][[1]][[1]]
  means[i] <- mean(agents[,"trade_type",1000])
}  














############################################################################################################
####################### SERIAL PARAMETER GENERIC ###########################################################

count <- 1
results <- list()
for(i in 1:nrow(seq_trade)){
  
  results[[count]] <- sim_dol_ngoods(seq_trade$nsim[i], seq_trade$ngen[i], seq_trade$n[i], seq_trade$groups[i], seq_trade$ngoods[i],
                                     seq_trade$mu[i], seq_trade$k[i], seq_trade$alpha[i], seq_trade$skew[i],
                                     seq_trade$yeild[i], seq_trade$rank[i], seq_trade$max_visits[i], seq_trade$relatedness[i],
                                     seq_trade$trade_cost[i], seq_trade$utility[i], seq_trade$gamma[i], seq_trade$rho[i], 
                                     seq_trade$beta[i],  seq_trade$omega[i], seq_trade$trading[i], seq_trade$sharing[i],
                                     seq_trade$specialization[i], seq_trade$fluctuations[i])
  print(count)
  count <- count+1  
}







sim_dol_ngoods(
  ngen = 1000, 
  n= 200, 
  groups = 20,
  ngoods = 10, 
  mu =0.001,
  k =10,
  alpha = 10,
  skew = 2,
  yeild = 100,
  rank = 1,
  max_visits = 20,
  relatedness = 0.125,
  trade_cost = 0,
  utility = "cobb",
  gamma = 1,
  rho = 1,
  beta = 1)



#broke at i = 297 j = 17
#check and see why this broke 
ngen = 1000 
n= 200 
groups = 2
ngoods = 10 
mu =0.001
k =10
alpha = 10
skew = 2
yeild = 100
rank = 1
max_visits = 20
relatedness = 0
trade_cost = 0
utility = "leontif"
gamma = 1
rho = 0.5
beta = 1
omega = 0
trading = T
sharing = F
specialization = T
fluctuations = F


















library(parallel)


count <- 1
sim_results <- list()
for(i in 1:nrow(seq_spec_trade)){
  sim_results[[count]] <- sim_dol_ngoods(sim_dol_ngoods(seq_spec_trade$nsim[i], seq_spec_trade$ngen[i], seq_spec_trade$n[i], seq_spec_trade$groups[i], seq_spec_trade$ngoods[i],
                                                        seq_spec_trade$mu[i], seq_spec_trade$k[i], seq_spec_trade$alpha[i], seq_spec_trade$skew[i],
                                                        seq_spec_trade$yeild[i], seq_spec_trade$rank[i], seq_spec_trade$max_visits[i], seq_spec_trade$relatedness[i],
                                                        seq_spec_trade$trade_cost[i], seq_spec_trade$utility[i], seq_spec_trade$gamma[i], seq_spec_trade$rho[i], 
                                                        seq_spec_trade$beta[i], seq_spec_trade$beta[i], seq_spec_trade$trading[i], seq_spec_trade$sharing[i], seq_spec_trade$specialization[i]))
  count <- count +1
  print(count)
}

saveRDS(sim_results, file = "sim_results")



sim_dol_ngoods(
  ngen = 1000, 
  n= 200, 
  groups = 20,
  ngoods = 10, 
  mu =0.001,
  k =10,
  alpha = 10,
  skew = 2,
  yeild = 100,
  rank = 1,
  max_visits = 20,
  relatedness = 0.125,
  trade_cost = 0,
  utility = "cobb",
  gamma = 1,
  rho = 1,
  beta = 1)

#############################
########Plotting ############





sim.plot <- function(df, g){
  m <- matrix(ncol = dim(df)[1], nrow =dim(df)[3])
  for (i in 1:dim(df)[3]){
    m[i,] <- as.integer(cut(df[, g, i], c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)))
  } 
  collist <- c( "yellow","gold", "orange","red","green",  "blue" ,"purple", "darkgrey","brown", "black")
  plot(NULL, xlim = c(0, dim(df)[3]), ylim = c(0, 1))
  for(i in 1:10){
    means <- apply(m==i, 1, mean)
    lines(1:ngen, means, col =collist[i] )
  }
}




#plot types
plot(NULL, xlim =c(0, 1000), ylim = c(1, 2))
means <- apply(temp$agents[, "trade_type", ], 2, mean)
lines(1:1000, means, col = "green")    
means <- apply(temp$agents[, "share_type", ], 2, mean)
lines(1:1000, means, col = "red")    

#plot prices
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
for(i in 1:dim(temp$prices)[2]){
  means <- apply(temp$prices[, i, ], 2, mean)
  lines(1:1000, means, col = collist[i])
}



#plot demand
plot(NULL, xlim =c(0, 1000), ylim = c(-5, 5))
for (i in 1:dim(temp$prices)[2]){
  means <- apply(temp$optima[,i,]-temp$pre_sale_inv[,i,], 2, mean)
  lines(1:1000, means, col = collist[i])    
}


#price to demand ratio
plot(NULL, xlim =c(0, 1000), ylim = c(0, 15))
for (i in 1:dim(temp$prices)[2]){
  means <- apply(temp$optima[,2,1000]-temp$pre_sale_inv[,2,1000]/temp$prices[,2,], 2, mean)
  lines(1:1000, means, col = collist[i])    
}



#plot prices variance
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
for(i in 1:dim(temp$prices)[2]){
  means <- apply(temp$prices[, i, ], 2, sd)
  lines(1:1000, means, col = collist[i])
}

#plot trades Proportion of trades 
plot(NULL, xlim =c(0, 1000), ylim = c(0, 1))
for(i in 1:dim(temp$prices)[2]){
  means <- apply(temp$sell[, i, ]>0, 2, mean)
  lines(1:1000, means, col = collist[i])
}


#plot profit
plot(NULL, xlim =c(0, 1000), ylim = c(-10, 20))
for(i in 1:dim(temp$prices)[2]){
  means <- apply(temp$sell[, i, ]*apply(temp$price[, i, ],2,mean), 2, mean)-apply(temp$buy[, i, ]*apply(temp$prices[, i, ],2,mean), 2,mean)
  lines(1:1000, means, col = collist[i])
}


#plot visits
plot(NULL, xlim =c(0, 1000), ylim = c(0, 20))
means <- apply(temp$visits, 2, mean)
lines(1:1000, means, col = "green")    


#plot inequality
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
means <- apply(temp$agents[, "utility", ], 2, sd)
lines(1:1000, means, col = "green")    


#plot need
plot(NULL, xlim =c(0, 1000), ylim = c(-10, 10))
for(i in 1:dim(temp$prices)[2]){
  means <- apply(abs(temp$need[, i, ]), 2, mean)
  lines(1:1000, means, col = collist[i])    
}



#plot need across difrent types
mean_traders <- rep(0, ngen)
mean_non <- rep(0, ngen)
for (t in 1:ngen){
  mean_traders[t] <- mean(abs(temp$need[, , t][temp$agents[, "trade_type", t]==2]))
  mean_non[t] <- mean(abs(temp$need[, , t][temp$agents[, "trade_type", t]==1]))
}
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
lines(1:1000, mean_traders, col = collist[1])
lines(1:1000, mean_non, col = collist[2])


#plot learning
plot(NULL, xlim =c(0, 1000), ylim = c(0, 1))
for(i in 1:dim(temp$prices)[2]){
  means <- apply(temp$learn[, i, ], 2, mean)
  lines(1:1000, means, col = collist[i])    
}

#plot utility
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
means <- apply(temp$agents[, "utility", ], 2, mean)
lines(1:1000, means, col = collist)    


#plot utility across difrent types
mean_traders <- rep(0, ngen)
mean_non <- rep(0, ngen)
for (t in 1:ngen){
  mean_traders[t] <- mean(temp$agents[, "utility", t][temp$agents[, "trade_type", t]==2])
  mean_non[t] <- mean(temp$agents[, "utility", t][temp$agents[, "trade_type", t]==1])
}
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
lines(1:1000, mean_traders, col = collist[1])
lines(1:1000, mean_non, col = collist[2])




#plot utility across difrent producer types
mean_traders <- rep(0, ngen)
mean_non <- rep(0, ngen)
for (t in 1:ngen){
  mean_traders[t] <- mean(temp$agents[, "utility", t][temp$p[, 1, t] <.5])
  mean_non[t] <- mean(temp$agents[, "utility", t][temp$p[, 2, t] <.5])
}
plot(NULL, xlim =c(0, 1000), ylim = c(0, 10))
lines(1:1000, mean_traders, col = collist[1])
lines(1:1000, mean_non, col = collist[2])




#######################################
##TEST 


#plot utility across difrent producer types
mean_traders <- rep(0, ngen)
mean_non <- rep(0, ngen)
for (t in 1:ngen){
  mean_traders[t] <-mean((temp$optima[,2,t][temp$p[, 1, t] >.5] -temp$pre_sale_inv[,2,t][temp$p[, 1, t] >.5])-temp$need[, 2, t][temp$p[, 1, t] >.5])
  mean_non[t] <- mean((temp$optima[,2,t][temp$p[, 2, t] >.5] -temp$pre_sale_inv[,2,t][temp$p[, 2, t] >.5])-temp$need[, 2, t][temp$p[, 2, t] >.5])
}
plot(NULL, xlim =c(0, 1000), ylim = c(-10, 10))
lines(1:1000, mean_traders, col = collist[1])
lines(1:1000, mean_non, col = collist[2])


sim_dol_ngoods(k = 1)
