#' An S4 class to represent exchange neighborhood for permutations
#'
#' @slot time_total A positive real number representing the total time available. It can be null, indicating that there is no time limit
#' 
setClass(
  Class="AlterNeighborhood", 
  representation=representation(base="numeric", 
                                position.list="matrix", 
                                random="logical",
                                id="numeric", 
                                max = "numeric"))

setValidity(
  Class="AlterNeighborhood", 
  method=function(object) {
   
    n <- length(object@base)
    ## Create all possible alterations
    alter <- matrix(NA,(n*2),2)
    for(i in 1:n){
      alter[i*2-1,1] <- i
      alter[i*2-1,2] <- 1
      alter[i*2,1] <- i
      alter[i*2,2] <- -1
    }
    ## Remobe all invalid alterations
    base <- sort(object@base)
    for(i in 1:(n-1))
      if((base[i]+1) == (base[i+1])){
        alter[i*2-1,2] <- NA
        alter[(i+1)*2,2] <- NA
      }
    if(base[1] == 1)
      alter[2,2] <- NA
    if(base[n] == object@max)
      alter[2*n-1,2] <- NA
    
    alter <- alter[complete.cases(alter),]

    if (!all(object@position.list %in% alter)) {
      stop("Some of the defined alterations are not adequate")
    }
    if (!all(alter %in% object@position.list)) {
      stop("Not all possible alterations are considered")
    }
    if (object@id!=1) {
      stop("The first id to be used has to be 1")
    }
    return (TRUE)
  }
)


# GENERIC METHODS ---------------------------------------------------------

setMethod(
  f="nextNeighbor", 
  signature="AlterNeighborhood", 
  definition=function(neighborhood) {
    if (hasMoreNeighbors(neighborhood)) {
      # We do not use directly the id, but the position in the postion list!!
      nxt <- neighborhood@base
      nxt[neighborhood@position.list[neighborhood@id, 1]] <- 
        nxt[neighborhood@position.list[neighborhood@id, 1]]+
        neighborhood@position.list[neighborhood@id, 2]

      # Update the object
      # obtain the global name of the variable to modify
      objectGlobalName <- deparse(substitute(neighborhood))
      neighborhood@id <- neighborhood@id + 1
      # assign the local variable to the global variable 
      assign(objectGlobalName, neighborhood, envir=parent.frame())  
    } else {
      nxt <- NULL
    }
    return(nxt)
  })


setMethod(
  f="hasMoreNeighbors", 
  signature="AlterNeighborhood", 
  definition=function(neighborhood) {
    return(neighborhood@id <= nrow(neighborhood@position.list))
  })


setMethod(
  f="resetNeighborhood", 
  signature="AlterNeighborhood", 
  definition=function(neighborhood, solution) {
    ## Create new neighborhood with a new solution
    # Order solution vector
    solution <- sort(solution)
    n <- length(solution)
    if(n != length(neighborhood@base)) {
      stop ("The new solution is not of the correct size")
    }
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(neighborhood))
    neighborhood@id <- 1
    neighborhood@base <- solution
    
    alter <- matrix(NA,(n*2),2)
    for(i in 1:n){
      alter[i*2-1,1] <- i
      alter[i*2-1,2] <- 1
      alter[i*2,1] <- i
      alter[i*2,2] <- -1
    }
    ## Remobe all invalid alterations
    
    for(i in 1:(n-1))
      if((solution[i]+1) == (solution[i+1])){
        alter[i*2-1,2] <- NA
        alter[(i+1)*2,2] <- NA
      }
    if(solution[1] == 1)
      alter[2,2] <- NA
    if(solution[n] == neighborhood@max)
      alter[2*n-1,2] <- NA
    
    alter <- alter[complete.cases(alter),]
    
    if (neighborhood@random){
      positions <- alter[sample(nrow(alter)), ]
    }else{
      positions <- alter
    }
    
    neighborhood@position.list <- positions
      
    # If the search is random, shuffle the position list
    if (neighborhood@random) {
      neighborhood@position.list <- neighborhood@position.list[sample(nrow(neighborhood@position.list)), ]
    }
    # assign the local variable to the global variable
    assign(objectGlobalName, neighborhood, envir=parent.frame())  
  })

# CONSTRUCTOR -------------------------------------------------------------

alterNeighborhood <- function(base, random=FALSE, max) {
  n <- length(base)
  ## Order base vector
  base <- sort(base)
  ## Create all possible alterations
  alter <- matrix(NA,(n*2),2)
  for(i in 1:n){
    alter[i*2-1,1] <- i
    alter[i*2-1,2] <- 1
    alter[i*2,1] <- i
    alter[i*2,2] <- -1
  }
  ## Remobe all invalid alterations

  for(i in 1:(n-1))
    if((base[i]+1) == (base[i+1])){
      alter[i*2-1,2] <- NA
      alter[(i+1)*2,2] <- NA
    }
  if(base[1] == 1)
    alter[2,2] <- NA
  if(base[n] == max)
    alter[2*n-1,2] <- NA
  
  alter <- alter[complete.cases(alter),]
  
  if (random){
    positions <- alter[sample(nrow(alter)), ]
  }else{
    positions <- alter
  }

  
  obj <- new("AlterNeighborhood", base=base, position.list=positions, 
             random=random, id=1,max=max)
  return(obj)
}

