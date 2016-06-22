GraphCommunicationProblem <- function (graph,k) {
  size <- length(V(graph))
  edges <- get.edgelist(graph)
  solution.size <- k
  
  evaluate <- function(solution) {
    if (!valid(solution)) {
      stop(paste("The solution is not valid: ", 
                 size))
    }
    itr <- 0
    border <- solution
    oldborder <- c()
    while(length(border) > 0 ){
      adj <- adjacent_vertices(graph,border)
      oldborder <- c(border, oldborder)
      newborder <- c()
      for(i in 1:length(border)) newborder <- c(newborder,setdiff(adj[[i]],oldborder))
      oldborder <- border
      border <- unique(newborder)
      itr <- itr +1
    }
    #Remove last iteration, since evaluted nodes are already communicated
    return (itr-1)
  };
  
  valid <- function(solution) {
    res <- TRUE
    if (length(solution) > size) {
      warning(paste("The solution has to be smaller than the list of nodes: ", 
                 size))
      res <- FALSE
    }
    if (length(solution) < 1) {
      warning(paste("The solucion must have at least one value: ", 
                    size))
      res <- FALSE
    }
    if (anyDuplicated(solution) != 0) {
      warning(paste("The solution can not have duplicated nodes: ", 
                 size))
      res <- FALSE
    }
    if(anyNA(solution)){
      warning(paste("The solution can not contain NA values: ", 
                    size))
      res <- FALSE
    }
    if( (min(solution) < 0) || (max(solution) > size) ){
      warning(paste("The solution can not contain values less than 0 or greater than size: ", 
                    size))
      res <- FALSE
    }
    return (res);
  };
  correct <-  function(solution){
    solution <- unique(solution)
    #solution <- solution[!is.na(solution)]
    solution <- solution[!is.finite(solution)]
    while(length(solution) < solution.size) {
      solution <- c(solution, floor(runif(solution.size-length(solution),1,size)))
      solution <- unique(solution)
    }
    return (solution)
  };
  
  debug <- function(solution){
    itr <- 0
    border <- solution
    oldborder <- c()

    #Initial color of vertexes
    col <- c(rep("red",size))

    
    while(length(border) > 0 ){
      
      cat(border)
      print("\n")
      #for graph visualization
      for(i in 1:length(border)){
        node <- border[[i]] 
        col[node] <- "yellow"
      } 
      for(i in 1:length(col)){
        if(col[i] == "green"){
          col[i] <- "blue"
        }
      }
      for(i in 1:length(oldborder)){
        node <- oldborder[[i]] 
        col[node] <- "green"
      }
      name <- paste("C:/Users/Mikel/Google Drive/Grado en Ing Informática/Cuarto/BH/ImagenesGrafos/graph",itr,".pdf")
      saveSnapshop(graph,name,col)

      #Get neighbours of border
      adj <- adjacent_vertices(graph,border)
      
      #Remove those already in border
      oldborder <- c(border, oldborder)
      newborder <- c()
      for(i in 1:length(border)) newborder <- c(newborder,setdiff(adj[[i]],oldborder))
      oldborder <- border
      border <- unique(newborder)
      #Count of stepts
      itr <- itr +1
    }
    #Remove last iteration, since evaluted nodes are already communicated
    return (itr-1)
  };
  cross <- function (sol.1, sol.2, cross.factor) {
    if(is.null(cross.factor)){
      cross.factor <- 0.5
    }
    if (length(sol.1) != length(sol.2) ){
      stop("The two vectos should be of the same size")
    }
    n <- length(sol.1)
    if (n < 2) {
      stop("The minimum length to use this crossover operator is 2")
    }
    if (class(sol.1[1]) != class(sol.1[2])) {
      stop("The elements in the two vectors should be of the same class")
    }
    else {
      only.in.1 <- setdiff(sol.1,sol.2)
      only.in.2 <- setdiff(sol.2,sol.1)
      if(length(only.in.1) < 2){ return(list(sol.1,sol.2))}
      common.nodes <- intersect(sol.1,sol.2)
      new.sol.1 <- setdiff(sol.1,only.in.1)
      new.sol.2 <- setdiff(sol.2,only.in.2)
      
      availeable.positions <- length(only.in.1)
      if(cross.factor<0 || cross.factor>1){
        warning("Crossing parameter must be between 0 and 1!")
        cross.factor <- 0.5
      }
      pos.to.swap <- floor(availeable.positions * cross.factor)
      if(pos.to.swap < 1){
        pos.to.swap <- 1
      }
      # Swaps pos.to.swap positions of the array randomly chosen from the availeable
      swap.positions <- floor(runif(pos.to.swap,1,availeable.positions))
      for(i in 1:length(swap.positions)){
        new.sol.1 <- c(new.sol.1, only.in.2[swap.positions[i]])
        new.sol.2 <- c(new.sol.2, only.in.1[swap.positions[i]])
      }
      new.sol.1 <- c(new.sol.1, sol.1[c(seq(length(new.sol.1)+1,n))])
      new.sol.2 <- c(new.sol.2, sol.2[c(seq(length(new.sol.2)+1,n))])
    }
    return(list(new.sol.1, new.sol.2))
  };
  mutation <- function (solution, mutation.rate, ...) {
    if (mutation.rate <= 0) {
      stop("The ratio has to be a strictly positive value")
    }
    n <- length(solution)
    id <- sample(n, 1)
    new.value <- sample(size, 1)
    while(new.value %in% solution){ new.value <- sample(size, 1)}
    solution[id] <- new.value
    return(solution)
  };
  environment <- function(old.solution){
    
    new.solution <- old.solution
    old.solutionVal <- evaluate(old.solution)
    
    for(i in 1:length(old.solution)){
      new.solution[i] <- old.solution[i]+1
      if(valid(new.solution)){
        if(evaluate(new.solution) < old.solutionVal)
          return(new.solution)
      }
      new.solution[i] <- old.solution[i]-1
      if(valid(new.solution)){
        if(evaluate(new.solution) < old.solutionVal)
          return(new.solution)
      }
      new.solution <- old.solution
    }
    return (new.solution)
  };
  saveSnapshop <- function(graph,name,col){
    pdf(name)
    layout <-layout.fruchterman.reingold(graph)
    plot(graph, layout=layout, vertex.color=col)
    dev.off()
    return()
  };

  return(list(evaluate = evaluate, valid = valid, debug = debug, 
              saveSnapshop = saveSnapshop, cross = cross, 
              mutation = mutation, correct = correct, environment = environment))
}