Eval <- function(A,S){
  # Objective function for the group balancing problem.
  #
  # arg1: vector of n positive integers.
  # arg2: vector containig either 0, 1 or -1 (0 means no group,
  #       1 for groupA, -1 for groupB).
  #
  # Returns: Positive integer value indicating the difference 
  #          between the weight of both groups.
  
  return(abs(sum(S*A)))
}

Solver <- function(A){
  # Orders and Divides the vector in two group of elements
  # based on the next heuristic; Add the next element of the 
  # list in the smallest group.
  #
  # arg1: vector of n positive integers.
  # 
  # Returns: Positive integer value indicating the difference 
  #          between the weight of both groups.
  
  sum.groupA <- 0
  sum.groupB <- 0
  S <- rep(0, length(A))
  
  for( i in 1:length(A)){
    if( sum.groupA <= sum.groupB){
      sum.groupA <- sum.groupA + A[i]
      S[i] <- 1
    }else{
      sum.groupB <- sum.groupB + A[i]
      S[i] <- (-1)
    }
  }
  
  return (S)
}

TestSolver <- function(P){
  E <- rep(0, length(P))
  for( i in 1:length(P)){
    A <- sort(P[[i]],decreasing = TRUE)
    E[i] <- Eval(A,Solver(A))
  }
  return (E)
}

