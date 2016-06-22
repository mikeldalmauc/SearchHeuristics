library("igraph")
library(metaheuR)
source("C:/Users/TOSHIBA C855 1WC/Desktop/UPV  EHU/15-16/16/BH/practica/GraphCommunicationProblem.R")

# The erdos.renyi.game function, when you set the probability to 1, generates complete graphs.
# g1 <<-erdos.renyi.game(50, 3/50)


# The degree.sequence.game function generates random graphs with a prescribed degree distribution.
# g2 <- degree.sequence.game( c(3,3,3,2,1,1,1), method="vl" )

# The watts.strogatz.game function generates small-world networks.
# dim
# Integer constant, the dimension of the starting lattice.
# size
# Integer constant, the size of the lattice along each dimension.
# nei
# Integer constant, the neighborhood within which the vertices of the lattice will be connected.
# p
# Real constant between zero and one, the rewiring probability.
# loops
# Logical scalar, whether loops edges are allowed in the generated graph.
# multiple
# Logical scalar, whether multiple edges are allowed int the generated graph
# Note that this function might create graphs with loops and/or multiple edges. You can use simplify to get rid of these.
#g1 <- watts.strogatz.game(1, 100, 5, 0.05)
#g1 <-simplify(g1,remove.multiple = TRUE,remove.loops = TRUE,edge.attr.comb = igraph_opt("edge.attr.comb") )
#write.graph(g1, file='C:/Users/Mikel/Google Drive/Grado en Ing Informática/Cuarto/BH/GrafosDePruebas/my_graph1.txt', format="edgelist")
#g4 <- read.graph("C:/Users/Mikel/Google Drive/Grado en Ing Informática/Cuarto/BH/GrafosDePruebas/my_graph1.txt", format="edgelist",directed = FALSE)
#as.undirected(g4, mode = c("collapse", "each", "mutual"), edge.attr.comb = igraph_opt("edge.attr.comb"))

localSearch <- function (evaluate, 
                         new.solution, 
                         #environment,
                         valid, ...) {
  
  new.solutionVal <- evaluate(new.solution)
  old.solutionVal <- 0
  
  while(new.solutionVal != old.solutionVal){
    new.solution <- environment(new.solution, evaluate, valid)
    new.solutionVal <- evaluate(new.solution)
    old.solutionVal <- new.solutionVal
  }
  return(new.solution)
}

environment <- function(old.solution, 
                       evaluate, 
                       valid, ...){
  
  new.solution <- old.solution
  old.solutionVal <- evaluate(old.solution)
  
  for(i in 1:length(old.solution)){
    new.solution[i] <- old.solution[i]+1
    if(valid(new.solution)){
      new.solutionVal <- evaluate(new.solution)
      if(new.solutionVal < old.solutionVal)
        return(new.solution)
    }
    new.solution[i] <- old.solution[i]-1
    if(valid(new.solution)){
      new.solutionVal <- evaluate(new.solution)
      if(new.solutionVal < old.solutionVal)
        return(new.solution)
    }
  }
  return (old.solution)
}
  
  