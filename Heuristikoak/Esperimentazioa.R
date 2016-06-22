# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Graph Communication Problem Experimentation       #
# ................................................. #
# Author: Mikel Dalmau (mdalmau002@ikasle.ehu.eus)  #
# Date:   April, 2016                               #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# REQUIRED PACKAGES AND INFORMATION --------------------------------------------
library(metaheuR)
library(igraph)
library(scmamp)
library(ggplot2)
library(multtest)
source("C:/Users/Mikel/Google Drive/Grado en Ing Informática/Cuarto/BH/BH/Praktika/GraphCommunicationProblem.R")
source("C:/Users/Mikel/Google Drive/Grado en Ing Informática/Cuarto/BH/BH/Praktika/basicMemeticAlgorithm.R")
source("C:/Users/Mikel/Google Drive/Grado en Ing Informática/Cuarto/BH/BH/Praktika/AlterNeighborhood.R")

# USER-DEFINED FUNCTIONS --------------------------------------------------


generateGraphAging <- function(size, max.arcs=1) {
  # Simple function to generate random evolving graphs
  # ARGS:
  #  size - Number of nodes in the graph
  #  max.arcs - Maximum number of arcs that the nodes will have
  # RETURNS:
  #  Object of type "igraph" representing the generated random evolving graph
  #
  rnd.graph <- aging.ba.game(n=size, pa.exp = 0.3, aging.exp=-0.4, m=max.arcs, 
                             directed=FALSE)
  # The above function can generate multiple arcs between two nodes
  # To remove this, we use the adjacency matrix
  adj.matrix <- get.adjacency(rnd.graph)!=0
  # Finally, the graph is re-generated with the information
  rnd.graph <- graph.adjacency(adj.matrix, mode="undirected")
  return(rnd.graph)
}

generateGraphErdos <- function(size, p) {
  # Simple function to generate random graphs
  # ARGS:
  #  size - Number of nodes in the graph
  #  p - Probablity to create an edge between two nodes
  # RETURNS:
  #  Object of type "igraph" representing the generated random graph
  #
  rnd.graph <-erdos.renyi.game(size, p)
  # The above function can generate multiple connected components
  # To remove this, we decompose the graph and take the top of the list
  rnd.graphs <- decompose.graph(rnd.graph, mode = c("weak", "strong"), 
                                max.comps = NA, min.vertices = 0)
  
  max.graph <- rnd.graphs[[1]]
  max.length <- length(V(max.graph))
  for(i in 1:length(rnd.graphs)){
    l <- length(V(rnd.graphs[[i]]))
    if( l > max.length ){
      max.graph <- rnd.graphs[[i]]
      max.length <- l
    }
  }  
  return(max.graph)
}

generateGraphWatts <- function(size, dim, nei, p) {
  # Simple function to generate random graphs
  # ARGS:
  #  dim	- Integer constant, the dimension of the starting lattice.
  #  size - Integer constant, the size of the lattice along each dimension.
  #  nei	- Integer constant, the neighborhood within which the vertices of 
  #         the lattice will be connected.
  #  p	- Real constant between zero and one, the rewiring probability.
  # RETURNS:
  #  Object of type "igraph" representing the generated random graph
  #
  rnd.graph <-watts.strogatz.game(dim,size,nei, p)
  # The above function can generate multiple edges between nodes and loops
  # To remove this, we use function simplify
  rnd.graphs <- simplify(rnd.graph, remove.multiple = TRUE, remove.loops = TRUE)
  # The above function can generate multiple connected components
  # To remove this, we decompose the graph and take the top of the list
  rnd.graphs <- decompose.graph(rnd.graph, mode = c("weak", "strong"), 
                                max.comps = NA, min.vertices = 0)
  
  max.graph <-rnd.graphs[[1]]
  max.length <- length(V(max.graph))
  for(i in 1:length(rnd.graphs)){
    l <- length(V(rnd.graphs[[i]]))
    if( l > max.length ){
      max.graph <- rnd.graphs[[i]]
      max.length <- l
    }
  }  
  return(max.graph)
}

runMemetic <- function (graph, max.time,k,populationSize,mrate) {
  # Function to run a simple memetic algorithm on a graph
  # ARGS:
  #  graph - Instance of the problem to solve
  #  max.time - Maximum time given to the algorithm (in seconds)
  # RETURNS:
  #  The evaluation function of the best found solution
  #
  n <- length(V(graph))
  
  # Initialization of neighbourhood object for local search
  argsn <- list()
  argsn$base <- sample(1:n,k,replace=F)
  argsn$random <- TRUE
  argsn$max <- n
  alt.ngh <- do.call(alterNeighborhood,argsn)
  
  
  createRndSolution <- function(x) {
    sol <- sample(1:n, k, replace=F)
    return(sol)
  }
  population <- lapply(1:populationSize, FUN=createRndSolution)
  
  gcp <- GraphCommunicationProblem(graph,k)
  args <- list()
  args$evaluate             <- gcp$evaluate
  args$initial.population   <- population
  args$selectSubpopulation  <- elitistSelection
  args$selection.ratio      <- 0.5
  args$selectCross          <- tournamentSelection
  args$cross                <- gcp$cross
  args$cross.factor         <- 0.5
  args$mutate               <- gcp$mutation
  args$mutation.rate        <- mrate
  args$neighbor.selector    <- firstImprovementSelector
  args$neighborhood         <- alt.ngh
  args$verbose              <- TRUE
  args$non.valid            <- "correct"
  args$valid                <- gcp$valid
  args$correct              <- gcp$correct
  args$resources            <- cResource(evaluations= 100000 ,time = max.time)
  
  res <- do.call(basicMemeticAlgorithm, args)
  
  plotProgress(res)
  name <- paste("memetic",populationSize, k,mrate,".png",sep="_")
  ggsave(name, width = 10, height = 10)
  
  return(getEvaluation(res))
}

runLocalSearch <- function (graph, max.time,k) {
  # Function to run a simple memetic algorithm on a graph
  # ARGS:
  #  graph - Instance of the problem to solve
  #  max.time - Maximum time given to the algorithm (in seconds)
  # RETURNS:
  #  The evaluation function of the best found solution
  #
  n <- length(V(graph))
  
  resources <- cResource(evaluations= 100000,time = max.time)
  # Initialization of neighbourhood object for local search
  argsn <- list()
  argsn$base <- sample(1:n,k,replace=F)
  argsn$random <- TRUE
  argsn$max <- n
  alt.ngh <- do.call(alterNeighborhood,argsn)
  
  best <- n
  
  while(!isFinished(resources)) {
    
  t0 <- Sys.time()
  
  gcp <- GraphCommunicationProblem(graph,k)
  args <- list()
  args$neighborhood         <- alt.ngh
  args$evaluate             <- gcp$evaluate
  args$selector              <- firstImprovementSelector
  args$non.valid            <- "correct"
  args$valid                <- gcp$valid
  args$correct              <- gcp$correct
  args$resources            <- cResource(evaluations= 100000,time = max.time)
  args$initial.solution     <- sample(1:n, k, replace=F)

  res <- do.call(basicLocalSearch, args)
  
  eval <- length(numeric(getEvaluation(res)))
  if(eval < best){
    best <- eval  
  }
    addConsumed(resources, t=as.numeric(Sys.time() - t0), 
              ev=1, it=1)
  }
  return(toString(best))
}

runRandom <- function (graph, max.time, k) {
 
  n <- length(V(graph))
  gcp <- GraphCommunicationProblem(graph,k)
  resources <- cResource(evaluations= 100000,time = max.time)
  best <- n
  
  while(!isFinished(resources)){
    t0 <- Sys.time()
    solution <- sample(1:n,k,replace=F)
    res <- gcp$evaluate(solution)
    if(res < best){
      best <- res  
    }
    addConsumed(resources, t=as.numeric(Sys.time() - t0), 
                ev=1)
  }
  return(toString(best))
}

# EXPERIMENTATION ---------------------------------------------------------

# Generation of graphs (small, medium and big)
generateSMALL <- FALSE
generateMEDIUM <- TRUE
generateBIG <- FALSE

# Execution of Algorithms on graphs (small, medium and big)
SMALL <- FALSE
MEDIUM <- FALSE
BIG <- FALSE

DOHISTO <- FALSE

# Number of each graphs to create
num.graphs <- 2
seed <- floor(runif(n=1,max=1000,min=1))

if(generateSMALL){
  set.seed(seed)
  aging.small.graphs <- lapply(rep(50, num.graphs), FUN=generateGraphAging, max.arcs=3)
  erdos.small.graphs <- lapply(rep(50, num.graphs), FUN=generateGraphErdos, p= 2/50)
  watts.small.graphs <- lapply(rep(50, num.graphs), FUN=generateGraphWatts, dim=1, nei=1, p=0.2)
 
}

if(generateMEDIUM){
  set.seed(seed)
  aging.medium.graphs <- lapply(rep(100, num.graphs), FUN=generateGraphAging, max.arcs=4)
  erdos.medium.graphs <- lapply(rep(100, num.graphs), FUN=generateGraphErdos, p= 2/100)
  watts.medium.graphs <- lapply(rep(100, num.graphs), FUN=generateGraphWatts, dim=1, nei=1, p=0.2)
}  

if(generateBIG){
  set.seed(seed)
  aging.big.graphs   <- lapply(rep(200, num.graphs), FUN=generateGraphAging, max.arcs=7)
  erdos.big.graphs   <- lapply(rep(200, num.graphs), FUN=generateGraphErdos, p= 3/200)
  watts.big.graphs   <- lapply(rep(200, num.graphs), FUN=generateGraphWatts, dim=1, nei=1, p=0.2)   
}  
  
# Application of the algorithms ---------------------------------------------------------

## SMALL ·······················································································································

if(SMALL){
  time        <- 60
  repetitions <- 5
  populationSize <- 20
  k <- 5
  
  memetic.aging.small <- lapply (1:length(aging.small.graphs),
                          FUN=function(i) {
                            graph <- aging.small.graphs[[i]]
                            res <- sapply(1:repetitions,
                                                  FUN=function(j) {
                                                    message("Running repetition ",j," of memetic algorithm in aging small graph#",i)
                                                    runMemetic(graph, max.time=time,k,populationSize)
                                                  })
                            return (cbind(Graph=paste("aging_small",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                          })
  
  local.aging.small <- lapply (1:length(aging.small.graphs),
                                 FUN=function(i) {
                                   graph <- aging.small.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                         FUN=function(j) {
                                                           message("Running repetition ",j," of Local Search algorithm in aging small graph#",i)
                                                           runLocalSearch(graph, max.time=time,k)
                                                         })
                                   return (cbind(Graph=paste("aging_small",i,sep="_"), Algorithm="Local", Evaluation=res))
                                 })
  
  random.aging.small <- lapply (1:length(aging.small.graphs),
                               FUN=function(i) {
                                 graph <- aging.small.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Random Search algorithm in aging small graph#",i)
                                                 runRandom(graph, max.time=time, k)
                                               })
                                 return (cbind(Graph=paste("aging_small",i,sep="_"), Algorithm="Random", Evaluation=res))
                               })
  
  
  memetic.erdos.small <- lapply (1:length(erdos.small.graphs),
                                 FUN=function(i) {
                                   graph <- erdos.small.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in Erdos small  graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("Erdos_small",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.erdos.small <- lapply (1:length(erdos.small.graphs),
                                FUN=function(i) {
                                  graph <- erdos.small.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Local Search algorithm in Erdos small graph#",i)
                                                  runLocalSearch(graph, max.time=time,k)
                                                })
                                  return (cbind(Graph=paste("Erdos_small",i,sep="_"), Algorithm="Local", Evaluation=res))
                                })
  
  random.erdos.small <- lapply (1:length(erdos.small.graphs),
                                FUN=function(i) {
                                  graph <- erdos.small.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in  Erdos small graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("Erdos_small",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
  
  memetic.watts.small <- lapply (1:length(watts.small.graphs),
                                 FUN=function(i) {
                                   graph <- watts.small.graphs[[1]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in watts small graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize,mrate[j])
                                                 })
                                   return (cbind(Graph=paste("watts_small",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.watts.small <- lapply (1:length(watts.small.graphs),
                               FUN=function(i) {
                                 graph <- watts.small.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in watts small graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("watts_small",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.watts.small <- lapply (1:length(watts.small.graphs),
                                FUN=function(i) {
                                  graph <- watts.small.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in watts small graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("watts_small",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
}
## MEDIUM ·······················································································································

if(MEDIUM){
  time        <- 180
  repetitions <- 5
  populationSize <- 30
  k <- 7
  
  memetic.aging.medium <- lapply (1:length(aging.medium.graphs),
                                 FUN=function(i) {
                                   graph <- aging.medium.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in aging medium graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("aging_medium",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.aging.medium <- lapply (1:length(aging.medium.graphs),
                               FUN=function(i) {
                                 graph <- aging.medium.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in aging medium graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("aging_medium",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.aging.medium <- lapply (1:length(aging.medium.graphs),
                                FUN=function(i) {
                                  graph <- aging.medium.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in aging medium graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("aging_medium",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
  
  
  memetic.erdos.medium <- lapply (1:length(erdos.medium.graphs),
                                 FUN=function(i) {
                                   graph <- erdos.medium.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in Erdos medium graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("Erdos_medium",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.erdos.medium <- lapply (1:length(erdos.medium.graphs),
                               FUN=function(i) {
                                 graph <- erdos.medium.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in Erdos medium graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("Erdos_medium",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.erdos.medium <- lapply (1:length(erdos.medium.graphs),
                                FUN=function(i) {
                                  graph <- erdos.medium.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in  Erdos medium graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("Erdos_medium",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
  
  memetic.watts.medium <- lapply (1:length(watts.medium.graphs),
                                 FUN=function(i) {
                                   graph <- watts.medium.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in watts medium graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("watts_medium",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.watts.medium <- lapply (1:length(watts.medium.graphs),
                               FUN=function(i) {
                                 graph <- watts.medium.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in watts medium graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("watts_medium",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.watts.medium <- lapply (1:length(watts.medium.graphs),
                                FUN=function(i) {
                                  graph <- watts.medium.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in watts medium graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("watts_medium",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
}
## BIG ·······················································································································

if(BIG){
  time        <- 200
  repetitions <- 5
  populationSize <- 30
  k <- 9
  
  memetic.aging.big <- lapply (1:length(aging.big.graphs),
                                 FUN=function(i) {
                                   graph <- aging.big.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in aging big graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("aging_big",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.aging.big <- lapply (1:length(aging.big.graphs),
                               FUN=function(i) {
                                 graph <- aging.big.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in aging big graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("aging_big",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.aging.big <- lapply (1:length(aging.big.graphs),
                                FUN=function(i) {
                                  graph <- aging.small.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in aging big graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("aging_big",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
  
  
  memetic.erdos.big <- lapply (1:length(erdos.big.graphs),
                                 FUN=function(i) {
                                   graph <- erdos.big.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in Erdos big graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("Erdos_big",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.erdos.big <- lapply (1:length(erdos.big.graphs),
                               FUN=function(i) {
                                 graph <- erdos.big.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in Erdos big graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("Erdos_big",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.erdos.big <- lapply (1:length(erdos.big.graphs),
                                FUN=function(i) {
                                  graph <- erdos.big.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in  Erdos big graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("Erdos_big",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
  
  memetic.watts.big <- lapply (1:length(watts.big.graphs),
                                 FUN=function(i) {
                                   graph <- watts.big.graphs[[i]]
                                   res <- sapply(1:repetitions,
                                                 FUN=function(j) {
                                                   message("Running repetition ",j," of memetic algorithm in watts big graph#",i)
                                                   runMemetic(graph, max.time=time,k,populationSize)
                                                 })
                                   return (cbind(Graph=paste("watts_big",i,sep="_"), Algorithm="Memetic", Evaluation=res))
                                 })
  
  local.watts.big <- lapply (1:length(watts.big.graphs),
                               FUN=function(i) {
                                 graph <- watts.big.graphs[[i]]
                                 res <- sapply(1:repetitions,
                                               FUN=function(j) {
                                                 message("Running repetition ",j," of Local Search algorithm in watts big graph#",i)
                                                 runLocalSearch(graph, max.time=time,k)
                                               })
                                 return (cbind(Graph=paste("watts_big",i,sep="_"), Algorithm="Local", Evaluation=res))
                               })
  
  random.watts.big <- lapply (1:length(watts.big.graphs),
                                FUN=function(i) {
                                  graph <- watts.big.graphs[[i]]
                                  res <- sapply(1:repetitions,
                                                FUN=function(j) {
                                                  message("Running repetition ",j," of Random Search algorithm in watts big graph#",i)
                                                  runRandom(graph, max.time=time, k)
                                                })
                                  return (cbind(Graph=paste("watts_big",i,sep="_"), Algorithm="Random", Evaluation=res))
                                })
}

# RESULTS ANALYSIS --------------------------------------------------------

## GROUPING FUNCTION FOR RESULTS  ·····························································
groupTable <- function(table){
  aux <- rbind(table[[1]],table[[2]])
  #for(i in 3:num.graphs ){
  # aux <- rbind(aux,table[[i]])
  #}
  return(aux)
}
## GROUPING FUNCTION FOR RESULTS ································································

if(DOHISTO){
group1 <- groupTable(memetic.aging.small)
group2 <- groupTable(local.aging.small)
group3 <- groupTable(random.aging.small)
group4 <- groupTable(memetic.erdos.small)
group5  <- groupTable(local.erdos.small)
group6 <- groupTable(random.erdos.small)
group7  <-groupTable(memetic.watts.small)
group8  <- groupTable(local.watts.small)
group9 <-  groupTable(random.watts.small)
AllTable <- rbind(group1,group2,group3,group4,group5,group6,group7,group8,group9)

Allframe <- data.frame(AllTable)



g <- ggplot(resss, aes(x=Graph, y=Evaluation, fill=Algorithm)) + 
      geom_bar(stat="summary", fun.y=mean, fun.ymin=function(x){mean(x)-sd(x)}, fun.ymax=function(x){mean(x)+sd(x)},  position="dodge") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", position="dodge")

ggsave("HistoSmall.png", width = 20, height = 10)

}