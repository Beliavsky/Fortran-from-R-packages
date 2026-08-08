#' Global optimization of a blackbox function given a finite budget of noisy evaluations,
#' via the Parallel Optimistic Optimization algorithm.
#' The knowledge of the function's smoothness is not required.
#' @title Parallel Optimistic Optimization
#' @param f function to maximize.
#' @param horizon maximum number of function evaluations.
#' @param nu scalar (> 0) assessing the complexity of the function, along with \code{rho} (see the near optimality definition in the reference below).
#' @param rhomax number of equidistant \code{rho} values in [0,1], that are used by the corresponding HOO subroutines, see Details.
#' @param noise.level scalar bound on the noise value.
#' @details Only 1-dimensional functions defined on [0, 1] are handled so far.
#' POO uses Hierarchical Optimistic Optimisation (HOO) as a subroutine, whose number is set by \code{rhomax}.
#' \code{POO} handles more difficult functions than \code{\link[OOR]{StoSOO}}.
#' @return Random point evaluated by the best HOO, in the form of a list with elements:
#' \itemize{
#'  \item par parameter value at this point,
#'  \item value noisy value at \code{par},
#'  \item best_rho best \code{rho} value.
#' }
#' @export
#' @importFrom methods new
#' @importFrom stats runif
#' @author
#' M. Binois (translation in R code), J.-B. Grill, M. Valko and R. Munos (Python code)
#' @references
#' J.-B. Grill, M. Valko and R. Munos (2015), Black-box optimization of noisy functions with unknown smoothness,
#'   \emph{NIPS}, 667-675 \url{https://inria.hal.science/hal-01222915}. Python code: \url{https://team.inria.fr/sequel/software/POO/}.
#' @examples
#' \dontrun{
#' #------------------------------------------------------------
#' # Maximization with POO
#' #------------------------------------------------------------
#' set.seed(10)
#' noise.level <- 0.05
## ' XX <- c(NA)
#'
#' ## Define and display objective
## ' ftest <- function(x){XX <<- c(XX, x); return(double_sine(x) + runif(1, min = -noise.level, max = noise.level))}
#' ftest <- function(x){return(double_sine(x) + runif(1, min = -noise.level, max = noise.level))}
#' xgrid <- seq(0, 1, length.out = 1000)
#' plot(xgrid, sapply(xgrid, double_sine), type = 'l', ylab = "double_sine(x)", xlab = 'x')
#'
#' ## Optimization
#' Sol <- POO(ftest, horizon = 1000, noise.level = noise.level)
#'
#' ## Display result
#' abline(v = Sol$par)
#'}
POO <- function(f, horizon = 100, noise.level, rhomax = 20, nu = 1){
  alpha <- log(horizon)*noise.level^2
  rhos <- seq(0, 1, length.out = rhomax)
  ypoo <- rep(NA, horizon)

  tree <- Tree(c(0,1), 0, 0, rhos, f)
  compt <- 0
  empp <- rep(0, rhomax)
  nsam <- rep(0, rhomax)

  while(compt <= horizon){
    for(k in 1:rhomax){
      smp <- Tree_sample(tree, alpha, nu, k)
      empp[k] <- empp[k] + smp$noisyvalue
      compt <- compt + smp$existed
      nsam[k] <- nsam[k] + 1
      if(smp$existed && compt <= horizon){
        best_k <- which.max(empp[which(nsam > 0)] / nsam[which(nsam > 0)])
      }
    }
  }
  result <- Tree_sample(tree, alpha, nu, best_k)
  return(list(par = result$evaluated, value = result$noisyvalue, best_rho = rhos[best_k]))
}


## ' Tree structure for POO
## ' @title Tree constructor
## ' @param support support of the leaf
## ' @param father scalar id of the father node
## ' @param depth scalar number of edge to the root of the tree
## ' @param rhos
## ' @param fun function to optimize
Tree <- function(support, father, depth, rhos, fun){
  tmp <- new.env()
  tmp$tree <- list(Leaf(support, father, depth, rhos))
  tmp$fun <- fun
  return(tmp)
}


# ## ----------------
# ## TREE CLASS definition
# ## ----------------
#
# #' S3 Tree Class
# #' @export
# setClass("Tree", representation(tree = "list", fun = "function"))

std_split <- function(support){
  m <- sum(support)/2
  return(list(c(support[1], m), c(m, support[2])))
}

std_rpoint <- function(support){
  return(runif(n = 1, min = support[1], max = support[2]))
}

## ----------------
## Add_children method
## ----------------

add_children <- function(Tree, id){
  child_supports <- std_split(Tree$tree[[id]]@support)
  for(i in 1:length(child_supports)){
    newChild <- Leaf(child_supports[[i]], id, Tree$tree[[id]]@depth + 1, Tree$tree[[id]]@rhos, Tree$tree[[id]]@bbox)
    Tree$tree <- c(Tree$tree, list(newChild))
    Tree$tree[[id]]@children <- c(Tree$tree[[id]]@children, length(Tree$tree))
  }
}

# if(!isGeneric("add_children")) {
#   setGeneric(name = "add_children",
#              def = function(object, ...) standardGeneric("add_children")
#   )
# }
#
# setMethod("add_children", "Tree",
#           function(object, id) {
#             Tree.add_children(object = object, id = id)
#           }
# )

## ----------------
## Explore method
## ----------------

explore <- function(Tree, id, k){
  if(!Tree$tree[[id]]@visited[k]) return(id)
  if(length(Tree$tree[[id]]@children) == 0){
    add_children(Tree, id)
    return(sample(Tree$tree[[id]]@children, size = 1))
  }
  idchild <- which.max(lapply(Tree$tree[Tree$tree[[id]]@children], function(leaf) return(leaf@bvalue[k])))
  idchild <- Tree$tree[[id]]@children[idchild]
  return(explore(Tree, idchild, k))
}

# if(!isGeneric("explore")) {
#   setGeneric(name = "explore",
#              def = function(object, ...) standardGeneric("explore")
#   )
# }
#
# setMethod("explore", "Tree",
#           function(object, id, k) {
#             Tree.explore(object = object, id = id, k = k)
#           }
# )

## ----------------
## Update_node method
## ----------------

update_node <- function(Tree, id, alpha, nu, k){
  mean <- Tree$tree[[id]]@rewards[k] / Tree$tree[[id]]@visited[k]
  ucb <- sqrt(2*alpha/Tree$tree[[id]]@visited[k])
  metric <- nu * Tree$tree[[id]]@rhos[k]^Tree$tree[[id]]@depth^k
  Tree$tree[[id]]@uvalue[k] <- mean + ucb + metric
  Tree$tree[[id]]@muvalue[k] <- mean - ucb - metric
}


# if(!isGeneric("update_node")) {
#   setGeneric(name = "update_node",
#              def = function(object, ...) standardGeneric("update_node")
#   )
# }
#
# setMethod("update_node", "Tree",
#           function(object, id, alpha, nu, k) {
#             Tree.update_node(object = object, id = id, alpha = alpha, nu = nu, k = k)
#           }
# )


## ----------------
## Update_path method
## ----------------

update_path <- function(Tree, id, reward, alpha, nu, k){
  Tree$tree[[id]]@rewards[k] <- Tree$tree[[id]]@rewards[k] + reward
  Tree$tree[[id]]@visited[k] <- Tree$tree[[id]]@visited[k] + 1
  update_node(Tree, id = id, alpha = alpha, nu = nu, k = k)
  if(length(Tree$tree[[id]]@children) == 0){
    Tree$tree[[id]]@bvalue[k] <- Tree$tree[[id]]@uvalue[k]
    Tree$tree[[id]]@mbvalue[k] <- Tree$tree[[id]]@muvalue[k]
  }else{
    Tree$tree[[id]]@bvalue[k] <- min(Tree$tree[[id]]@uvalue[k],
                                       max(unlist(lapply(Tree$tree[Tree$tree[[id]]@children],
                                                  function(leaf) return(leaf@bvalue[k])))))
    Tree$tree[[id]]@mbvalue[k] <- max(Tree$tree[[id]]@muvalue[k],
                                       max(unlist(lapply(Tree$tree[Tree$tree[[id]]@children],
                                                  function(leaf) return(leaf@mbvalue[k])))))
  }
  if(Tree$tree[[id]]@father != 0) update_path(Tree, Tree$tree[[id]]@father, reward, alpha, nu, k)
}

# if(!isGeneric("update_path")) {
#   setGeneric(name = "update_path",
#              def = function(object, ...) standardGeneric("update_path")
#   )
# }
#
# setMethod("update_path", "Tree",
#           function(object, id, reward, alpha, nu, k) {
#             Tree.update_path(object = object, id = id, reward = reward, alpha = alpha, nu = nu, k = k)
#           }
# )


## ----------------
## Sample method
## ----------------


Tree_sample <- function(Tree, alpha, nu, k){
  idleaf <- explore(Tree = Tree, id = 1, k = k)
  existed <- FALSE
  if(is.na(Tree$tree[[idleaf]]@noisyvalue)){
    x <- std_rpoint(Tree$tree[[idleaf]]@support)
    Tree$tree[[idleaf]]@evaluated <- x
    Tree$tree[[idleaf]]@noisyvalue <- Tree$fun(x)
    existed <- TRUE
  }
  update_path(Tree = Tree, id = idleaf, reward = Tree$tree[[idleaf]]@noisyvalue, alpha = alpha, nu = nu, k = k)
  return(list(evaluated = Tree$tree[[idleaf]]@evaluated, noisyvalue = Tree$tree[[idleaf]]@noisyvalue, existed = existed))
}

# if(!isGeneric("sample")) {
#   setGeneric(name = "sample",
#              def = function(object, ...) standardGeneric("sample")
#   )
# }
#
# setMethod("sample", "Tree",
#           function(object, alpha, nu, k) {
#             Tree.sample(object = object, alpha = alpha, nu = nu, k = k)
#           }
# )




## ----------------
## LEAF CLASS definition
## ----------------


## ' S3 Leaf Class
## ' @importFrom methods new
## ' @importFrom stats runif
setClass("Leaf",
         representation(
           noisyvalue = "numeric",
           evaluated = "numeric",
           # meanvalue = "numeric",
           visited = "numeric",
           rewards = "numeric",
           bvalue = "numeric",
           uvalue = "numeric",
           mbvalue = "numeric",
           muvalue = "numeric",
           # bestempiricalnoised = "numeric",
           # bestempiricalmean = "numeric",
           support = "numeric",
           father = "integer",
           depth = "integer",
           rhos = "numeric",
           # bbox = "list",
           children = "integer"
         )#,
         # prototype(
         #   noisyvalue = NA_real_,
         #   evaluated = FALSE,
         #   meanvalue = NA_real_,
         #   bestempiricalnoised = -Inf,
         #   bestempiricalmean = -Inf
         # )
)




`Leaf` <- function(support, father, depth, rhos, bbox){
  newLeaf <- new("Leaf")
  newLeaf@noisyvalue <- NA_real_
  newLeaf@evaluated <- NA_real_
  # newLeaf@meanvalue <- NA_real_
  newLeaf@visited <- rep(0, length(rhos))
  newLeaf@rewards <- rep(0, length(rhos))
  newLeaf@bvalue <- rep(Inf, length(rhos))
  newLeaf@uvalue <- rep(Inf, length(rhos))
  newLeaf@mbvalue <- rep(-Inf, length(rhos))
  newLeaf@muvalue <- rep(-Inf, length(rhos))
  # newLeaf@bestempiricalnoised = -Inf
  # newLeaf@bestempiricalmean = -Inf
  newLeaf@rhos <- rhos
  newLeaf@depth <- as.integer(depth)
  newLeaf@support <- support
  newLeaf@father <- as.integer(father)

  return(newLeaf)
}


# #' S3 Box Class
# #' @export
# setClass("Box",
#          representation(
#            support = "numeric",
#            fmax = "numeric",
#            center = "numeric",
#            rpoint = "numeric",
#            split = "matrix"#,
#            # f_noised = "function",
#            # f_mean = "function"
#          )
# )
#
# `Box` <- function(fmax){
#   newBox <- new("Box")
#   newBox <- fmax
#   return(newBox)
# }
#
# # std_part.Box <- function(){
# #
# # }
# #
# # # s: one variable per row, 2 columns for min and max
# # std_center <- function(s){
# #   m <- rowMeans(s)
# #   return(cbind(s[,1]))
# # }
#
# std_split <- function(s){
#
# }



