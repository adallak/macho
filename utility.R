###############
### Generate data
###############
#' This function generates  (n xp) data. Borrowed from
#' EqVarDAG package
#'
#' @param n - number of observations
#' @param p - number of features
#' @param pc - sparsity probability
#' @param type - type of the undirected graph
#' @param err - error distribution
#'
#' @return - X (n x p) data
#'
#' @export

gen_DAGdata<-function(n, p, pc, type=c("er","hub","chain"),err= c("g", "nong"), Bmin = 0.5,
                      eq.var = TRUE, lower.sd = 0.7, upper.sd = 1.7, seed = 1234){
      set.seed(seed)
      type = match.arg(type)
      err = match.arg(err)
      if (type=='hub'){
            D<-randomDAG2_hub(p,pc)
      } else if (type=='chain') {
            D<-randomDAG2_chain(p,pc)
      } else {
            D<-randomDAG_mine(p,pc) #randomDAG2_er(p,pc) #
      }
      DAG<-D$DAG
      order<-D$TO
      if (err == 'nong'){
         if (isTRUE(eq.var)){
            stop("currently not supported")
            #errs <- matrix((rbinom(p * n,1,0.5)*2-1)*sqrt(0.8), nrow = p, ncol = n)
         } else {
            sample_prob = rbinom(n * p, 1, 0.4)
            errs1 = rnorm(p * n, mean = 0, sd = lower.sd)
            errs2 = rnorm(p *n, mean = 3, sd = upper.sd)
            errs = matrix(sample_prob * errs1 + (1 - sample_prob) * errs2,
                          nrow = p, ncol = n)
         }
      } else {
            if (isTRUE(eq.var)){
                  errs <- matrix(rnorm(p * n), nrow = p, ncol = n)
            } else {
                  errs = matrix(0, nrow = p, ncol = n)
                  for (i in 1:p){
                        sd = runif(1, lower.sd, upper.sd)
                        errs[i,] = rnorm(n,mean = 0, sd = sd)
                  }
            }
      }
      B<-DAG
      B[B==1]<-runif(sum(DAG),Bmin,1)*(2*rbinom(sum(DAG),1,0.5)-1)
      X <- solve(diag(rep(1, p)) - B) %*% errs##solve(diag(rep(1, p)) - B, errs)
      X <- t(X)
      return(list(DAG=DAG, B=B,
                  X=X, order=order))
}

cv.covEstimation <- function(X, k.list = NULL, n.k = 30,
                             k.min = 0.01, nfolds = 5,
                             method = c("thresholding", "banding", "tapering"),
                             thr.method = c("adaptive", "soft", "hard"),
                             demean = FALSE, scale = FALSE, k.max = NULL) {
   n <- nrow(X)
   p <- ncol(X)
   S = cov(scale(X,center = demean, scale = scale))
   folds = makefolds(n, nfolds)
   method = match.arg(method)
   thr.method = match.arg(thr.method)
   if (is.null(k.list)){
      if (method == "thresholding")
      {
         if (is.null(k.max))
         {
            k.max = k_max(S)
         }
         else{
            k.max = k.max
         }
         k.list = gen.list(n.k, k.max, k.min, S)
      }
      else{
         k.list = seq(1:p)
         n.k = length(k.list)
      }
   }
   errs_fit <- matrix(NA, n.k, nfolds)
   est_cov = array(NA, dim = c(p,p,n.k))
   for (i in 1:nfolds){
      x_tr <- X[-folds[[i]],]
      meanx <- colMeans(x_tr)
      x_tr <- scale(x_tr, center = meanx, scale = scale)
      S_tr <- cov(x_tr)
      iter = 1
      for (k in k.list){
         est_cov[,,iter] = covEstimation(x_tr, k, method = method,
                                         thr.method = thr.method,
                                         demean = demean, scale = scale)$Sigma
         iter = iter + 1
      }
      x_te <- X[folds[[i]], ]
      x_te <- scale(x_te, center = meanx, scale = scale)
      S_te <- cov(x_te)

      for (j in 1:n.k) {
         errs_fit[j, i] <- est.error(est_cov[,,j], S_te)
      }
   }
   ibest_fit <- which.min( rowMeans(errs_fit))

   sigma_fit = covEstimation(X, k.list[ibest_fit], method = method, thr.method = thr.method,
                             demean = demean, scale = scale)$Sigma
   return(list(errs_fit = errs_fit, folds = folds, k.list = k.list,
               cv.k = k.list[ibest_fit], Sigma_fit = sigma_fit, S = S))
}
# ER graph
randomDAG_mine <- function(p,probConnect)
{
      # This function is modified from randomDAG2 function by Jonas Peters
      DAG <- diag(rep(0,p))
      edges = matrix(rbinom(p^2,1,probConnect), p, p)
      edges = lower.tri(edges) * edges
      DAG = edges
      causalOrder = sample(p) #c(1:p)
      P = diag(p)[, causalOrder]
      DAG = P %*% DAG %*% t(P)
      return(list(DAG=DAG,TO=causalOrder))
}

# ER graph
randomDAG2_er <- function(p,probConnect)
{
      # This function is modified from randomDAG2 function by Jonas Peters
      DAG <- diag(rep(0,p))
      causalOrder <- sample(p)
      #cat("causal order is :", causalOrder, "\n")
      for(i in 3:(p))
      {
            node <- causalOrder[i]
            possibleParents <- causalOrder[1:(i-1)]
            Parents <- possibleParents[rbinom(length(possibleParents),1,probConnect)==1]
            DAG[node,Parents] <- rep(1,length(Parents))
      }
      node <- causalOrder[p-1]
      ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
      DAG[causalOrder[2],causalOrder[1]] <- 1
      # causalOrder = c(1,2,3,4,5)
      # DAG = matrix(c(0,0,0,0,0,
      #                1,0,0,0,0,
      #                0,1,0,0,0,
      #                1,0,1,0,0,
      #                1,1,0,1,0), nrow = 5, byrow = TRUE)
      return(list(DAG=DAG,TO=causalOrder))
}
# Chain graph
randomDAG2_chain <- function(p,probConnect)
{
      # This function is modified from randomDAG2 function by Jonas Peters
      DAG <- diag(rep(0,p))
      causalOrder <- sample(p)
      DAG[causalOrder[2],causalOrder[1]] <- 1
      for(i in 3:(p))
      {
            node <- causalOrder[i]
            possibleParents <- causalOrder[1:(i-1)]
            possibleParents <- possibleParents[which(rowSums(DAG[possibleParents,])<4)]
            if (length(possibleParents)>0){
                  Parents <- sample(possibleParents,min(length(possibleParents),2))
                  DAG[node,Parents] <- rep(1,length(Parents))
            }
            DAG[node, causalOrder[i-1]]<-1
      }
      return(list(DAG=DAG,TO=causalOrder))
}
# Hub-and-chain graph
randomDAG2_hub <- function(p,probConnect)
{
      # This function is modified from randomDAG2 function by Jonas Peters
      DAG <- diag(rep(0,p))
      causalOrder <- sample(p)
      Z<-10
      for(i in 1:(p))
      {
            node <- causalOrder[i]
            DAG[node, causalOrder[i-1]]<-1
            if (i>2){DAG[node, causalOrder[sample(seq(min(i-1,Z)),2)]]<-1}
      }
      DAG[causalOrder[2], causalOrder[1]] <- 1
      return(list(DAG=DAG,TO=causalOrder))
}


######################################
#' Compares the estimated and true adjacency matrices
#' @param estAdj - estimated adjacency matrix
#' @param trueAdj - true adjacency matrix
#'
#' @return True Positive, False Positive and True Discovery rates
#'
#' @export
compareGraph <- function (estAdj, trueAdj)
{

   #### This function compares esitmated
   ### adjacency matrix of A with the true Adj matrix
   ### of coefficient matrix A
   ml <- estAdj
   mt <- trueAdj
   p <- dim(ml)[2]
   mt[mt != 0] <- rep(1, sum(mt != 0))
   ml[ml != 0] <- rep(1, sum(ml != 0))
   diffm <- ml - mt
   nmbTrueGaps <- (sum(mt == 0) - p)/2
   fpr <- if (nmbTrueGaps == 0)
      1
   else (sum(diffm > 0)/2)/nmbTrueGaps
   diffm2 <- mt - ml
   nmbTrueEdges <- (sum(mt == 1)/2)
   tpr <- if (nmbTrueEdges == 0)
      0
   else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
   trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
   tdr <- if (sum(ml == 1) == 0) {
      if (trueEstEdges == 0)
         1
      else 0
   }
   else trueEstEdges/(sum(ml == 1)/2)
   return(list(tpr = tpr, fpr = fpr, tdr = tdr))
}

hammingDistance <- function(G1,G2)
   # hammingDistance(G1,G2)
   #
   # Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1!!!!
   #
   # INPUT:  G1, G2     adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
   #
   # OUTPUT: hammingDis Hamming Distance between G1 and G2
   #
   # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
   # All rights reserved.  See the file COPYING for license terms.
{
   allMistakesOne <- FALSE
   if(allMistakesOne)
   {
      Gtmp <- (G1+G2)%%2
      Gtmp <- Gtmp + t(Gtmp)
      nrReversals <- sum(Gtmp == 2)/2
      nrInclDel <- sum(Gtmp == 1)/2
      hammingDis <- nrReversals + nrInclDel
   } else
   {
      hammingDis <- sum(abs(G1 - G2))
      # correction: dist(-,.) = 1, not 2
      hammingDis <- hammingDis - 0.5*sum(G1 * t(G1) * (1-G2) * t(1-G2) +
                                            G2 * t(G2) * (1-G1) * t(1-G1))
   }
   return(hammingDis)
}

linear_experiment <- function(n, p, nsim = 50, type=c("er","hub","chain"),
                              err = c("g", "nong"), eq.var = TRUE,
                              est.method = "chol", threshold = 0.1,
                              method = c("macho", "lingam", "eqvar_BU", "eqvar_TD"),
                              lower.sd = 0.7, upper.sd = 1.2, pc = 0.4,
                              Bmin = 0.5, seed = 1234) {
   type = match.arg(type)
   err = match.arg(err)
   method = match.arg(method)
   store = matrix(0, nrow = nsim, ncol = 4)
   colnames(store) = c("tpr", "fpr", "tdr", "shd")
   true.order.list = matrix(0, nrow = nsim, ncol = p)
   est.order.list  = matrix(0, nrow = nsim, ncol = p)
   for (i in 1:nsim){
      seed = seed + 250*i
      df = gen_DAGdata(n,p, pc = pc, lower.sd = lower.sd,
                       upper.sd = upper.sd, type = type,
                       err = err, eq.var = eq.var, seed = seed)
      X = df$X
      true.dag = df$DAG
      true.order = df$order
      if ( method == "macho"){
         est.order = macho(X, method = "chol", find.first = "TRUE")
         est.dag = DAG_from_Ordering(X, est.order, mtd = est.method, threshold = threshold)
      } else if (method == "eqvar_BU"){
         est_BU = EqVarDAG::EqVarDAG_BU(X, mtd = est.method, threshold = threshold)
         est.order = est_BU$TO
         est.dag = t(est_BU$adj * 1)
      } else if (method == "eqvar_TD"){
         est_TD = EqVarDAG::EqVarDAG_TD(X, mtd = est.method, threshold = threshold)
         est.order = est_TD$TO
         est.dag = t(est_TD$adj * 1)
      } else {
         est_lingam = pcalg::lingam(X)
         est.dag =t(as(est_lingam, "amat") * 1)
         #est.dag[abs(est.dag)> 0] = 1
         g <- igraph::graph_from_adjacency_matrix(est.dag, mode = "directed", diag = FALSE)
         est.order = igraph::topo_sort(g, mode = "out")
      }
      res = compareGraph(est.dag, true.dag)
      store[i, 1] = res$tpr
      store[i, 2] = res$fpr
      store[i, 3] = res$tdr
      store[i, 4] = hammingDistance(est.dag,true.dag)
      true.order.list[i,] = true.order
      est.order.list[i,]  = est.order
   }
   return(list("result" = store, true.order = true.order.list,
               est.order = est.order.list))

}
