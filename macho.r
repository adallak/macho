###############
### Find a topological ordering
###############
#' Infer  DAG using topological ordering. Adopted from EqVarDAG package.
#' @param X: data in n x p matrix
#' @param method: ordering estimation methods; default is chol
#' @param first.entry: possible index for the first entry (optional).
#' @param find.first: whether to use specific method to find an index for the first entry(optional).
#' @return topological ordering
#'
macho <- function(X, method = c("chol","det"), first.entry = NULL, find.first = FALSE) {
      p = dim(X)[2]
      n = dim(X)[1]
      S = var(X)
      current.chol = NULL
      method = match.arg(method)
      search = 1:p
      if (!is.null(first.entry) && isTRUE(find.first)) {
         stop("only one of first.entry and find.first should be specified")
      }
      if (!is.null(first.entry)) {
         if(first.entry> p){
            stop("first.entry should be less than the number of dimensions")
         }
         if(first.entry<0){
            stop("first.entry should be positive integer")
         }
      }
      if (is.null(first.entry)){
         order = which.min(diag(S))
      } else {
         order = first.entry
      }
      search = search[-order]
      while(length(search)>0){
            if(isFALSE(find.first)){
                  if (method == "det"){
                     find.min = find.min.det(S, order, search)
                  } else {
                         find.min = find.min.chol(S, order, search, current.chol)
                         current.chol = find.min$current.chol
                  }
            } else {
               if (method == "det"){
                  find.min = find.min.det.ext(S, order, search)
               } else {
                  find.min = find.min.chol.ext(S, order, search, current.chol)
                  current.chol = find.min$current.chol
               }
            }
            order = find.min$current
            search = find.min$search
      }
      return(order)
}


macho_HD <- function(X, method = c("chol","det"), first.entry = NULL, find.first = FALSE,
                     lambda = NULL, L = "NULL", crit = c("bic", "ebic"), start = NULL) {
   p = dim(X)[2]
   n = dim(X)[1]
   crit = match.arg(crit)
   S = covglasso::covglasso(X, L = L, lambda = lambda, crit = crit, start = start)
   current.chol = NULL
   method = match.arg(method)
   search = 1:p
   if (!is.null(first.entry) && isTRUE(find.first)) {
      stop("only one of first.entry and find.first should be specified")
   }
   if (!is.null(first.entry)) {
      if(first.entry> p){
         stop("first.entry should be less than the number of dimensions")
      }
      if(first.entry<0){
         stop("first.entry should be positive integer")
      }
   }
   if (is.null(first.entry)){
      order = which.min(diag(S))
   } else {
      order = first.entry
   }
   search = search[-order]
   while(length(search)>0){
      if(isFALSE(find.first)){
         if (method == "det"){
            find.min = find.min.det(S, order, search)
         } else {
            find.min = find.min.chol(S, order, search, current.chol)
            current.chol = find.min$current.chol
         }
      } else {
         if (method == "det"){
            find.min = find.min.det.ext(S, order, search)
         } else {
            find.min = find.min.chol.ext(S, order, search, current.chol)
            current.chol = find.min$current.chol
         }
      }
      order = find.min$current
      search = find.min$search
   }
   return(order)
}

find.min.det <- function(S, current, search){
      current.min = c()
      ind = 1
      for (i in search){
            if (length(current) != 1){
               denom = det(S[current, current])
            } else {
               denom = S[current,current]
            }
            current.min[ind] = det(S[c(current,i), c(current,i)]) / denom
            ind = ind + 1
      }
      min.ind = which.min(current.min)
      return(list(current = c(current,search[min.ind]), search = search[-min.ind]))
}

find.min.chol <- function(S, current, search, current.chol = NULL){
   current.min = c()
   ind = 1
   if (length(current) == 1){
      current.chol = sqrt(S[current,current])
   }
   #   else {
 #     current.chol = t(chol(S[current,current]))
 #  }
   L.store = array(0, c(length(current) + 1, length(current) + 1, length(search) ))
   for (i in search){
      est.min = find.chol.row(S[c(current,i), c(current,i)], current.chol)
      current.min[ind] = est.min$diag.sum #find.chol.row(S[c(current,i), c(current,i)], current.chol)
      L.store[,,ind] = est.min$current.chol
      ind = ind + 1
   }
   min.ind = which.min(current.min)
   current.chol = L.store[,,min.ind]
   return(list(current = c(current,search[min.ind]), search = search[-min.ind], current.chol = current.chol))
}


find.chol.row <- function(S, current.chol) {
   if (length(current.chol) == 1)
   {
      p = 1
   } else {
      p = dim(current.chol)[1]
   }
   i = p + 1
   L = matrix(0, i, i)
   L[1:p, 1:p] = current.chol
   L[i,1] = S[i,1]/ L[1,1]
   for (j in 2:i){
      sum = 0
      for (k in 1:(j - 1)) {
         sum = sum + L[i, k] * L[j,k]
      }
      if(i  == j){
         L[i,j] = sqrt(S[i,i] - sum)
      } else {
         L[i,j] = (1/L[j,j]) * (S[i,j] - sum)
      }
   }
   return(list(diag.sum = sqrt(sum(diag(L)^2)), current.chol = L))
}

find.min.chol.ext <- function(S, current, search, current.chol = NULL){
   current.min = c()
   ind = 1
   if (length(current) == 1){
      current.chol = sqrt(S[current,current])
   }
   #else {
   #   current.chol = t(chol(S[current,current]))
   #}
   L.store = array(0, c(length(current) + 1, length(current) + 1, length(search) ))
   for (i in search){
      est.min = find.chol.row(S[c(current,i), c(current,i)], current.chol)
      current.min[ind] = est.min$diag.sum #find.chol.row(S[c(current,i), c(current,i)], current.chol)
      L.store[,,ind] = est.min$current.chol
      ind = ind + 1
   }

   min.ind = which.min(current.min)

   if(sqrt(sum(diag(chol(S[c(current,search[min.ind]),c(current,search[min.ind])]))^2)) >
      sqrt(sum(diag(chol(S[c(search[min.ind], current),c(search[min.ind], current)]))^2))) {

      current = c(search[min.ind], current)
   }
   #      }
   else {
      current = c(current,search[min.ind])
   }
   search = search[-min.ind]
   current.chol = L.store[,,min.ind]
   return(list(current = current, search = search, current.chol = current.chol))
}


find.min.det.ext <- function(S, current, search){
      current.min = c()
      ind = 1
      for (i in search){
            if (length(current) != 1){
               denom = det(S[current, current])
            } else {
               denom = S[current,current]
            }
            current.min[ind] = det(S[c(current,i), c(current,i)]) / denom
            ind = ind + 1
      }
      min.ind = which.min(current.min)
      curr.len = length(current)
       if(sqrt(sum(diag(chol(S[c(current,search[min.ind]),c(current,search[min.ind])]))^2)) >
           sqrt(sum(diag(chol(S[c(search[min.ind], current),c(search[min.ind], current)]))^2))) {

               current = c(search[min.ind], current)
             }
#      }
      else {
         current = c(current,search[min.ind])
      }
      search = search[-min.ind]
      return(list(current = current, search = search))
}


###############
### Fit DAG using topological ordering
###############
#' Infer  DAG using topological ordering. Adopted from EqVarDAG package.
#' @param X: data in n x p matrix
#' @param TO: topological ordering
#' @param alpha: desired selection significance level
#' @param mtd: methods for learning DAG from topological orderings.
#'  "ztest": (p<n) [Multiple Testing and Error Control in Gaussian Graphical Model Selection. Drton and Perlman.2007]
#'  "rls": (p<n) fit recursive least squares using ggm package and threshold the regression coefs
#'  "chol": (p<n) perform cholesky decomposition and threshold the regression coefs
#'  "dlasso": debiased lasso (default with FCD=True and precmtd="sqrtlasso");
#'   "lasso": lasso with fixed lambda from [Penalized likelihood methods for estimation of sparse high-dimensional directed acyclic graphs. Shojaie and Michailidis. 2010];
#'   "adalasso": adaptive lasso with fixed lambda from [Shojaie and Michailidis. 2010];
#'   "cvlasso": cross-validated lasso from glmnet;
#'    "scallasso": scaled lasso.
#' @param threshold: only used in rls and chol. the hard threshold level.
#' @param FCD: only used in debiased lasso,  the FCD procedure [False Discovery Rate Control via Debiased Lasso. Javanmard and Montanari. 2018]
#' or use individual tests to select support.
#' @param precmtd: only used in debiased lasso, how to compute debiasing matrix
#'               "cv": node-wise lasso w/ joint 10 fold cv
#'               "sqrtlasso": square-root lasso(no tune, default)
#' @return Adjacency matrix with ADJ[i,j]!=0 iff i->j
DAG_from_Ordering<-function(X, order, mtd="ztest", alpha=0.05,
                            threshold=1e-1, FCD=NULL, precmtd=NULL){
   n=dim(X)[1]
   p=dim(X)[2]
   if (p!=length(order)){stop("length mismatch")}
   if (mtd=="ztest"){
      # sidak
      C=cor(X)
      adj=matrix(0,p,p)
      for (i in 2:p){
         u=order[i]
         for (j in 1:(i-1)){
            v = order[j]
            s = setdiff(order[seq(i-1)],v)
            pval = 1-(2*pnorm(abs(pcalg::zStat(u,v,s,C=C,n=n)))-1)^(p*(p-1)/2)
            adj[v,u]=ifelse(pval<alpha,1,0)
         }
      }
      return((adj!=0)*1)
   }
   if (mtd=="chol"){
      Sigma=cov(X)
      B = solve(t(chol(Sigma[order,order]))[order(order),order(order)])
      gm = diag(p)-B%*%diag(1/diag(B))
      # if (isTRUE(unbias.chol)){
      #    N = n - 1
      #    gm = getL.tilda(gm, N)
      # }
      return((gm*(abs(gm)>threshold)!=0)*1)
   }
   if (mtd=="rls"){
      gm = upper.tri(matrix(0,p,p))[order(order),order(order)]
      colnames(gm)=rownames(gm)=colnames(X)
      return(abs(t(ggm::fitDag(gm,cov(X),dim(X)[1])$Ahat))-diag(p)>threshold)
   } else {
      # dblasso
      if (is.null(FCD)){FCD="T"}
      if (is.null(precmtd)){precmtd="sqrtlasso"}
      gm = matrix(0,p,p)
      gm[order[1],order[2]]=anova(lm(X[,order[2]]~X[,order[1]]))$`Pr(>F)`[1]<alpha
      if(p==2){return(gm * 1)}
      for (i in 3:p){
         gm[order[1:(i-1)],order[i]]=
            EqVarDAG::vselect(X[,order[1:i-1]],X[,order[i]],alpha=alpha,p_total = p,
                              selmtd = mtd,FCD = FCD,precmtd = precmtd)$selected
      }
      return((gm!=0)*1)
   }
}

getL.tilda <- function(L, n) {
   p = nrow(L)
   Lambda = diag(n + 1- (1:p), nrow = p)
   return(L %*% Lambda)
}



