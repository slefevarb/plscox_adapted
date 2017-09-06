
#########
#   
#     Function my_plsRcox() 
#       . adapted from package plsRcox version 1.7.2
#       . Bertrand F, Philippe B, Meyer N, Maumy-Bertrand M. plsRcox: Cox- Models in a high dimensional setting in R, UseR 2014, Los Angeles. USA: 2014.
# 
# Required packages 
#   . Survival
#
# Inputs
#   . X - The matrix of standardized predictors dim(n x p)
#   . time1 - For right censored data, this is the follow up time. 
#             For interval-censored data, the first argument is the starting time for the interval. (see Survival::Surv())
#   . time2 - For interval-censored data, ending time of the interval. 
#             Intervals are assumed to be open on the left and closed on the right, (start, end]. (see Survival::Surv())
#   . status - A numeric vector of 0 (alive/healthy) and 1 (dead/ill) (see Survival::Surv())
#   . adjust - A vector (1 covariate) or a matrix (>1 covariate) of confounders to adjust for
#            => The defaut is no adjustment (equivalent to plsRcox)
#   . ncomp - The number of PLS component to compute
#
# Outputs
#
#   . tt - The matrix of ncomp PLS components dim(n x ncomp)
#   . pp - The matrix of Pearson correlation coefficients dim(p x ncomp)
#   . cox - Results from the fitted Cox proportional hazards regression model
#
#########

my_plsRcox <- function(X, time1, time2, status, adjust=NULL, ncomp)
  
{
  
  require(survival)
  
  names_var = colnames(X)
  n = nrow(X)
  p = ncol(X)
  
  t_mat = matrix(nrow = n, ncol = ncomp, 
                 dimnames = list(c(1:n), c(1:ncomp)))
  p_mat = matrix(nrow = p, ncol = ncomp)
  
  names_var = colnames(X)
  p = ncol(X)
  n = nrow(X)
  
  
  if(is.null(dim(adjust))) { adjustment = 'adjust' } else { adjustment = paste(names(adjust), collapse = "+") }
  if(is.null(adjust)) { adjustment = 0 }
  
  
  ##################
  # LOOP
  ##################
  
  for (comp in seq(ncomp)){
    
    SurvObj = Surv(time = time1, time2 = time2, event = status)
    #========================================================================
    # COMPUTATION OF THE FIRST PLS COMPONENT h=1
    #======================================================================== 
    if (comp == 1) {   
      
      #1) Run separate cox regressions of y on each x_j
      
      a = rep(NA,ncol(X)) # a - vector of p regression coef a_j
      for (xi in 1:ncol(X)){
        a[xi] = coxph(as.formula(paste("SurvObj~X[,xi]+", adjustment)), 
                      data=data.frame(cbind(X,adjust)))$coef[1]
      }
      
      w = a/sqrt(sum(a^2)) # w - normalized vector of p regression coef w_j
      
      
      #2) Compute the linear combination of x_j weighted by w_j
      
      t = as.matrix(X)%*%w # tt - the PLS component
      t_mat[,comp] = t
      
      
      #3) Run separate linear regressions of each x_j on t
      
      res = matrix(nrow = nrow(X), ncol = ncol(X))
      coef = rep(NA, ncol(X))
      
      for (xi in 1:ncol(X)){
        mod = lm(X[,xi]~scale(t))
        res[,xi] = mod$residuals # res is the matrix of p residuals x_hj
        coef[xi] = mod$coef[length(mod$coef)] 
      }
      
      p_mat[,comp] = coef # p_mat - matrix of vectors of correlation p_hj
    }
    
    
    #========================================================================
    # COMPUTATION OF THE NEXT PLS COMPONENT h=(2,...,ncomp)
    #========================================================================
    if (comp > 1){ 
      
      #1) Run separate cox regressions of y on each x_j and t1,...,t_h-1
      
      a = rep(NA,ncol(X))
      
      for (xi in 1:ncol(X)){
        a[xi] = coxph(as.formula(
          paste("SurvObj~X[,xi]+", 
                paste(c(names(data.frame(t_mat))[1:(comp-1)],
                        adjustment), collapse = "+"))), 
          data=data.frame(cbind(t_mat,adjust)))$coef[1]
      }
      
      w = a/sqrt(sum(a^2))
      
      #2) Compute the linear combination of x_j weighted by w_j
      
      t = as.matrix(res)%*%w
      t_mat[,comp] = t
      
      #3) Run separate linear regressions of each x_j on t
      
      res = matrix(nrow = nrow(X), ncol = ncol(X))
      coef = rep(NA, ncol(X))
      
      for (xi in 1:ncol(X)){
        mod = lm(as.formula(
          paste("X[,xi]~",
                paste(names(data.frame(t_mat))[1:comp], collapse = "+"))), 
          data=data.frame(scale(t_mat)))
        res[,xi] = mod$residuals 
        coef[xi] = mod$coef[length(mod$coef)]
      }
      
      p_mat[,comp] = coef 
    }
  }
  
  ##################
  
  mod.cox = coxph(as.formula(
    paste("SurvObj~", 
          paste(names(data.frame(t_mat)), collapse = "+"))), 
    data=data.frame(scale(t_mat)))
  rownames(p_mat)=names_var
  object = list("tt" = t_mat, "pp" = p_mat, "cox"=mod.cox)
  return(object)
  
}


