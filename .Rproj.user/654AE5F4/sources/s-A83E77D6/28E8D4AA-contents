#' algoleverage Function
#'
#' This function implements algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows.
#' @param method The method for algorithmic leveraging, choose from "unif","leverage".
#' @param X dataset X
#' @param Y dataset Y
#' @param n sample size of data set
#' @param r size of a randomly selected subset
#' @param Pi sampling probability
#' @param replication replication, how many draws of subsample
#' @param beta_ols the linear regression coefficient from ols
#'
#' @return  return the absolute error between coefficient from leveraging algorithm and coefficient from ols for each draw of subsample
#'
#' @export
#' @examples
#' set.seed(65206520)
#' n=500
#' X=rt(n, 6)
#' eps=rnorm(n,0,1)
#' Y=-1*X+eps
#' beta_ols=lm(Y ~ 0 + X )$coefficients
#' algoleverage("unif",X,Y,n,r=10,replication=500,beta_ols=beta_ols)


algoleverage <- function(method,X,Y,n,r,Pi=NULL,replication=500,beta_ols) {
  if (method=="unif"){
    Pi=rep(1/n,n)
    res=unif_fun(X,Y,n,r,Pi,replication,beta_ols)

  }
  else if (method=="leverage"){
    res=blev_fun(X,Y,n,r,Pi,replication,beta_ols)

  }
  else
    print("unsupported method type")

  return(res)
}


unif_fun <- function(X,Y,n,r,Pi,replication=500,beta_ols) {
  beta_unif=rep(NA,replication)
  for (i in 1:replication){
    index_1=sample(c(1:n), size = r, replace = T,prob=Pi)
    #Phi_1=diag(Pi_1[index_1])
    X_star_1=X[index_1]
    Y_star_1=Y[index_1]
    Phi_1=1/ (Pi[index_1]^{1/2})
    model_1 <- lm(Y_star_1 ~ 0 + X_star_1 , weights=Phi_1)

    beta_hat_1=model_1$coefficients
    beta_unif[i]=abs(beta_hat_1-beta_ols)
  }
  return(beta_unif)
}

blev_fun <- function(X,Y,n,r,Pi,replication=500,beta_ols) {
  beta_unif=rep(NA,replication)
  for (i in 1:replication){
    index_2=sample(c(1:n), size = r, replace = T,prob=Pi)

    X_star_2=X[index_2]
    Y_star_2=Y[index_2]
    Phi_2=1/ (Pi[index_2]^{1/2})
    #Phi_2=1/ (Pi_2[index_2]^{1/2})
    model_2 <- lm(Y_star_2 ~ 0 + X_star_2 , weights=Phi_2)

    beta_hat_2=model_2$coefficients
    beta_unif[i]=abs(beta_hat_2-beta_ols)
  }
  return(beta_unif)
}



