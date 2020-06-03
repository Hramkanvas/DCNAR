t <- 100
s <- 3
m <- 4
n <- 10
x_hist <- c(1,2,3)
#для примера возьмем произвольный вектор параметров theta
#предположим, что вектор psi(u) =(1, u1, u2, u3)
theta <- c(-0.5, -0.3, -0.2, -0.4)
psi <- function(x){
  return(c(1, x[1], x[2], x[3]))
}

update_hist <- function(x, add_var){
  new_x <- x[-1]
  new_x <- append(new_x, add_var)
  return(new_x)
}

psi <- function(x){
  return(c(1, x[1], x[2], x[3]))
}
sigmoid <- function(x){
  return(1 / (1 + exp(-x)))
}
prob <- function(x, theta, psi){
  return(sigmoid(sum(theta*psi(x))))
}

#t - size
#n - value range
#x - history
#' @export
BiCNAR <- function(x, t, n, theta, psi){
  result = x
  for(i in seq(0,t)){
    add_var = rbinom(1,n, prob(x, theta, psi))
    result = append(result, add_var)
    x = update_hist(x, add_var)
  }
  return(result)
}



# PCNAR
lam_pcnar <- function(x, theta, psi){
  return(exp(sum(theta*psi(x))))
}

#' @export
PCNAR <- function(x, t, theta, psi){
  result = x
  for(i in seq(0,t)){
    add_var = rpois(1, lam_pcnar(x, theta, psi))
    result = append(result, add_var)
    x = update_hist(x, add_var)
  }
  return(result)
}

#INGARCH

lam_ingarch <- function(alpha_0, alpha, beta, x, lam, p, q){
  first_sum = 0
  second_sum = 0
  result = alpha_0
  for(i in seq(1,p)){
    first_sum = first_sum + alpha[i]*x[length(x) - i + 1]
  }
  for(j in seq(1,q)){
    second_sum = second_sum + beta[j]*lam[length(lam) - j + 1]
  }
  result = result + first_sum + second_sum
  return(result)
}

#' @export
INGARCH <- function(x, t, p, q, alpha_0, alpha, beta, lam){
  result = x
  for( i in seq(1, t)){
    add_lam = lam_ingarch(alpha_0, alpha, beta, x, lam, p, q)
    add_var = rpois(1, add_lam)
    result = append(result, add_var)
    x = update_hist(x, add_var)
    lam = update_hist(lam, add_lam)
  }
  return(result)
}

# SETINAR
#' @export
SETINAR <- function(x, t, lam, R, alpha_1, alpha_2){
  result = x
  for(i in seq(0, t)){
    z = rpois(1,lam)
    sum = 0
    if(x <= R){
      for( i in seq(x)){
        sum = sum + rbinom(1,1,alpha_1)
      }
    } else {
      for(i in seq(x)){
        sum = sum + rbinom(1,1,alpha_2)
      }
    }
    x = z + sum
    result = append(result, x)
  }
  return(result)
}

#BCNAR

#' @export
BCNAR <- function(x, t, theta, psi){
  result = x
  for(i in seq(0,t)){
    add_var = rbinom(1,1, prob(x, theta, psi))
    result = append(result, add_var)
    x = update_hist(x, add_var)
  }
  return(result)
}


# JL

discreteGenerator <- function(probVector){
  drop = runif(1)
  index = 1
  sum = 0
  prev = 0
  result = 0
  for(i in probVector){
    prev = sum
    sum = sum + i
    if (prev <= drop && drop <= sum) {
      result = index
    }
    index = index + 1
  }
  return(result)
}
#' @export
JL <- function(pi, lambda, ro, t, x_hist){
  result = x_hist
  for(i in seq(1,t)){
    ksi = discreteGenerator(pi)
    nu = discreteGenerator(lambda)
    mu = discreteGenerator(c(1-ro, ro)) - 1
    add_var = mu*x_hist[length(x_hist) + 1 - nu] + (1 - mu)*ksi
    result = append(result, add_var)
    x_hist = update_hist(x_hist, add_var)
  }
  return(result)
}

#Markcov chain

#' @export
MC <- function(P, t, x_hist) {

  numStates <- nrow(P)
  states <- numeric(t)
  states[1] <- x_hist[length(x_hist)]

  for(i in 2:t) {
    p  <- P[states[i-1], ]
    states[i] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}

#PCNAR Estimation

count_ar <- function(vec, item, test, s){
  counter = c(0,0)
  ctr = s + 1
  for(i in vec){
    if (all(i == item)){
      counter[1] =  counter[1] + 1
      counter[2] =  counter[2] + test[ctr]
    }
    ctr = ctr + 1
  }
  return(counter)
}
uniq <- function(vec, test, s){
  result = list()
  index <- 1
  for( i in vec){
    if (count_ar(result, i, test, s)[1] == 0){
      result[[index]] <- i
      index = index + 1
    }
  }
  return(result)
}

#[array([1, 2, 3, 4]), [2, 2]] - так выглядит элемент, 1элемент второго массива - #частота, 2элемент - сумма x
get_lambda <- function(info, m){
  lmbda = c()

  for(i in 1:m){
    if (1.*info[[i]][[2]][2]/info[[i]][[2]][1] == 0){
      lmbda = append(lmbda, 10**(-6))
    }
    else{

      lmbda = append(lmbda, 1.*info[[i]][[2]][2]/info[[i]][[2]][1])
    }
  }
  return(lmbda)
}

#' @export
PCNAR_Estimation <- function(global_results, s, m, psi ){
  data = list()
  jndex <- 1
  for(test in global_results){
    length_n_test = length(test) - s
    n_test = vector(mode = "list", length = length_n_test)
    for(i in 1:length_n_test){
      n_test[[i]] <- test[i:(i+s-1)]
    }

    unique_test = uniq(n_test, test, s)
    test_with_freq <- vector(mode = "list", length = length(unique_test))
    index <- 1
    for (i in unique_test){
      test_with_freq[[index]] <- list(i, count_ar(n_test, i, test, s))
      index = index + 1
    }
    id <- order(sapply(test_with_freq,function(i)-i[[2]][1]))
    test_with_freq = test_with_freq[id]
    lmbda = get_lambda(test_with_freq, m)
    psi_matr <- NULL
    for (i in test_with_freq[1:m]){
      psi_matr <- cbind(psi_matr, psi(i[[1]]))
    }
    element <- log(lmbda) %*% solve(psi_matr)
    data[[jndex]] <- element
    jndex = jndex + 1
  }
  return(data)
}

#BiCNAR Extimation

reverseSigmoid <- function(x){
  result <- c()
  for(i in x){
    result = append(result,-log(1./i-1))
  }
  return(result)
}

#' @export
BiCNAR_Estimation <- function(global_results, s, m, psi,n ){
  data = list()
  jndex <- 1
  for(test in global_results){
    length_n_test = length(test) - s
    n_test = vector(mode = "list", length = length_n_test)
    for(i in 1:length_n_test){
      n_test[[i]] <- test[i:(i+s-1)]
    }

    unique_test = uniq(n_test, test, s)
    test_with_freq <- vector(mode = "list", length = length(unique_test))
    index <- 1
    for (i in unique_test){
      test_with_freq[[index]] <- list(i, count_ar(n_test, i, test, s))
      index = index + 1
    }
    id <- order(sapply(test_with_freq,function(i)-i[[2]][1]))
    test_with_freq = test_with_freq[id]
    lmbda = get_lambda(test_with_freq, m)/n
    psi_matr <- NULL
    for (i in test_with_freq[1:m]){
      psi_matr <- cbind(psi_matr, psi(i[[1]]))
    }
    element <- reverseSigmoid(lmbda) %*% solve(psi_matr)
    data[[jndex]] <- element
    jndex = jndex + 1
  }
  return(data)
}


