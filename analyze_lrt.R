library(Hmisc)
library(tidyverse)




compute_lik <- function (n, xbar) {
    xbar^(n*xbar) * (1- xbar)^(n*(1-xbar))
}

compute_lik_rat <- function (n,xbar, nj, xbarj) {
    numerator <- compute_lik(n, xbar)
    denominator <- reduce(map2(nj, xbarj, function(nn, bar){compute_lik(nn,bar)}),`*`)
    return( numerator / denominator)
}

compute_np <- function (df, task) {
    p <- mean(df[[task]], na.rm=TRUE)
    n <- nrow(df %>% na.omit())
    nj <- df %>% 
        na.omit() %>%
        group_by(group) %>%
        summarise(nj = n()) %>%
        pull(nj)
    return(list("p"=p, "n"=n, "nj"=nj))
}

dist_lik_rat <- function (n, p, nj) {
    # Simulate 10000 values of the likelihood ratio for this p
    # Assume the null hypothesis is true
    liks <- c()
    for (x in 1:20000) {
        samples <- rbinom(n, 1, prob=p)
        group_samples <- partition.vector(samples, sep=nj)
        xbarj <- map( group_samples, mean)
        xbar <- mean(samples)
        liks <- append(liks, compute_lik_rat(n,xbar,nj,xbarj))
    }
    return (liks)
}

lik_cdf <- function (x, liks) {
    total <- length(liks)
    cdf <- data.frame(liks) %>% 
            filter(liks <= x) %>%
            group_by(liks) %>%
            summarise(probs = n()/total) %>%
            pull(probs) %>%
            sum
    return( cdf)
}

glrt <- function(df, task) {
    data.in <- df %>% select(all_of(c(task,"group")))
    np <- compute_np(data.in, task)
    liks <- dist_lik_rat(np$n, np$p, np$nj)
    xbar <- np$p
    xbarj <- data.in %>%
        group_by(group) %>%
        summarise(xbarj = mean(.data[[task]],na.rm=TRUE)) %>%
        pull(xbarj)
    lik <- compute_lik_rat(np$n, xbar, np$nj, xbarj)
    return(lik_cdf(lik,liks))
}