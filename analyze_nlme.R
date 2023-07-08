library(tidyverse)
library(nlme)
library(broom)
library(gtools)



###########################################################
### Find starting points for nlme by fitting individual nls
###########################################################

fit_nls <- function(df) {
     try(
          nls( skill ~ SSlogis( progress, asym, xmid, scal), data=df),
          silent = TRUE 
     )
}
tidy_up <- function(model){
     try(
          broom::tidy(model),
          silent = TRUE
     )
}
is_tbl <- function(alist){
     "tbl" %in% alist
}
find_starts <- function (df) {
    starts <- df %>%
        filter(session != 4) %>%
        pivot_longer(
            cols = c(integrated_skill,ls_skill,mc_skill,de_skill),
            names_to = "task",
            values_to = "skill") %>%
        select(group,id,task,progress,skill) %>%
        group_by(group, task) %>%
        nest() %>%
        mutate(models = map( map(data, fit_nls), tidy_up))%>%
        filter(map_vec(map(models,class), is_tbl)) %>%
        unnest(cols=c(models)) %>%
        ungroup() %>%
        select(group,task,term,estimate) %>%
        pivot_wider(names_from=term, values_from=estimate)
    return (starts)
}
find_bounds <- function(starts) {
    bounds <- starts %>%
        group_by(task) %>%
        summarise(
            asym_low = min(asym)*0.8,
            asym_high = max(asym)*1.2,
            scal_low = min(scal)*0.8,
            scal_high = max(scal)*1.2,
            xmid_start = mean(xmid)
        )
    return (bounds)
}


###############################################
### Functions for fitting a single nlme
###############################################

new_data <- function(groups) {
     group <- list_c(map(groups, (function (name){rep(name,times=30)})))
     progress <- rep(1:30 / 30, times= length(groups))
     return (data.frame(group=group,progress=progress))
}

fit_nlme <- function(df, skill, start) {
    control <- nlmeControl(
        opt="nlminb",
        pnlsMaxIter = 7,
        msMaxIter=100,
        maxIter=200,
        minScale=0.0001,
        pnlsTol=.01,
        tolerance=1e-8,
        rel.tol=1e-8)
    formula <- as.formula(paste0(skill,"~ SSlogis(progress, asym, xmid, scal)"))
    model <- do.call( "nlme", list( # Ensures the formula is evaluated in the right environment
        model = formula,
        data=df, # Ensures df is not printed in the output
        fixed = c(asym  + scal ~  group, xmid  ~ 1),
        random = asym  ~ 1 | id,
        start = start,
        control = control,
        na.action = na.omit
        ))
    fits <- new_data(unique(df$group))
    fits$predictions <- predict(object=model, newdata=fits, level=0)
    return(list(model=model, fits=fits))
}

extract_residuals <- function(fitted,skill) {
    residual_data <- data.frame(
        progress = df %>% 
            filter(session!=4) %>%
            select(progress, .data[[skill]]) %>%
            na.omit() %>% 
            select(progress), 
        res_fixed = fitted$model$residuals[,'fixed'],
        res_rand = fitted$model$residuals[,'id']) %>%
        mutate(res = res_fixed + res_rand)
    return(residual_data)
}


#############################################################
## Permutation testing 
#############################################################


# Samples contains a named list of previously extracted
# Model coefficients.  Append new coefficients from model
extract_coeff <- function (model, samples) {
     new_samples <- c()
     if (length(samples) == 0) {
          return(map(model$coefficients$fixed, c))
     }
     for (name in names(model$coefficients$fixed)) {
          new_samples[[name]] <- append(model$coefficients$fixed[[name]], samples[[name]])
     }
     return (new_samples)
}


# Function to reformat fixed coeffs with the 2U1D_unlocked group as the
# intercept.
change_intercept <- function(coeffs) {
    adj_coeffs <- c()
    adj_coeffs[["asym.(Intercept)"]] <- 
        coeffs[["asym.(Intercept)"]] + coeffs[["asym.group2U1D_unlocked"]]
    adj_coeffs[["scal.(Intercept)"]] <-
        coeffs[["scal.(Intercept)"]] + coeffs[["scal.group2U1D_unlocked"]]
    adj_coeffs[["asym.group1U1D"]] <- - coeffs[["asym.group2U1D_unlocked"]]
    adj_coeffs[["scal.group1U1D"]] <- - coeffs[["scal.group2U1D_unlocked"]]
    adj_coeffs[["xmid"]] <- coeffs[["xmid"]]
    for (param in names(coeffs)) {
        if (grepl(".group", param) & !grepl("unlocked",param)) {
            if (grepl("asym", param)) {
                adj_coeffs[[param]] <- 
                     coeffs[[param]] - coeffs[["asym.group2U1D_unlocked"]] 
            } else{
                adj_coeffs[[param]] <-
                    coeffs[[param]] - coeffs[["scal.group2U1D_unlocked"]] 
            }

        }
    }
    return(adj_coeffs)
}

# Performs a permutation test on the groups
dist_params <- function(df, n, skill) {
    samples <- c() # Initialize empty vector of samples for each parameter
    for (x in 1:n) {
        try({
            # Reassign each subject to a random group
            scrambled_df <- df %>%
                nest(.by=c(id,group)) %>%
                mutate(group = sample(group)) %>%
                unnest(cols=c(data))
            # Attempt to fit the scrambled assignments
            scrambled_fit <- fit_nlme(
                scrambled_df, 
                skill,
                start = c(
                    asym= rep(1.3192336,times=length(unique(df$group))),
                    scal= rep(0.1537297,times=length(unique(df$group))),
                    xmid = 0.10) )
            # print("fit success")
            # print(scrambled_fit$model$coefficients$fixed)
            samples <- extract_coeff(scrambled_fit$model,samples)},
            silent = FALSE
        )
    }
    # How many models were successfully fit
    print(paste0("Model fit success rate: ", length(samples[[1]]) / n))
    if (any(grepl("2U1D_unlocked",names(samples)))) {
        samples <- change_intercept(samples)
    }
    return(samples)
}

params_cdf <- function (x, samples) {
    return( 
        integrate(
            lower = abs(x),
            upper = Inf,
            f = approxfun(density(abs(samples), from=0), yright=0))
    )
}

permutation_test <- function(df, task, true_params) {
    samples <- dist_params(df, 1000, task)
    pvals <- c()
    for (param in names(samples)) {
        if (grepl(".group", param, fixed=TRUE)) {
            pvals[[param]] <- params_cdf(
                true_params[[param]], 
                samples[[param]])$value
        }
    }
    return ( list("pvals"=pvals, "samples"=samples))
}



######################################################
#### Code for autofitting with poor starting guess
######################################################
# This code is not used in the notebooks
# In its place, a simpler model was used.
# The more complex models required this approach
# But it wasn't effective enough to use for now.
# In case further attempts are made, the code is here.

# Similar to above, but expects coeffs as a list already
extract_coeff_auto <- function(coeffs, samples) {
     new_samples <- c()
     if (length(samples) == 0) {
          return(map(coeffs, c))
     }
     for (name in names(coeffs)) {
          new_samples[[name]] <- append(coeffs[[name]], samples[[name]])
     }
     return (new_samples)  
}

# Generate a named vector of perturbed starting points
rstart <- function(df) {
    groups <- unique(df$group)
    start_asym <- 1.3192336
    start_scal <- 0.1537297
    start_xmid <- 0.10
    #  start_asym <- 1.84165007
    #  start_scal <- 0.515466
    #  start_xmid <- 0.1790322
    #  start_asym <- 1.773407
    #  start_scal <- 0.3242479
    #  start_xmid <- 0.1681572

     asym <- rep(start_asym + rnorm(1,0,0.2), times=length(groups))
     names(asym) <- groups
     scal <- rep(start_scal + rnorm(1,0,0.4766), times=length(groups))
     names(scal) <- groups
     xmid <- start_xmid + rnorm(1,0,0.001)

     start <- c( asym = asym, scal=scal, xmid = xmid) 
     return(start) 
}

# Attempt ten fits using randomized starting points
rfit <- function (df, skill) {
    # Model has not yet been fit
    succeeded <- FALSE
    # We have attempted zero times
    n <- 0
    # Until success or ten attempts
    while (!succeeded && n < 2) {
        try( # Perform a random fit, catch errors and continue
            {
            model <- fit_nlme(df,skill,rstart(df))
            succeeded <- TRUE
            },
            silent=TRUE # fail silently
        )
        n <- n + 1 # increment attempt counter
    }
    # Model fit successfully
    if (succeeded) {
        return (model)
    }
    # Failed to fit model, return FALSE
    return (FALSE)
}

### Functions used in the autofitting algorithm below
# Generate the form of a fresh guess, as required by nlme package
fresh_start <- function(asym_guess, scal_guess, xmid_guess) {
     c(
          'asym.(Intercept)'=asym_guess, 
          'scal.(Intercept)'=scal_guess,
          'xmid'=xmid_guess)  
}

# Add a new group to start
add_to_start <- function(current_guess, new_group) {
     start <- current_guess
     start[paste0('asym.group',new_group)] <- 0
     start[paste0('scal.group',new_group)] <- 0
     return(start)
}

# Remove a list of groups from a start
remove_from_start <- function(current_guess, bad_groups) {
     start <- current_guess
     start_groups <- groups_in_start(start)
     for (group in bad_groups) {
          if (group %in% start_groups) {
               start <- start[ !names(start) %in% c( paste0('asym.group',group), paste0('scal.group',group)) ]
          }
     }
     return (start)
}

# Return a list of the groups in start
groups_in_start <- function(current_guess){
     if (length(current_guess)==3){
          return (list())
     }
     return( 
          str_subset(
               unique(
                    str_split(
                         names(current_guess), 
                         ".group", 
                         simplify=TRUE)[,2]
               ),
               ".+"
          )
     )
}
# Handler for errors in autofitting
handle_errors <- function(group_idx, perms, perm_idx, res) {
        return( list(
            "skip_to"=group_idx, 
            "perm_prev"=perms[perm_idx,], 
            "failed"=TRUE,
            "current_guess"=res$current_guess))
}

# Takes poor starting guesses and iteratively improves them until convergence
# 1. Constructs all permutations of the groups
# 2. For each permutation, try fitting the first two groups only
# 3. If failure, perturb the starting guesses slightly and try again
# 4. If success, start at results of first fit and add a third group.
# 5. Continue until all four groups are fit
# 6. If failure, try a different permutation of the groups and repeat
# Perform auto_fitting
auto_fit <- function(df, task, asym_guess, scal_guess, xmid_guess) {

    groups <- unique(df$group) 
    # Matrix of all possible permutations of group names
    # Perms change last element to first.  
    perms <- permutations(length(groups),length(groups),levels(groups))
    # Index of last group that failed to fit
    skip_to <- 1
    # Permutation before the current
    perm_prev <- perms[1,]
    # Guess for parameter values
    current_guess <- fresh_start(asym_guess, scal_guess, xmid_guess)
    # Result object to return
    res <- list("current_guess"=current_guess, "perm_prev"=perm_prev, 'skip_to'=skip_to)

    # Iterate over all possible permutations
    for (perm_idx in 1:length(perms[,1])) {
        # If all parameters have been fit
        if (length(res$current_guess) == 2*length(groups)+1) {
            return(res$current_guess)
        }
        # Detect first element of new perm that is different from prev
        diff_idx <- which.min(res$perm_prev==perms[perm_idx,])
        
        # If the first element is different, need a fresh start
        if (diff_idx == 1) {
            res$current_guess <- fresh_start(asym_guess, scal_guess, xmid_guess)
        } else { 
            # If new perm is the same up to the last failure, get the next one
            if (diff_idx > res$skip_to) { 
                next 
            } 
            # If new perm only matches the first few of prev, remove guesses
            # for groups that are no longer needed
            else { 
                res$current_guess <- 
                        remove_from_start(
                            res$current_guess,
                            res$perm_prev [diff_idx : length(groups) ]
                            )
            }
        }
        # Reset failed flag.  The new perm hasn't failed the fitting test
        res$failed <- FALSE
        # Iterate over the groups in the new perm, adding one at a time
        for (group_idx in max(c(2,diff_idx)):length(groups)) {
            # If this perm can't be fit, skip to the next
            if (res$failed) {
                next
            }
            failed <- FALSE ############### To erase

            # Try a fit, catch failures and continue the program
            res<- tryCatch(
                {
                # If this is the first attempt for this perm
                # Randomize the starting point.
                if (group_idx == 2) {
                        fit <- rfit(
                            df %>% filter(session != 4, group %in% perms[perm_idx,1:group_idx]),
                            task)
                } else{ # Use the parameters from the previous successful fit
                # Add the new group to the list of starting parameters
                start <-add_to_start(res$current_guess, perms[perm_idx,group_idx]);
                # Perform the fit
                fit <- do.call( "fit_nlme", list(
                        substitute(df %>% filter(session != 4, group %in% perms[perm_idx,1:group_idx])),
                        skill=task,
                        start = start
                ))};
                # Only executed if fit is successful
                # Return the new current guesses
                list("current_guess" = fit$model$coefficients$fixed, "failed"=FALSE)
                },
                # If there was an error, call the handler
                error = function(cond) {handle_errors(group_idx, perms, perm_idx, res)},
                # If there was a warning, call the handler
                warning = function(cond) {handle_errors(group_idx, perms, perm_idx, res)}
            )
        }
    }
    # If all permutations have been exhausted,
    # The algorithm failed.
    return(FALSE)
}