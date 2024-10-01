# The WrappedmlegpGP class
#
# This is a wrapper class used to contain a GP instance from the mlegp library and provide
# the interface expected by GPTree and/or GPNode. See CreateWrappedGP.R for more information.
#

msgprefix <- paste(packageName(), ": ", sep = "")

#' R6 class WrappedmlegpGP
#'
#' @description
#' Contains the GP created by [mlegp::mlegp] from the \code{mlegp} package
#' @details
#' This package is by default not able to include individual uncertainties for input points. For this reason, all fields related to \code{y_var} are not used when updating the GP. No covariance kernel can be specified either. This implementation also assumes a vector for \code{y} (and not a matrix with multiple columns). Moreover, since no parameters can be specified for the GP, we will only update the GP parameters due to internal dependencies, but not use \code{init_covpars}.
#' 
WrappedmlegpGP <- R6::R6Class("WrappedmlegpGP",
  public = list(
    
    #' @field gp The mlegp GP object ([mlegp::mlegp] in the \code{mlegp} manual)
    gp = NULL,
    
    #' @field X_buffer Buffer matrix to collect x points until first GP can be trained
    X_buffer = NULL,
    
    #' @field y_buffer Buffer vector to collect y points until first GP can be trained
    y_buffer = NULL,
    
    #' @field y_var_buffer Buffer vector to collect variance of y points until first GP can be trained
    y_var_buffer = NULL,
    
    #' @field add_y_var Small additional variance used to keep the covariance matrix condition number under control
    add_y_var = 0,
    
    #' @field n_points_train_limit Number of points needed before we can create the GP
    n_points_train_limit = NULL,
    
    #' @field n_points The number of collected points belonging to this GP
    n_points = NULL,

    #' @field x_dim Dimensionality of input points
    x_dim = NULL,
    
    #' @field gp_control A list of GP implementation-specific options, passed directly to the wrapped GP implementation
    gp_control = NULL,
    
    #' @field init_covpars The initial covariance parameters when training the mlegp GP object in self@gp
    init_covpars = NULL,
    
    #' @field estimate_covpars If TRUE, the parameters are estimated by the package. Otherwise, the parameters from init_covpars are taken
    estimate_covpars = NULL,

    #' @field retrain_buffer_length Only retrain after this many new points have been added to the buffer
    retrain_buffer_length = NULL,
    
    #' @field retrain_buffer_counter Counter for the number of new points added since last retraining
    retrain_buffer_counter = NULL,
    
    #' @field add_buffer_in_prediction If TRUE, points in the data buffers are added to the GP before prediction. They are added into a temporarily created GP which contains the not yet included points. The GP in the node is not yet updated.
    add_buffer_in_prediction = NULL,
    
    #' @field X_shared Matrix with x points that this GP shares with the GP in the sibling node
    X_shared = NULL,
    
    #' @field y_shared Vector of y points that this GP shares with the GP in the sibling node
    y_shared = NULL,
    
    #' @field y_var_shared Vector of y_var points that this GP shares with the GP in the sibling node
    y_var_shared = NULL,
    
    #' @field n_shared_points The number of own points shared with the GP in the sibling node
    n_shared_points = NULL,
    
    
    #' @description
    #' Create a new WrappedmlegpGP object
    #' @param X Input data matrix with x_dim columns and at maximum Nbar rows. Is used to create the first iteration of the local GP.
    #' @param y Value of target variable at input point x; has to be a one-dimensional matrix or a vector; any further columns will be ignored
    #' @param y_var Variance of the target variable; has to be a one-dimensional matrix or vector
    #' @param gp_control A list of GP implementation-specific options, passed directly to the wrapped GP implementation
    #' @param init_covpars Initial covariance parameters of the local GP
    #' @param retrain_buffer_length Only retrain when the number of buffer points or collected points exceeds this value
    #' @param estimate_covpars If TRUE, the parameters are estimated by the package. Otherwise, the parameters from init_covpars are taken
    #' @param add_buffer_in_prediction If TRUE, points in the data buffers are added to the GP before prediction. They are added into a temporarily created GP which contains the not yet included points. The GP in the node is not yet updated.
    #' @param X_shared Matrix with x points that this GP shares with the GP in the sibling node
    #' @param y_shared Vector of y points that this GP shares with the GP in the sibling node
    #' @param y_var_shared Vector of y_var points that this GP shares with the GP in the sibling node
    #' @return A new WrappedmlegpGP object. Besides the local GP, information on the shared points and those stored in the buffer are collected. For more information on the GP, consult the method [mlegp::mlegp] in the \code{mlegp} package.
    #' @export
    initialize = function(X, y, y_var, gp_control, init_covpars, retrain_buffer_length, add_buffer_in_prediction, estimate_covpars=TRUE, X_shared=NULL, y_shared=NULL, y_var_shared=NULL) {
      
      self$x_dim <- ncol(X)
      self$init_covpars <- init_covpars
      self$retrain_buffer_length <- max(1, retrain_buffer_length)  # Cannot have retrain_buffer_length < 1
      self$estimate_covpars <- estimate_covpars
      self$add_buffer_in_prediction <- add_buffer_in_prediction
      
      self$n_points_train_limit <- max(3, self$x_dim + 1)
      
      # Set some default values for mlegp options to avoid that arguments set to NULL are not initialized correctly
      self$gp_control = list(
        constantMean = 1,
        nugget = NULL,
        nugget.known = 0, 
        min.nugget = 0,
        param.names = NULL,
        gp.names = NULL,
        PC.UD = NULL,
        PC.num = NULL,
        PC.percent = NULL, 
        simplex.ntries = 1,
        simplex.maxiter = 500,
        simplex.reltol = 1e-8,  
        BFGS.maxiter = 500,
        BFGS.tol = 0.01,
        BFGS.h = 1e-10,
        seed = 42, 
        verbose = 0,
        parallel = FALSE
      )
      
      # Remove all NULL elements from gp_control, to avoid element deletion 
      # when assigning values to self$gp_control below
      gp_control <- gp_control[!vapply(gp_control, is.null, logical(1))]
      if ("control" %in% names(gp_control)) {
        gp_control$control <- gp_control$control[!vapply(gp_control$gp_control, is.null, logical(1))]
      }
      
      # Assign remaining gp_control elements to the corresponding elements of self$gp_control
      for(name in names(gp_control)) {
        self$gp_control[[name]] <- gp_control[[name]]
      }
      for(name in names(gp_control$control)) {
        self$gp_control$control[[name]] <- gp_control$control[[name]]
      }
      
      # Initialize remaining variables
      self$X_buffer <- X
      self$y_buffer <- y
      self$y_var_buffer <- y_var
      self$n_points <- nrow(X)
      
      self$retrain_buffer_counter <- self$n_points
      
      self$X_shared <- X_shared
      self$y_shared <- y_shared
      self$y_var_shared <- y_var_shared
      self$n_shared_points <- 0
      if (!is.null(X_shared)) {
        self$n_shared_points <- nrow(X_shared)
      }
      
      # Try creating the GP right away
      self$train(do_buffer_check = FALSE)
    },
    
    
    #' @description
    #' Stores the initial covariance parameters (length-scales, standard deviation and trend coefficients) of the GP in the field \code{init_covpars}
    update_init_covpars = function() {
      self$init_covpars <- list(lengthscales=self$gp$beta,
                                sd2=self$gp$sig2,
                                # Assuming constant mean, we only return the first mean
                                trend.coef=self$gp$mu[1]
      )
    },
    
    #' @description
    #' Retrieves the length-scales of the kernel of the local GP
    get_lengthscales = function() {
      return(self$gp$beta)
    },
    
    
    #' @description
    #' Retrieves the design matrix X
    #' @param include_shared If TRUE, shared points between this GP and its sibling GP are included
    get_X_data = function(include_shared = FALSE) {
      X <- self$X_buffer
      if (include_shared) {
        X <- rbind(X, self$X_shared)
      }
      return(X)
    },
    
    
    #' @description
    #' Retrieves the response
    #' @param include_shared If TRUE, shared points between this GP and its sibling GP are included
    get_y_data = function(include_shared = FALSE) {
      y <- self$y_buffer
      if (include_shared) {
        y <- rbind(y, self$y_shared)
      }
      return(y)
    },
    
    
    #' @description
    #' Retrieves the individual variances from the response
    #' @param include_shared If TRUE, shared points between this GP and its sibling GP are included
    get_y_var_data = function(include_shared = FALSE) {
      y_var <- self$y_var_buffer
      if (include_shared) {
        y_var <- rbind(y_var, self$y_var_shared)
      }
      return(y_var)
    },
    
    
    #' @description
    #' Retrieves the covariance matrix
    #' @return the covariance matrix
    get_cov_mat = function() {
      if (is.null(self$gp)) {
        return(matrix())
      }
      ## computes the inverse of the inverse covariance matrix
      return(solve(self$gp$invVarMatrix, diag(nrow(self$gp$X))))
    },
    
    
    #' @description
    #' Method for updating add_y_var based on a bound for the covariance matrix condition number, based on \href{https://arxiv.org/abs/1602.00853}{this paper}, Section 5.4
    #' @param max_cond_num Max allowed condition number
    update_add_y_var = function(max_cond_num) {
      cond_num <- 0
      if (!is.null(self$gp)) {
        cov_mat <- self$get_cov_mat()
        
        eigvals <- eigen(cov_mat, symmetric=TRUE, only.values=TRUE)$values
        eigval_max <- eigvals[1]
        eigval_min <- eigvals[length(eigvals)]
        cond_num <- eigval_max / eigval_min
        message(msgprefix, "Covariance matrix condition number: ", cond_num)
        
        self$add_y_var <- max(0, (eigval_max - max_cond_num * eigval_min) / (max_cond_num - 1.0))
      }
      return(cond_num)
    },
    
    
    #' @description
    #' Stores a new point into the respective buffer method
    #' @param x Single input data point from the data stream; has to be a vector or row matrix with length equal to x_dim
    #' @param y Value of target variable at input point x; has to be a one-dimensional matrix or a vector; any further columns will be ignored
    #' @param y_var Variance of the target variable; has to be a one-dimensional matrix or vector
    #' @param shared If TRUE, this point is shared between this GP and its sibling GP
    #' @param remove_shared If TRUE, the last of the shared points is removed
    store_point = function(x, y, y_var, shared = FALSE, remove_shared = TRUE) {
      
      # Add the point to the correct buffers
      if (shared) {
        
        self$X_shared <- rbind(x, self$X_shared)  # Note the order of (x, self$X_shared) -- we're adding the new point at the beginning
        self$y_shared <- rbind(c(y), self$y_shared)
        self$y_var_shared <- rbind(c(y_var), self$y_var_shared)
        self$n_shared_points <- self$n_shared_points + 1
        
      } else {
        
        self$X_buffer <- rbind(self$X_buffer, x)
        self$y_buffer <- rbind(self$y_buffer, y)
        self$y_var_buffer <- rbind(self$y_var_buffer, y_var)
        self$n_points <- self$n_points + 1
        self$retrain_buffer_counter <- self$retrain_buffer_counter + 1
        
      }
      
      if (remove_shared) {
        if (self$n_shared_points > 0) {
          n_keep <- self$n_shared_points - 1
          self$X_shared <- head(self$X_shared, n_keep)  # Remove the last of the shared points
          self$y_shared <- head(self$y_shared, n_keep)
          self$y_var_shared <- head(self$y_var_shared, n_keep)
          self$n_shared_points <- n_keep
        }
      }
    },
    
    
    #' @description
    #' Method for clearing the buffers
    delete_buffers = function() {
      # Get rid of the buffered data if the GP was created
      self$X_buffer <- matrix(ncol = self$x_dim, nrow = 0, dimnames = list(c(), colnames(self$X_buffer)))
      self$y_buffer <- numeric(0)
      self$y_var_buffer <- numeric(0)
    },
    
    
    #' @description
    #' Method for (re)creating / (re)training the GP
    #' @param do_buffer_check If TRUE, only train the GP if the number of stored points is larger than retrain_buffer_length
    #' @return TRUE if training was performed, otherwise FALSE
    train = function(do_buffer_check = TRUE) {
      
      ## Do not train if the n_points < n_points_train_limit
      n_total <- self$n_points + self$n_shared_points
      if (n_total < self$n_points_train_limit) {
        message(msgprefix, "Will not train the GP since the number of collected points (", n_total, ") is below the minimum needed (", self$n_points_train_limit, ").")
        return(FALSE)
      }
      
      ## Do not train if the buffers are empty or shorter than self$retrain_buffer_length
      if (do_buffer_check) {
        if (self$retrain_buffer_counter < self$retrain_buffer_length) {
          message(msgprefix, "Will not train the GP since retrain_buffer_counter (", self$retrain_buffer_counter, ") or number of collected points (", self$n_points, ") are below retrain_buffer_length (", self$retrain_buffer_length, ").")
          return(FALSE)
        }
      }
      
      # Construct the total set of input points
      X <- NULL
      y <- NULL
      y_var <- NULL
      
      X <- rbind(X, self$X_buffer)
      y <- rbind(y, self$y_buffer)
      y_var <- rbind(y_var, self$y_var_buffer)
      
      # Turn X_shared into matrix if X is one-dimensional
      if(!is.null(self$X_shared)){
        self$X_shared <- as.matrix(self$X_shared, ncol = self$x_dim)
      }
      
      X <- rbind(X, self$X_shared)
      y <- rbind(y, self$y_shared)
      y_var <- rbind(y_var, self$y_var_shared)
      
      # Recreate the self$gp object
      self$delete_gp()
      self$create_mlegp_gp(X, y, y_var + self$add_y_var)
      
      # Reset the retrain_buffer_counter
      self$retrain_buffer_counter <- 0
      
      return(TRUE)
    },
    
    
    #' @description
    #' Method for prediction
    #' @param x Single data point for which the predicted mean (and standard deviation) is computed; has to be a vector or row matrix with length equal to x_dim
    #' @param return_std If TRUE, the standard error is returned in addition to the prediction
    #' @return Prediction for input point x
    predict = function(x, return_std=TRUE) {
      
      predict_result <- NULL
      
      if (is.null(self$gp)) {
        
        message(msgprefix, "The GP has not been created yet. Returning 0 +- Inf as prediction.")
        y_pred_with_err <- c(0., Inf)
        if (return_std) {
          return(y_pred_with_err)
        } else {
          return(y_pred_with_err[1])
        }
        
      } else if (self$add_buffer_in_prediction & self$retrain_buffer_counter > 0) {
        
        # Construct the total set of input points
        X <- NULL
        y <- NULL
        y_var <- NULL
        
        X <- rbind(X, self$X_buffer)
        y <- rbind(y, self$y_buffer)
        y_var <- rbind(y_var, self$y_var_buffer)
        
        X <- rbind(X, self$X_shared)
        y <- rbind(y, self$y_shared)
        y_var <- rbind(y_var, self$y_var_shared)
        
        temp_gp <- mlegp::mlegp(
          X, y,
          constantMean = self$gp_control$constantMean,
          nugget = self$gp_control$nugget,
          nugget.known = self$gp_control$nugget.known
        )
        
        # Get prediction using our temp_gp object
        predict_result <- self$call_mlegp_predict(matrix(x, ncol = self$x_dim), use_gp = temp_gp)
        
      } else {
        
        # Get prediction using self$gp
        predict_result <- self$call_mlegp_predict(matrix(x, ncol = self$x_dim))
        
      }
      
      # Put prediction and error in a two-element vector
      y_pred_with_err <- c(predict_result$mean, predict_result$sd)
      
      if (return_std) {
        return(y_pred_with_err)
      } else {
        return(y_pred_with_err[1])
      }
      
    },
    
    
    #' @description
    #' Method to delete the GP object in self$gp
    delete_gp = function() {
      self$gp <- NULL
    },
    
    
    #' @description
    #' Method for calling the 'mlegp' function in mlegp to create a GP object, stored in self$gp
    #' @param X Input data matrix with x_dim columns and at maximum Nbar rows for the local GP.
    #' @param y Value of target variable at input point x; has to be a one-dimensional matrix or a vector; any further columns will be ignored
    #' @param y_var Variance of the target variable; has to be a one-dimensional matrix or vector
    #' @return TRUE
    create_mlegp_gp = function(X, y, y_var) {
      
      message(msgprefix, "Creating/training the GP using ", length(y), " data points.")
      # Create the GP
      self$gp <- mlegp::mlegp(
        X, y,
        constantMean = self$gp_control$constantMean,
        nugget = self$gp_control$nugget,
        nugget.known = self$gp_control$nugget.known,
        min.nugget = self$gp_control$min.nugget,
        param.names = self$gp_control$param.names,
        gp.names = self$gp_control$gp.names,
        PC.UD = self$gp_control$PC.UD,
        PC.num = self$gp_control$PC.num,
        PC.percent = self$gp_control$PC.percent,
        simplex.ntries = self$gp_control$simplex.ntries,
        simplex.maxiter = self$gp_control$simplex.maxiter,
        simplex.reltol = self$gp_control$simplex.reltol,
        BFGS.maxiter = self$gp_control$BFGS.maxiter,
        BFGS.tol = self$gp_control$BFGS.tol,
        BFGS.h = self$gp_control$BFGS.h,
        seed = self$gp_control$seed,
        verbose = self$gp_control$verbose,
        parallel = self$gp_control$parallel
      )
      
      # Register the new covariance parameters
      self$update_init_covpars()
      message(msgprefix, "New covariance parameters: length scales = ", paste(self$init_covpars$lengthscales, collapse=", "), "  sd2 = ", self$init_covpars$sd2, "  trend.coef = ", self$init_covpars$trend.coef)
      
      return(TRUE)
    },
    
    
    #' @description
    #' Method for calling the 'predict' function in mlegp
    #' @param x Single data point for which the predicted mean (and standard deviation) is computed; has to be a vector with length equal to x_dim
    #' @param use_gp Optional user-defined GP which is evaluated instead of the local GP
    #' @return The predictions for x from the specified GP, by default the local GP. The output needs to be a list with fields mean and sd for the prediction and prediction error, respectively.
    call_mlegp_predict = function(x, use_gp = NULL) {
      if (!is.null(use_gp)) {
        prediction <- mlegp::predict.gp(
          use_gp, newData = x,
          se.fit = TRUE
        )
        # Renaming the list so that the prediction is stored under $mean and the prediction error under $sd
        prediction <- list(mean = prediction$fit, sd = prediction$se.fit)
        return(prediction)
      } else {
        prediction <- mlegp::predict.gp(
          self$gp, newData = x,
          se.fit = TRUE
        )
        # Renaming the list so that the prediction is stored under $mean and the prediction error under $sd
        prediction <- list(mean = prediction$fit, sd = prediction$se.fit)
        return(prediction)
      }
    }
    
  )
)
