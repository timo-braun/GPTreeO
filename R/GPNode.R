# The GPNode class

msgprefix <- paste(packageName(), ": ", sep = "")

#' R6 Class for the nodes / leaves in the GPTree tree
#'
#' @description
#' The nodes contain the local GP if they are leaves (at the end of a branch). Nodes that are just nodes contain information on how the input space was split. They are responsible for computing and updating the splitting probabilities. Also, the tree interacts with the local GPs through the nodes.
#' 
#' Currently, GPs from the \code{DiceKriging} package (\link{WrappedDiceKrigingGP}) and \code{mlegp} package (\link{WrappedmlegpGP}) are implemented. The user can create their own wrapper using \link{WrappedGP}.
#' 
#' @seealso [GPTree()] for the main methods
#' 
GPNode <- R6::R6Class("GPNode",
  public = list(

    # Start user inputs
    #' @field key A string like "0110100" to identify the node in the binary tree
    key = NULL,
    
    #' @field x_dim Dimensionality of input points. It is set once the first point is received through the \link{GPTree} method \code{update}. It needs to be specified if \code{min_ranges} should be different from default.
    x_dim = NULL,
    
    #' @field theta Overlap ratio between two leafs in the split direction. The default value is 0.
    theta = NULL,
    
    #' @field split_direction_criterion A string that indicates which spitting criterion to use. The options are:
    #' * \code{"max_spread"}: Split along the direction which has the largest data spread.
    #' * \code{"min_lengthscale"}: split along the direction with the smallest length-scale hyperparameter from the local GP.
    #' * \code{"max_spread_per_lengthscale"}: Split along the direction with the largest data spread relative to the corresponding GP length-scale hyperparameter.
    #' * \code{"max_corr"}: Split along the direction where the input data is most strongly correlated with the target variable.
    #' * \code{"principal_component"}: Split along the first principal component.
    #' 
    #' The default value is \code{"max_spread_per_lengthscale"}.
    split_direction_criterion = NULL,
    
    #' @field split_position_criterion A string indicating how the split position along the split direction should be set. Possible values are (\code{"mean"} and \code{"median"}). The default is \code{"mean"}.
    split_position_criterion = NULL,
    
    #' @field shape_decay A string specifying how the probability function for a point to be assigned to the left leaf should fall off in the overlap region. The available options are a linear shape (\code{"linear"}), an exponential shape (\code{"exponential"}) or a Gaussian shape (\code{"gaussian"}). Another option is to select no overlap region. This can be achieved by selecting \code{"deterministic"} or to set \code{theta} to 0. The default is \code{"linear"}.
    shape_decay = NULL,
    
    #' @field prob_min_theta Minimum probability after which the overlap shape gets truncated (either towards 0 or 1). The default value is 0.01.
    prob_min_theta = NULL,
    
    #' @field Nbar Maximum number of data points for each GP in a leaf before it is split. The default value is 1000.
    Nbar = NULL,

    #' @field min_ranges Smallest allowed input data spread (per dimension) before node splitting stops. It is set to its default \code{min_ranges = rep(0.0, x_dim)} once the first point is received through the \code{update} method. \code{x_dim} needs to be specified by the user if it should be different from the default.
    min_ranges = NULL,
    
    #' @field is_leaf If TRUE, this node a leaf, i.e the last node on its branch
    is_leaf = NULL,
    # End user inputs
    
    # Start auxiliary fields
    #' @field wrapped_gp An instance of the WrappedGP type
    wrapped_gp = NULL,
    
    #' @field can_split If TRUE for a given dimension, the leaf can be split along that dimension
    can_split = NULL,

    #' @field rotation_matrix A rotation matrix, used for transforming the data
    rotation_matrix = NULL,

    #' @field shift A shift, used for transforming the data
    shift = NULL,

    #' @field use_pc_transform TRUE if principal components transformation is used for node splitting
    use_pc_transform = NULL,

    #' @field x_spread Vector of data spread for each dimension
    x_spread = NULL,

    #' @field split_index Index for the split dimension
    split_index = NULL,

    #' @field position_split Position of the split along dimension split_index
    position_split = NULL,

    #' @field width_overlap Width of overlap region along dimension split_index
    width_overlap = NULL,

    #' @field point_ids IDs of the points assigned to this node
    point_ids = NULL,

    #' @field residuals Vector of residuals
    residuals = NULL,

    #' @field pred_errs Vector of prediction uncertainties
    pred_errs = NULL,

    #' @field error_scaler Scaling factor for the prediction error to ensure desired coverage
    error_scaler = NULL,

    #' @field use_n_residuals Number of past residuals to use in calibrating the \code{error_scaler}
    use_n_residuals = NULL,
    # End auxiliary fields
    

    #' @description
    #' Create a new node object
    #' @param key A string like "0110100" to identify the node in the binary tree
    #' @param x_dim Dimensionality of input points. It is set once the first point is received through the \link{GPTree} method \code{update}. It needs to be specified if \code{min_ranges} should be different from default.
    #' @param theta Overlap ratio between two leafs in the split direction. The default value is 0.
    #' @param split_direction_criterion A string that indicates which spitting criterion to use. The options are:
    #' * \code{"max_spread"}: Split along the direction which has the largest data spread.
    #' * \code{"min_lengthscale"}: split along the direction with the smallest length-scale hyperparameter from the local GP.
    #' * \code{"max_spread_per_lengthscale"}: Split along the direction with the largest data spread relative to the corresponding GP length-scale hyperparameter.
    #' * \code{"max_corr"}: Split along the direction where the input data is most strongly correlated with the target variable.
    #' * \code{"principal_component"}: Split along the first principal component.
    #' 
    #' The default value is \code{"max_spread_per_lengthscale"}.
    #' @param split_position_criterion A string indicating how the split position along the split direction should be set. Possible values are (\code{"mean"} and \code{"median"}). The default is \code{"mean"}.
    #' @param shape_decay A string specifying how the probability function for a point to be assigned to the left leaf should fall off in the overlap region. The available options are a linear shape (\code{"linear"}), an exponential shape (\code{"exponential"}) or a Gaussian shape (\code{"gaussian"}). Another option is to select no overlap region. This can be achieved by selecting \code{"deterministic"} or to set \code{theta} to 0. The default is \code{"linear"}.
    #' @param prob_min_theta Minimum probability after which the overlap shape gets truncated (either towards 0 or 1). The default value is 0.01.
    #' @param Nbar Maximum number of data points for each GP in a leaf before it is split. The default value is 1000.
    #' @param wrapper A string that indicates which GP implementation should be used. The current version includes wrappers for the packages \code{"DiceKriging"} and \code{"mlegp"}. The default setting is \code{"DiceKriging"}.
    #' @param gp_control A \code{list} of control parameter that is forwarded to the wrapper. Here, the covariance function is specified. \code{DiceKriging} allows for the following kernels, passed as string: \code{"gauss"}, \code{"matern5_2"}, \code{"matern3_2"}, \code{"exp"}, \code{"powexp"} where \code{"matern3_2"} is set as default.
    #' @param retrain_buffer_length Size of the retrain buffer. The buffer for a each node collects data points and holds them until the buffer length is reached. Then the GP in the node is updated with the data in the buffer. For a fixed \code{Nbar}, higher values for \code{retrain_buffer_length} lead to faster run time (less frequent retraining), but the trade-off is a temporary reduced prediction accuracy. We advise that the choice for \code{retrain_buffer_length} should depend on the chosen \code{Nbar}. By default \code{retrain_buffer_length} is set equal to \code{Nbar}.
    #' @param n_points_train_limit Number of points at which a GP is created in the leaf
    #' @param add_buffer_in_prediction If TRUE, points in the data buffers are added to the GP before prediction. They are added into a temporarily created GP which contains the not yet included points. The GP in the node is not yet updated. The default is \code{FALSE}.
    #' @param min_ranges Smallest allowed input data spread (per dimension) before node splitting stops. It is set to its default \code{min_ranges = rep(0.0, x_dim)} once the first point is received through the \link{GPTree} method \code{update}. \code{x_dim} needs to be specified by the user if it should be different from the default.
    #' @param is_leaf If TRUE, this node a leaf, i.e the last node on its branch.
    #' @return A new GPNode object. Contains the local GP in the field \code{wrapped_gp}, and information used for and related to splitting the node. If the node has been split, the local GP is removed.
    #' @md
    #' @export
    initialize = function(key, x_dim, theta, split_direction_criterion, split_position_criterion, shape_decay, prob_min_theta, Nbar, wrapper, gp_control, retrain_buffer_length, add_buffer_in_prediction, min_ranges = NULL, is_leaf = TRUE) {

      # Empty data frames to initialize GP object
      X <- matrix(ncol = x_dim, nrow = 0)
      colnames(X) <- as.character(paste("x", 1:x_dim, sep = "")) # names "x1", "x2", ...

      y <- numeric(0)
      y_var <- numeric(0)

      # Initialize the gp field with a new GP instance of the requested type
      self$wrapped_gp <- CreateWrappedGP(wrapper, X, y, y_var, gp_control, NULL, retrain_buffer_length, add_buffer_in_prediction)


      # Initialize other fields
      self$key <- key
      self$Nbar <- Nbar
      self$split_direction_criterion <- split_direction_criterion
      self$shape_decay <- shape_decay
      self$theta <- theta
      self$x_dim <- ncol(self$wrapped_gp$get_X_data())
      self$is_leaf <- is_leaf

      self$min_ranges <- min_ranges
      self$can_split <- rep(TRUE, x_dim)

      self$x_spread <- rep(NA, x_dim)
      self$split_index <- NA
      self$split_position_criterion <- split_position_criterion
      self$position_split <- NA
      self$width_overlap <- NA
      self$prob_min_theta <- prob_min_theta
      self$point_ids <- c()

      self$rotation_matrix <- diag(1, nrow = x_dim)
      self$shift <- rep(0, x_dim)
      self$use_pc_transform <- FALSE
      if (self$split_direction_criterion == "principal_component") {
        self$use_pc_transform <- TRUE
      }

      self$error_scaler <- 10
      self$use_n_residuals <- 25

      self$residuals <- rep(0, self$use_n_residuals)
      self$pred_errs <- rep(Inf, self$use_n_residuals)

      # Check settings consistency
      if (!(self$shape_decay %in% c("linear", "exponential", "gaussian", "deterministic"))) {
        stop("The option shape_decay = '", self$shape_decay, "' is not implemented. Please use one of these: linear, exponential, gaussian, deterministic.")
      }

      if (!(self$split_direction_criterion %in% c("max_spread", "min_lengthscale", "max_spread_per_lengthscale", "max_corr", "principal_component"))) {
        stop("The option split_direction_criterion = '", self$split_direction_criterion, "' is not implemented. Please use one of these: max_spread, min_lengthscale, max_spread_per_lengthscale, max_corr.")
      }

      if (!(self$split_position_criterion %in% c("mean", "median"))) {
        stop("The option split_direction_criterion = '", self$split_direction_criterion, "' is not implemented. Please use one of these: mean, median")
      }


    },


    #' @description
    #' Method to transform input data through a shift and a rotation. IS EXPECTED TO NOT BE CALLED BY THE USER
    #' @param X Matrix with x points
    #' @return The transformed X matrix
    transform = function(X) {

      # Shift 
      X <- X - matrix(self$shift, nrow=nrow(X), ncol=self$x_dim, byrow=TRUE)

      # Rotate
      X <- X %*% self$rotation_matrix

      return(X)
    },


    #' @description
    #' Method to update the probability parameters (x_spread, can_split, split_index, position_split, width_overlap). IS EXPECTED TO NOT BE CALLED BY THE USER
    update_prob_pars = function() {

      # Get the a local copy of the X data for this node
      X_node <- self$wrapped_gp$get_X_data()

      # If applicable, update self$rotation_matrix and self$shift 
      # and transform X_node data
      if (self$use_pc_transform) {

        # Perform PCA
        pc_output <- prcomp(X_node, center = TRUE, scale. = FALSE)

        # Update the parameters used by self$transform
        self$rotation_matrix <- pc_output$rotation
        self$shift <- pc_output$center

        # Transform local data copy
        X_node <- self$transform(X_node)
      }

      # Update 
      # - self$x_spread with the data spread in each direction
      # - self$can_split with TRUE/FALSE for each direction
      for (index_dim in 1:self$x_dim) {
        spread <- max(X_node[, index_dim]) - min(X_node[, index_dim])
        self$x_spread[[index_dim]] <- spread

        if (spread < self$min_ranges[[index_dim]]) {
          self$can_split[[index_dim]] <- FALSE
        }
      }

      # If we can't split the node further, we can't update the probability parameters
      if (!any(self$can_split)) {
        return(FALSE)
      }

      # Update self$split_index with the index for the direction we are 
      # splitting along.
      if (self$split_direction_criterion == "max_spread") {
        self$split_index <- NA
        largest_spread <- 0
        for (index_dim in 1:self$x_dim) {
          if (self$can_split[[index_dim]] & self$x_spread[[index_dim]] > largest_spread) {
            self$split_index <- index_dim
            largest_spread <- self$x_spread[[index_dim]]
          }
        }
      } else if (self$split_direction_criterion == "min_lengthscale") {
        self$split_index <- NA
        lengthscales <- self$wrapped_gp$get_lengthscales()
        smallest_lengthscale <- .Machine$double.xmax
        for (index_dim in 1:self$x_dim) {
          if (self$can_split[[index_dim]] & lengthscales[[index_dim]] < smallest_lengthscale) {
            self$split_index <- index_dim
            smallest_lengthscale <- lengthscales[[index_dim]]
          }
        }
      } else if (self$split_direction_criterion == "max_spread_per_lengthscale") {
        self$split_index <- NA
        lengthscales <- self$wrapped_gp$get_lengthscales()
        spread_per_lengthscale <- self$x_spread / lengthscales
        largest_spread <- 0
        for (index_dim in 1:self$x_dim) {
          if (self$can_split[[index_dim]] & spread_per_lengthscale[[index_dim]] > largest_spread) {
            self$split_index <- index_dim
            largest_spread <- spread_per_lengthscale[[index_dim]]
          }
        }
      } else if (self$split_direction_criterion == "max_corr") {
        y_data_node <- self$wrapped_gp$get_y_data()
        self$split_index <- NA
        largest_cor <- 0
        for (index_dim in 1:self$x_dim) {
          abs_cor <- abs(cor(X_node[, index_dim], y_data_node))
          if (self$can_split[[index_dim]] & abs_cor > largest_cor) {
            self$split_index <- index_dim
            largest_cor <- abs_cor
          }
        }
      } else if (self$split_direction_criterion == "principal_component") {
        # Set split_index to the first of the principal components for which can_split is TRUE
        self$split_index <- NA
        for (index_dim in 1:self$x_dim) {
          if (self$can_split[[index_dim]]) {
            self$split_index <- index_dim
            break
          }
        }
      }
        
      # Update self$width_overlap
      self$width_overlap <- self$theta * self$x_spread[[self$split_index]]

      # Update self$position_split
      # This method should only be called when the node is full,
      # i.e. when the number of data points equals the cap size (Nbar).
      if (self$split_position_criterion == "mean") {
        self$position_split <- mean(X_node[, self$split_index])
      } else if (self$split_position_criterion == "median") {
        self$position_split <- median(X_node[, self$split_index])
      }

      return(TRUE)
    },


    #' @description
    #' Method to compute the probability that a point x should go to child 1. IS EXPECTED TO NOT BE CALLED BY THE USER
    #' @param x Single data point for which probability is computed; has to be a vector with length equal to x_dim
    #' @return The probability that a point x should go to child 1
    get_prob_child_1 = function(x) {

      # Apply prior transformation and get the relevant component of the input point
      if (self$use_pc_transform) {
        x <- as.vector(self$transform(matrix(x, ncol=self$x_dim)))
      }

      x_split_direction <- x[[self$split_index]]
      x_normalized <- (x_split_direction - self$position_split) / self$x_spread[self$split_index]

      if (self$theta <= 0 | self$shape_decay == "deterministic") {

        if (x_split_direction <= self$position_split) { probability <- 0.0 }
        else if (self$position_split < x_split_direction) { probability <- 1.0 }

      } else if (self$shape_decay == "linear") {

        probability <- x_normalized / self$theta + 0.5
        if (probability < 0) { probability <- 0 }
        if (probability > 1) { probability <- 1 }

      } else if (self$shape_decay == "exponential") {

        x_temp <- x_normalized / as.double((self$theta / 2) / (-log(self$prob_min_theta * 2)))
        probability <- 1 - (0.5 + 0.5 * sign(x_temp) * (exp(-abs(x_temp)) - 1))
        if (probability < self$prob_min_theta) { probability <- 0 }
        if (probability > 1 - self$prob_min_theta) { probability <- 1 }

      } else if (self$shape_decay == "gaussian") {

        x_temp <- (x_normalized / as.double((self$theta / 2) / (sqrt(-2 * log(self$prob_min_theta * 2)))))
        probability <- 1 - (0.5 + 0.5 * sign(x_temp) * (exp(-0.5 * x_temp**2) - 1))
        if (probability < self$prob_min_theta) { probability <- 0 }
        if (probability > 1 - self$prob_min_theta) { probability <- 1 }

      } else {
        stop(msgprefix, "Unknown value for the option shape_decay. This is a bug as it should have been caught already in GPNode$initialize.")
      }

      return(probability)
    },


    #' @description
    #' Method to register prediction performance 
    #' @param x Most recent single input data point from the data stream; has to be a vector with length equal to x_dim
    #' @param y Target variable which has to be a one-dimensional matrix or a vector; any further columns will be ignored
    register_residual = function(x, y) {

      # Get prediction w/ uncertainty
      pred_with_err <- self$wrapped_gp$predict(x, return_std=TRUE)
      y_pred <- pred_with_err[1]
      y_pred_std <- pred_with_err[2]

      self$residuals <- c(rep(NA, times=1), self$residuals[1:(self$use_n_residuals-1)])
      self$residuals[1] <- y - y_pred

      self$pred_errs <- c(rep(NA, times=1), self$pred_errs[1:(self$use_n_residuals-1)])
      self$pred_errs[1] <- y_pred_std
    },


    #' @description
    #' Method for updating the empirical error parameters
    update_empirical_error_pars = function() {

      target_coverage <- 0.68
      res <- uniroot(
        function(x) sum(abs(self$residuals) < x * self$pred_errs) / self$use_n_residuals - target_coverage,
        c(0.2 * self$error_scaler, 5 * self$error_scaler),
        extendInt = "yes",
        maxiter = 100
      )

      self$error_scaler <- res$root
    },


    #' @description
    #' Method to delete the GP. IS EXPECTED TO NOT BE CALLED BY THE USER
    delete_gp = function() {
      self$wrapped_gp$delete_gp()
      self$wrapped_gp <- NULL
    }
  )
)
