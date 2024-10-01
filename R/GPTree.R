# The GPTree class

msgprefix <- paste(packageName(), ": ", sep = "")

#' Tree structure storing all nodes containing local GPs
#'
#' @description
#' The base class which contains and where all parameters are set. Here, all information on how and when the splitting is carried out is stored.
#' \code{wrapper} and \code{gp_control} specify the Gaussian process (GP) implementation and its parameters. Moreover, minimum errors and calibration of the predictions are specified here, too. 
#' 
#' **Essential methods**
#' 
#' 
#' The following three methods are essential for the package. The remaining ones are mostly not expected to be called by the user.
#' * \href{#method-GPTree-new}{\code{GPTree$new()}}: Creates a new tree with specified parameters
#' * \href{#method-GPTree-update}{\code{GPTree$update()}}: Adds the information from the input point to the tree and updates local GPs
#' * \href{#method-GPTree-joint_prediction}{\code{GPTree$joint_prediction()}}: Computes the joint prediction for a given input point
#' 
#' @section Brief package functionality overview:
#' The tree collects the information from all \link{GPNode}s which in turn contain the local GP. Currently, GPs from the \code{DiceKriging} package (\link{WrappedDiceKrigingGP}) and \code{mlegp} package (\link{WrappedmlegpGP}) are implemented. The user can create their own wrapper using \link{WrappedGP}.
#' 
#' 
#' @import R6
#' @import hash
#' @import DiceKriging
#' @import mlegp
#' 
GPTree <- R6::R6Class("GPTree",
  public = list(

    # Start user inputs
    #' @field Nbar Maximum number of data points for each GP in a leaf before it is split. The default value is 1000.
    Nbar = NULL,

    #' @field retrain_buffer_length Size of the retrain buffer. The buffer for a each node collects data points and holds them until the buffer length is reached. Then the GP in the node is updated with the data in the buffer. For a fixed \code{Nbar}, higher values for \code{retrain_buffer_length} lead to faster run time (less frequent retraining), but the trade-off is a temporary reduced prediction accuracy. We advise that the choice for \code{retrain_buffer_length} should depend on the chosen \code{Nbar}. By default \code{retrain_buffer_length} is set equal to \code{Nbar}.
    retrain_buffer_length = NULL,
    
    #' @field gradual_split If TRUE, gradual splitting is used for splitting. The default value is TRUE.
    gradual_split = NULL,
    
    #' @field theta Overlap ratio between two leafs in the split direction. The default value is 0.
    theta = NULL,
    
    #' @field wrapper A string that indicates which GP implementation should be used. The current version includes wrappers for the packages \code{"DiceKriging"} and \code{"mlegp"}. The default setting is \code{"DiceKriging"}.
    wrapper = NULL,
    
    #' @field gp_control A \code{list} of control parameter that is forwarded to the wrapper. Here, the covariance function is specified. \code{DiceKriging} allows for the following kernels, passed as string: \code{"gauss"}, \code{"matern5_2"}, \code{"matern3_2"}, \code{"exp"}, \code{"powexp"} where \code{"matern3_2"} is set as default.
    gp_control = NULL,
    
    #' @field split_direction_criterion A string that indicates which spitting criterion to use. The options are:
    #' * \code{"max_spread"}: Split along the direction which has the largest data spread.
    #' * \code{"min_lengthscale"}: split along the direction with the smallest length-scale hyperparameter from the local GP.
    #' * \code{"max_spread_per_lengthscale"}: Split along the direction with the largest data spread relative to the corresponding GP length-scale hyperparameter.
    #' * \code{"max_corr"}: Split along the direction where the input data is most strongly correlated with the target variable.
    #' * \code{"principal_component"}: Split along the first principal component.
    #' 
    #' The default value is \code{"max_spread_per_lengthscale"}.
    split_direction_criterion = NULL,

    #' @field split_position_criterion A string indicating how the split position along the split direction should be set. Possible values are (\code{"median"} and \code{"mean"}). The default is \code{"median"}.
    split_position_criterion = NULL,

    #' @field shape_decay A string specifying how the probability function for a point to be assigned to the left leaf should fall off in the overlap region. The available options are a linear shape (\code{"linear"}), an exponential shape (\code{"exponential"}) or a Gaussian shape (\code{"gaussian"}). Another option is to select no overlap region. This can be achieved by selecting \code{"deterministic"} or to set \code{theta} to 0. The default is \code{"linear"}.
    shape_decay = NULL,

    #' @field use_empirical_error If TRUE, the uncertainty is calibrated using recent data points. The default value is TRUE.
    #'
    #' The most recent 25 observations are used to ensure that the prediction uncertainty yields approximately 68 % coverage. This coverage is only achieved if \code{theta = 0} (also together with \code{gradual_split = TRUE}) is used. Nevertheless, the coverage will be closer to 68 % than it would be without calibration. The prediction uncertainties at the beginning are conservative and become less conservative with increasing number of input points.
    use_empirical_error = NULL,

    #' @field use_reference_gp If TRUE, the covariance parameters determined for the GP in node 0 will be used for all subsequent GPs. The default is \code{FALSE}.
    use_reference_gp = NULL,
     
    #' @field min_abs_y_err Minimum absolute error assumed for y data. The default value is 0.
    min_abs_y_err = NULL,

    #' @field min_rel_y_err Minimum relative error assumed for y data. The default value is \code{100 * .Machine$double.eps}.
    min_rel_y_err = NULL,

    #' @field min_abs_node_pred_err Minimum absolute error on the prediction from a single node. The default value is 0.
    min_abs_node_pred_err = NULL,

    #' @field min_rel_node_pred_err Minimum relative error on the prediction from a single node. The default value is \code{100 * .Machine$double.eps}.
    min_rel_node_pred_err = NULL,

    #' @field prob_min_theta Minimum probability after which the overlap shape gets truncated (either towards 0 or 1). The default value is 0.01.
    prob_min_theta = NULL,
    
    #' @field add_buffer_in_prediction If TRUE, points in the data buffers are added to the GP before prediction. They are added into a temporarily created GP which contains the not yet included points. The GP in the node is not yet updated. The default is \code{FALSE}.
    add_buffer_in_prediction = NULL,
    
    #' @field x_dim Dimensionality of input points. It is set once the first point is received through the \href{#method-GPTree-update}{\code{update()}} or \href{#method-GPTree-joint_prediction}{\code{joint_prediction()}} method. It needs to be specified if \code{min_ranges} should be different from default.
    x_dim = NULL,
    
    #' @field min_ranges Smallest allowed input data spread (per dimension) before node splitting stops. It is set to its default \code{min_ranges = rep(0.0, x_dim)} once the first point is received through the \href{#method-GPTree-update}{\code{update()}} method. \code{x_dim} needs to be specified by the user if it should be different from the default.
    min_ranges = NULL,
    
    #' @field max_cond_num Add additional noise if the covariance matrix condition number exceeds this value. The default is \code{NULL}.
    max_cond_num = NULL,

    #' @field max_points The maximum number of points the tree is allowed to store. The default value is \code{Inf}.
    #' 
    #' End of the user-defined input fields.
    max_points = NULL,
    # End user inputs
    
    # Start auxiliary fields
    #' @field nodes A hash to hold the GP tree, using string keys to identify nodes and their position in the tree  ("0", "00", "01", "000", "001", "010", "011", etc.)
    nodes = NULL,

    #' @field leaf_keys Stores the keys ("0", "00", "01", "000", "001", "010", "011", etc.) for the leaves
    leaf_keys = NULL,

    #' @field n_points Number of points in the tree
    n_points = NULL,

    #' @field n_fed Number of points fed to the tree
    n_fed = NULL,
    # End auxiliary fields
    
    #' @param Nbar Maximum number of data points for each GP in a leaf before it is split. The default value is 1000.
    #' @param retrain_buffer_length Size of the retrain buffer. The buffer for a each node collects data points and holds them until the buffer length is reached. Then the GP in the node is updated with the data in the buffer. For a fixed \code{Nbar}, higher values for \code{retrain_buffer_length} lead to faster run time (less frequent retraining), but the trade-off is a temporary reduced prediction accuracy. We advise that the choice for \code{retrain_buffer_length} should depend on the chosen \code{Nbar}. By default \code{retrain_buffer_length} is set equal to \code{Nbar}.
    #' @param gradual_split If TRUE, gradual splitting is used for splitting. The default value is TRUE.
    #' @param theta Overlap ratio between two leafs in the split direction. The default value is 0.
    #' @param wrapper A string that indicates which GP implementation should be used. The current version includes wrappers for the packages \code{"DiceKriging"} and \code{"mlegp"}. The default setting is \code{"DiceKriging"}.
    #' @param gp_control A \code{list} of control parameter that is forwarded to the wrapper. Here, the covariance function is specified. \code{DiceKriging} allows for the following kernels, passed as string: \code{"gauss"}, \code{"matern5_2"}, \code{"matern3_2"}, \code{"exp"}, \code{"powexp"} where \code{"matern3_2"} is set as default.
    #' @param split_direction_criterion A string that indicates which spitting criterion to use. The options are:
    #' * \code{"max_spread"}: Split along the direction which has the largest data spread.
    #' * \code{"min_lengthscale"}: split along the direction with the smallest length-scale hyperparameter from the local GP.
    #' * \code{"max_spread_per_lengthscale"}: Split along the direction with the largest data spread relative to the corresponding GP length-scale hyperparameter.
    #' * \code{"max_corr"}: Split along the direction where the input data is most strongly correlated with the target variable.
    #' * \code{"principal_component"}: Split along the first principal component.
    #' 
    #' The default value is \code{"max_spread_per_lengthscale"}.
    #' @param split_position_criterion A string indicating how the split position along the split direction should be set. Possible values are (\code{"median"} and \code{"mean"}). The default is \code{"median"}.
    #' @param shape_decay A string specifying how the probability function for a point to be assigned to the left leaf should fall off in the overlap region. The available options are a linear shape (\code{"linear"}), an exponential shape (\code{"exponential"}) or a Gaussian shape (\code{"gaussian"}). Another option is to select no overlap region. This can be achieved by selecting \code{"deterministic"} or to set \code{theta} to 0. The default is \code{"linear"}.
    #' @param use_empirical_error If TRUE, the uncertainty is calibrated using recent data points. The default value is TRUE.
    #'
    #' The most recent 25 observations are used to ensure that the prediction uncertainty yields approximately 68 % coverage. This coverage is only achieved if \code{theta = 0} (also together with \code{gradual_split = TRUE}) is used. Nevertheless, the coverage will be closer to 68 % than it would be without calibration. The prediction uncertainties at the beginning are conservative and become less conservative with increasing number of input points.
    #' @param use_reference_gp If TRUE, the covariance parameters determined for the GP in node 0 will be used for all subsequent GPs. The default is \code{FALSE}.
    #' @param min_abs_y_err Minimum absolute error assumed for y data. The default value is 0.
    #' @param min_rel_y_err Minimum relative error assumed for y data. The default value is \code{100 * .Machine$double.eps}.
    #' @param min_abs_node_pred_err Minimum absolute error on the prediction from a single node. The default value is 0.
    #' @param min_rel_node_pred_err Minimum relative error on the prediction from a single node. The default value is \code{100 * .Machine$double.eps}.
    #' @param prob_min_theta Minimum probability after which the overlap shape gets truncated (either towards 0 or 1). The default value is 0.01.
    #' @param add_buffer_in_prediction If TRUE, points in the data buffers are added to the GP before prediction. They are added into a temporarily created GP which contains the not yet included points. The GP in the node is not yet updated. The default is \code{FALSE}.
    #' @param x_dim Dimensionality of input points. It is set once the first point is received through the \code{update} method. It needs to be specified if \code{min_ranges} should be different from default.
    #' @param min_ranges Smallest allowed input data spread (per dimension) before node splitting stops. It is set to its default \code{min_ranges = rep(0.0, x_dim)} once the first point is received through the \code{update} method. \code{x_dim} needs to be specified by the user if it should be different from the default.
    #' @param max_cond_num Add additional noise if the covariance matrix condition number exceeds this value. The default is \code{NULL}.
    #' @param max_points The maximum number of points the tree is allowed to store. The default value is \code{Inf}.
    #' @return A new GPTree object. Tree-specific parameters are listed in this object. The field \code{nodes} contains a \link[hash]{hash} with all \link{GPNode}s and information related to nodes. The nodes in turn contain the local GPs. Nodes that have been split no longer contain a GP.
    #' @md
    #' @examples
    #' set.seed(42)
    #' ## Use the 1d toy data set from Higdon (2002)
    #' X <- as.matrix(sample(seq(0, 10, length.out = 31)))
    #' y <- sin(2 * pi * X / 10) + 0.2 * sin(2 * pi * X / 2.5)
    #' y_variance <- rep(0.1**2, 31)
    #' 
    #' ## Initialize a tree with Nbar = 15, retrain_buffer_length = 15, use_empirical_error = FALSE,
    #' ## and default parameters otherwise
    #' gptree <- GPTree$new(Nbar = 15, retrain_buffer_length = 15, use_empirical_error = FALSE)
    #' 
    #' ## For the purpose of this example, we simulate the data stream through a simple for loop.
    #' ## In actual applications, the input stream comes from e.g. a differential evolutionary scanner.
    #' ## We follow the procedure in the associated paper, thus letting the tree make a prediction
    #' ## first before we update the tree with the point.
    #' for (i in 1:nrow(X)) {
    #' y_pred_with_err = gptree$joint_prediction(X[i,], return_std = TRUE)
    #' ## Update the tree with the true (X,y) pair
    #' gptree$update(X[i,], y[i], y_variance[i])
    #' }
    #' 
    #' ## In the following, we go over different initializations of the tree
    #' ## 1. The same tree as before, but using the package mlegp:
    #' ## Note: since the default for gp_control is gp_control = list(covtype = "matern3_2"),
    #' ## we set gp_control to an empty list when using mlegp.
    #' gptree <- GPTree$new(Nbar = 15, retrain_buffer_length = 15, use_empirical_error = FALSE,
    #' wrapper = "mlegp", gp_control = list())
    #' 
    #' ## 2. Minimum working example:
    #' gptree <- GPTree$new()
    #' 
    #' ## 3. Fully specified example corresponding to the default settings
    #' ## Here, we choose to specify x_dim and min_ranges so that they correspond to the default values.
    #' ## If we do not specifiy them here, they will be automatically specified once
    #' ## the update or predict method is called.
    #' gptree <- GPTree$new(Nbar = 1000, retrain_buffer_length = 1000,
    #' gradual_split = TRUE, theta = 0, wrapper = "DiceKriging",
    #' gp_control = list(covtype = "matern3_2"),
    #' split_direction_criterion = "max_spread_per_lengthscale", split_position_criterion = "mean",
    #' shape_decay = "linear", use_empirical_error = TRUE, 
    #' use_reference_gp = FALSE, min_abs_y_err = 0, min_rel_y_err = 100 * .Machine$double.eps,
    #' min_abs_node_pred_err = 0, min_rel_node_pred_err = 100 * .Machine$double.eps,
    #' prob_min_theta = 0.01, add_buffer_in_prediction = FALSE, x_dim = ncol(X),
    #' min_ranges = rep(0.0, ncol(X)), max_cond_num = NULL, max_points = Inf)

    #' @export
    initialize = function(Nbar = 1000, retrain_buffer_length = Nbar, gradual_split = TRUE, theta = 0, wrapper = "DiceKriging", gp_control = list(covtype = "matern3_2"), split_direction_criterion = "max_spread_per_lengthscale", split_position_criterion = "median", shape_decay = "linear", use_empirical_error = TRUE, use_reference_gp = FALSE, min_abs_y_err = 0, min_rel_y_err = 100 * .Machine$double.eps, min_abs_node_pred_err = 0, min_rel_node_pred_err = 100 * .Machine$double.eps, prob_min_theta = 0.01, add_buffer_in_prediction = FALSE, x_dim = 0, min_ranges = NULL, max_cond_num = NULL, max_points = Inf) {
      self$Nbar <- Nbar
      self$split_direction_criterion <- split_direction_criterion
      self$split_position_criterion <- split_position_criterion
      self$shape_decay <- shape_decay
      self$prob_min_theta <- prob_min_theta
      self$theta <- theta
      self$x_dim <- x_dim
      self$min_abs_y_err <- min_abs_y_err
      self$min_rel_y_err <- min_rel_y_err
      self$min_abs_node_pred_err <- min_abs_node_pred_err
      self$min_rel_node_pred_err <- min_rel_node_pred_err
      self$use_empirical_error <- use_empirical_error
      self$retrain_buffer_length <- retrain_buffer_length
      self$add_buffer_in_prediction <- add_buffer_in_prediction
      self$min_ranges <- min_ranges
      self$max_cond_num <- max_cond_num
      self$max_points <- max_points

      self$wrapper <- wrapper
      self$gp_control <- gp_control

      self$use_reference_gp <- use_reference_gp
      
      self$gradual_split <- gradual_split

      self$nodes <- hash::hash()
      self$leaf_keys <- vector(mode = "character")

      self$n_points <- 0
      self$n_fed <- 0

      message(msgprefix, "Created a new GPTree (wrapper='", wrapper, "')")
      if (self$theta > 0 && self$use_empirical_error){
        message(msgprefix, "Warning: The coverage achieved through the uncertainty calibration cannot be guaranteed while theta > 0.")
      }
      
      # Add root node if x_dim and min_ranges have been specified
      if(self$x_dim != 0 && !(is.null(self$min_ranges))){
        message(msgprefix, "The user has specified these inputs for x_dim and min_ranges: x_dim: ", x_dim, ", min_ranges: ", paste(self$min_ranges, collapse=" "))
        # Check that min_ranges has the correct length
        if (length(self$min_ranges) != self$x_dim){
          stop(msgprefix, "The length of min_ranges does not match x_dim.")
        }
        
        # Create the root node
        self$add_node("0")
        
        # Get the wrapper-dependent n_points_train_limit from the WrappedGP instance in the root node
        n_points_train_limit <- self$nodes[["0"]]$wrapped_gp$n_points_train_limit
        
        # Require that Nbar >= 2 * n_points_train_limit
        if (self$Nbar < 2 * n_points_train_limit) {
          stop(msgprefix, "With the given wrapper and number of input dimensions, the lowest possible Nbar setting is ", 2 * n_points_train_limit, ".")
        }
        
        if (self$Nbar < 4 * n_points_train_limit) {
          message(msgprefix, "Nbar (", self$Nbar, ") is close to 2 * n_points_train_limit (2 * ", n_points_train_limit, " = ", 2 * n_points_train_limit, "). Expect that the auto-balancer activates.")
        }
      }
    },


    #' @description 
    #' Add a new GPNode to the tree. IS EXPECTED TO NOT BE CALLED BY THE USER
    #' @param key Key of the new leaf
    add_node = function(key) {
      self$leaf_keys[[length(self$leaf_keys) + 1]] <- key
      message(msgprefix, "Creating a new node at ", key, ".")

      self$nodes[[key]] <- GPNode$new(
        key,
        self$x_dim,
        self$theta,
        self$split_direction_criterion,
        self$split_position_criterion,
        self$shape_decay,
        self$prob_min_theta,
        self$Nbar,
        self$wrapper,
        self$gp_control,
        self$retrain_buffer_length,
        self$add_buffer_in_prediction,
        min_ranges = self$min_ranges,
        is_leaf = TRUE
      )
    },


    #' @description 
    #' Marginal probability for point x to belong to node with given key. IS EXPECTED TO NOT BE CALLED BY THE USER
    #' @param x Single input data point from the data stream; has to be a vector with length equal to x_dim
    #' @param key Key of the node
    #' @return Returns the marginal probability for point x to belong to node with given key
    get_marginal_point_prob = function(x, key) {
      p <- 1.0
      depth <- nchar(key)

      # The key identifies a unique path through the tree.
      # Walk down this path and multiply the probabilities to get the total probability
      # that the point x should have ended up in node with the given key
      current_level <- 2
      while (p > 0. && current_level <= depth) {
        current_key <- substr(key, 1, current_level - 1)
        child_index <- substr(key, current_level, current_level)

        new_p <- self$nodes[[current_key]]$get_prob_child_1(x)

        if (child_index == "0") {
          p <- p * (1.0 - new_p)
        } else if (child_index == "1") {
          p <- p * new_p
        }

        current_level <- current_level + 1
      }

      return(p)
    },
    
    
    
    #' @description 
    #' Assigns the given input point x with target variable y and associated variance y_var to a node and updates the tree accordingly
    #' @details The methods takes care of both updating an existing node and splitting the parent node into two child nodes. It ensures that the each child node has at least \code{n_points_train_limit} in each GP. Further handling of duplicate points is also done here.
    #' @param x Most recent single input data point from the data stream; has to be a vector with length equal to x_dim
    #' @param y Value of target variable at input point x; has to be a one-dimensional matrix or a vector; any further columns will be ignored
    #' @param y_var Variance of the target variable; has to be a one-dimensional matrix or vector
    #' @param retrain_node If TRUE, the GP node will be retrained after the point is added.
    update = function(x, y, y_var = 0., retrain_node = TRUE) {
      message(msgprefix, "Got (x,y) point: x: ", paste(x,collapse=" "), " y: ", y, " y_var: ", y_var)

      # Finish the initialization if x_dim has not yet been specified
      # Set x_dim and min_ranges based on the dimensions of the first input point
      if (self$x_dim == 0 || is.null(self$min_ranges)){
        self$x_dim <- length(x)
        self$min_ranges = rep(0.0, self$x_dim)
        
        # Create the root node
        self$add_node("0")
        
        # Get the wrapper-dependent n_points_train_limit from the WrappedGP instance in the root node
        n_points_train_limit <- self$nodes[["0"]]$wrapped_gp$n_points_train_limit
        
        # Require that Nbar >= 2 * n_points_train_limit
        if (self$Nbar < 2 * n_points_train_limit) {
          stop(msgprefix, "With the given wrapper and number of input dimensions, the lowest possible Nbar setting is ", 2 * n_points_train_limit, ".")
        }
        
        if (self$Nbar < 4 * n_points_train_limit) {
          message(msgprefix, "Nbar (", self$Nbar, ") is close to 2 * n_points_train_limit (2 * ", n_points_train_limit, " = ", 2 * n_points_train_limit, "). Expect that the auto-balancer activates.")
        }
      }
      self$n_fed <- self$n_fed + 1
      point_id <- self$n_fed

      # Check if max_points has been reached.
      if (self$n_points >= self$max_points)
      {
        message(msgprefix, "The max number of stored points has been reached (max_points = ", self$max_points, "). Will not add more points.")
        return(FALSE)
      }

      # Update the nodes.
      # Use the chain of node probabilities to pick a leaf node
      # for the input point. We'll refer to the chosen node as current_node.
      key <- "0"
      current_node <- self$nodes[[key]]
      while (!current_node$is_leaf) {
        p <- current_node$get_prob_child_1(x)
        child_index <- as.character(rbinom(1, 1, p))
        key <- paste(key, child_index, sep = "")
        current_node <- self$nodes[[key]]
      }

      # If this node is not allowed to split further, don't add more points to it.
      if (!any(current_node$can_split)) {
        message(msgprefix, "Will not add more points to node ", current_node$key, ". Its input space is smaller than the limit in min_ranges.")
        return(FALSE)
      }

      # If the given y_var is (0 or) too small, set it to the minimum allowed error
      y_err_limit <- max(self$min_rel_y_err * abs(y), self$min_abs_y_err)
      y_var <- max(y_var, y_err_limit**2)

      # Register residual if gradual splitting is activated
      # Register only if a GP has already been created
      if (self$use_empirical_error && !(is.null(current_node$wrapped_gp$gp))) {
        current_node$register_residual(x, y)
        current_node$update_empirical_error_pars()
      }

      # Add point to the current node
      current_node$wrapped_gp$store_point(x, y, y_var, shared=FALSE, remove_shared=TRUE)
      current_node$point_ids <- append(current_node$point_ids, point_id)

      # Root node is only trained when it is full, regardless of n_points_train_limit
      if (current_node$key == "0") {
        retrain_node <- FALSE
      }

      # If the current node is *not* full yet, add the new point and retrain.
      # If the current node *is* full, we first split it into two child nodes,
      # add the new point to one of the child nodes, and retrain both the new nodes.
      if (current_node$wrapped_gp$n_points < current_node$Nbar) {

        # Retrain the node?
        if (retrain_node) {

          did_train <- current_node$wrapped_gp$train()

          if (did_train) {
            message(msgprefix, "Trained GP in node ", current_node$key, ".")

            # Update add_y_var in the wrapped GP
            if (!is.null(self$max_cond_num)) {
              current_node$wrapped_gp$update_add_y_var(self$max_cond_num)
            }
          }
        }
      } else { # The current node is full, let's split it and train the new child nodes

        # For the very first node, make sure to retrain just before splitting
        if (current_node$key == "0") {
          current_node$wrapped_gp$train(do_buffer_check = FALSE)
        }

        # First we need to compute the probability parameters for the current node.
        # This decides how points will be distributed to the two child nodes
        # we'll create below
        did_update <- current_node$update_prob_pars()
        if (!did_update) {
          message(msgprefix, "Will not split the node ", current_node$key, ". Its input space is smaller than the limit in min_ranges.")
          return(FALSE)
        }

        # The current node is no longer a leaf
        current_node$is_leaf <- FALSE

        # Remove the key from the list of leaf keys
        self$leaf_keys <- self$leaf_keys[!self$leaf_keys %in% c(key)]

        # Create two child nodes
        child_0_key <- paste(key, "0", sep = "")
        child_1_key <- paste(key, "1", sep = "")
        self$add_node(child_0_key)
        self$add_node(child_1_key)

        # Get references to the new child nodes to avoid having
        # to look them up repeatedly
        child_0 <- self$nodes[[child_0_key]]
        child_1 <- self$nodes[[child_1_key]]

        # The child nodes inherit the empirical error parameters
        if (self$use_empirical_error) {
          child_0$residuals <- current_node$residuals
          child_0$pred_errs <- current_node$pred_errs
          child_0$error_scaler <- current_node$error_scaler

          child_1$residuals <- current_node$residuals
          child_1$pred_errs <- current_node$pred_errs
          child_1$error_scaler <- current_node$error_scaler
        }
        
        # Copy over some parameters from the GP of the current node to
        # GPs of the child nodes
        child_0$wrapped_gp$init_covpars <- current_node$wrapped_gp$init_covpars
        child_1$wrapped_gp$init_covpars <- current_node$wrapped_gp$init_covpars
        child_0$wrapped_gp$add_y_var <- current_node$wrapped_gp$add_y_var
        child_1$wrapped_gp$add_y_var <- current_node$wrapped_gp$add_y_var
        
        # If we are reusing the covariance parameters from the GP in the parent node
        # (going all the way back to node 0), tell the child GPs that they are not
        # allowed to estimate these parameters.
        if (self$use_reference_gp) {
          child_0$wrapped_gp$estimate_covpars <- FALSE
          child_1$wrapped_gp$estimate_covpars <- FALSE
        }
        
        # Get a table with info on how the data points in the current node
        # should be distributed between the two child nodes.
        data_split_table <- self$get_data_split_table(current_node)

        # Get all the data from the current node
        X_data <- current_node$wrapped_gp$get_X_data()
        y_data <- current_node$wrapped_gp$get_y_data()
        y_var_data <- current_node$wrapped_gp$get_y_var_data()
        point_ids <- current_node$point_ids

        # Check how many points are assigned to each child node
        n_node_0 <- length(which(data_split_table[, 4] == 0))
        n_node_1 <- length(which(data_split_table[, 4] == 1))

        # Now add the points to the child nodes.

        # First, if we are doing gradual splitting, provide each child the points of the other
        if (self$gradual_split) {

          split_dimension <- current_node$split_index

          # Give the child_0 points to child_1:
          for (i in which(data_split_table[, 4] == 0)[1:n_node_0]) {
            child_1$wrapped_gp$store_point(X_data[i, ], y_data[i], y_var_data[i], shared=TRUE, remove_shared=FALSE)
          }
          ordering <- order(child_1$wrapped_gp$X_shared[,split_dimension], decreasing=TRUE)
          child_1$wrapped_gp$X_shared <- child_1$wrapped_gp$X_shared[ordering,]
          child_1$wrapped_gp$y_shared <- matrix(child_1$wrapped_gp$y_shared[ordering], ncol=1)
          child_1$wrapped_gp$y_var_shared <- matrix(child_1$wrapped_gp$y_var_shared[ordering], ncol=1)

          # Give the child_1 points to child_0:
          for (i in which(data_split_table[, 4] == 1)[1:n_node_1]) {
            child_0$wrapped_gp$store_point(X_data[i, ], y_data[i], y_var_data[i], shared=TRUE, remove_shared=FALSE)
          }
          ordering <- order(child_0$wrapped_gp$X_shared[,split_dimension], decreasing=FALSE)
          child_0$wrapped_gp$X_shared <- child_0$wrapped_gp$X_shared[ordering,]
          child_0$wrapped_gp$y_shared <- matrix(child_0$wrapped_gp$y_shared[ordering], ncol=1)
          child_0$wrapped_gp$y_var_shared <- matrix(child_0$wrapped_gp$y_var_shared[ordering], ncol=1)

        }

        # Now provide each child node with its "own" points and create the GPs.
        for (i in which(data_split_table[, 4] == 0)[1:n_node_0]) {
          child_0$wrapped_gp$store_point(X_data[i, ], y_data[i], y_var_data[i], shared=FALSE, remove_shared=FALSE)
          child_0$point_ids <- append(child_0$point_ids, point_ids[i])
        }
        child_0$wrapped_gp$train(do_buffer_check = FALSE)
        
        for (i in which(data_split_table[, 4] == 1)[1:n_node_1]) {
          child_1$wrapped_gp$store_point(X_data[i, ], y_data[i], y_var_data[i], shared=FALSE, remove_shared=FALSE)
          child_1$point_ids <- append(child_1$point_ids, point_ids[i])
        }

        # If we're using gradual split, child_1 starts out as identical to child_0,
        # so we can use the same covariance parameters
        if (self$gradual_split) {
          child_1$wrapped_gp$estimate_covpars <- FALSE
          child_1$wrapped_gp$init_covpars <- child_0$wrapped_gp$init_covpars
          child_1$wrapped_gp$train(do_buffer_check = FALSE)
          child_1$wrapped_gp$estimate_covpars <- child_0$wrapped_gp$estimate_covpars
        } else {
          child_1$wrapped_gp$train(do_buffer_check = FALSE)
        } 
        
        # Update the add_y_var parameters in the new wrapped GPs
        if (!is.null(self$max_cond_num)) {
          child_0$wrapped_gp$update_add_y_var(self$max_cond_num)
          child_1$wrapped_gp$update_add_y_var(self$max_cond_num)
        }

        # Finished splitting

        # Print info about the x data ranges of the new nodes
        child_0_range_str <- paste("x data ranges for node ", child_0_key, ":", sep = "")
        child_1_range_str <- paste("x data ranges for node ", child_1_key, ":", sep = "")

        for (dim_index in 1:self$x_dim) {
          col_name <- colnames(child_0$wrapped_gp$get_X_data())[dim_index]

          child_0_xmin <- min(child_0$wrapped_gp$get_X_data()[, dim_index])
          child_0_xmax <- max(child_0$wrapped_gp$get_X_data()[, dim_index])

          child_1_xmin <- min(child_1$wrapped_gp$get_X_data()[, dim_index])
          child_1_xmax <- max(child_1$wrapped_gp$get_X_data()[, dim_index])

          child_0_range_str <- paste(child_0_range_str, " ", col_name, ": [", child_0_xmin, ",", child_0_xmax, "]", sep = "")
          child_1_range_str <- paste(child_1_range_str, " ", col_name, ": [", child_1_xmin, ",", child_1_xmax, "]", sep = "")
        }

        message(msgprefix, child_0_range_str)
        message(msgprefix, child_1_range_str)

        # Make sure that both child GPs have been created
        if (is.null(child_0$wrapped_gp$gp) | is.null(child_1$wrapped_gp$gp)) {
          stop(msgprefix, "One or both child GPs were not created. This should never happen.")
        }
        
        # Now we can delete the old parent GP
        current_node$delete_gp()
      }

      # Register the new point as included in the tree
      self$n_points <- self$n_points + 1

      return(TRUE)
    },


    #' @description 
    #' Generates a table used to distribute data points from a node to two child nodes
    #' @param current_node The GPNode whose data should be distributed
    #' @return A matrix object
    get_data_split_table = function(current_node) {

      # We must distribute the data points between the two child nodes in such a way
      # that both child nodes have a sufficient number of points. To determine a
      # data division that works, we create a table that for each data point contains
      # the point index, the point coordinate along the (transformed) split dimension, 
      # the probability for it to belong to child 1, and the key for the child node it
      # is currently assigned to. Then we modify this table as needed until we have an
      # acceptable point distribution that we can use to create the two child nodes.

      X_data <- current_node$wrapped_gp$get_X_data()
      X_data_transformed <- X_data
      if (current_node$use_pc_transform) {
        X_data_transformed <- current_node$transform(X_data)
      }

      nrows <- nrow(X_data)  # (This should equal Nbar)

      split_dimension <- current_node$split_index
      child_1_probs <- apply(X_data, 1, current_node$get_prob_child_1)
      data_split_table <- matrix( c(1:nrows,
                                   X_data_transformed[,split_dimension],
                                   child_1_probs,
                                   apply(as.matrix(child_1_probs), 1, rbinom, n = 1, size = 1)),
                                 ncol=4,
                               )
      
      # Ensure that there are at least n_points_train_limit points assigned to each node.
      n_node_0 <- length(which(data_split_table[, 4] == 0))
      n_node_1 <- length(which(data_split_table[, 4] == 1))
      if (n_node_0 < current_node$wrapped_gp$n_points_train_limit) {
        message(msgprefix, "Not enough points in child node 0. Will assign it points from child node 1.")
        # Sort first by the probabilities, then break ties by sorting by coordinate values
        data_split_table <- (data_split_table[order(data_split_table[, 3], data_split_table[, 2]), ])
        # Find the closest points in the other leaf and move them over
        i <- 1
        while (n_node_0 < current_node$wrapped_gp$n_points_train_limit) {
          if (data_split_table[i, 4] != 0) {
            data_split_table[i, 4] <- 0
            n_node_0 <- n_node_0 + 1
            n_node_1 <- n_node_1 - 1
          }
          i <- i + 1
        }
      } else if (n_node_1 < current_node$wrapped_gp$n_points_train_limit) {
        message(msgprefix, "Not enough points in child node 0. Will assign it points from child node 0.")
        # Sort first by the probabilities, then break ties by sorting by coordinate value
        data_split_table <- (data_split_table[order(data_split_table[, 3], data_split_table[, 2], decreasing = TRUE), ])
        # Find the closest points in the other leaf and move them over
        i <- 1
        while (n_node_1 < current_node$wrapped_gp$n_points_train_limit) {
          if (data_split_table[i, 4] != 1) {
            data_split_table[i, 4] <- 1
            n_node_0 <- n_node_0 - 1
            n_node_1 <- n_node_1 + 1
          }
          i <- i + 1
        }
      }

      # Re-order the table for the next step
      data_split_table <- (data_split_table[order(data_split_table[, 1], decreasing = FALSE), ])

      return(data_split_table)
    },


    #' @description 
    #' Compute the joint prediction from all relevant leaves for an input point x
    #' @details We follow Eqs. (5) and (6) in \href{https://arxiv.org/abs/2006.09446}{this paper}
    #' @param x Single data point for which the predicted joint mean (and standard deviation) is computed; has to be a vector with length equal to x_dim
    #' @param return_std If TRUE, the standard error of the prediction is returned
    #' @return The prediction (and its standard error) for input point x from this tree
    joint_prediction = function(x, return_std = TRUE) {

      message(msgprefix, "Got x point: ", paste(x, collapse=" "))
      
      # Finish the initialization if x_dim has not yet been specified
      # Set x_dim and min_ranges based on the dimensions of the first input point
      if (self$x_dim == 0 || is.null(self$min_ranges)){
        self$x_dim <- length(x)
        self$min_ranges = rep(0.0, self$x_dim)
        
        # Create the root node
        self$add_node("0")
        
        # Get the wrapper-dependent n_points_train_limit from the WrappedGP instance in the root node
        n_points_train_limit <- self$nodes[["0"]]$wrapped_gp$n_points_train_limit
        
        # Require that Nbar >= 2 * n_points_train_limit
        if (self$Nbar < 2 * n_points_train_limit) {
          stop(msgprefix, "With the given wrapper and number of input dimensions, the lowest possible Nbar setting is ", 2 * n_points_train_limit, ".")
        }
        
        if (self$Nbar < 4 * n_points_train_limit) {
          message(msgprefix, "Nbar (", self$Nbar, ") is close to 2 * n_points_train_limit (2 * ", n_points_train_limit, " = ", 2 * n_points_train_limit, "). Expect that the auto-balancer activates.")
        }
      }

      joint_y_pred <- 0
      joint_y_pred_var <- 0

      # Preallocate two vectors, one for leaf keys and one for marginal probabilities.
      # These will be filled when the nested function collect_leaves below is run.
      n_leaves_tot <- length(self$leaf_keys)
      pred_leaf_keys <- character(n_leaves_tot)
      pred_leaf_probs <- numeric(n_leaves_tot)

      # Keep track of the number of contributing leaves, so we can truncate
      # pred_leaf_keys and pred_leaf_probs when we're done collecting leaves.
      n_collected <- 0
      sum_probs <- 0
      collection_done <- FALSE

      # Recursive function to collect contributing leaves. This function writes to
      # the variables pred_leaf_keys, pred_leaf_probs, n_collected, sum_probs and
      # collection_done. It always returns NULL.
      collect_leaves <- function(x, current_key, current_prob) {

        if(collection_done | current_prob <= 0) {
          return(NULL)
        }

        # Get the current node
        current_node <- self$nodes[[current_key]]

        # Return if we have reached a leaf node
        if(current_node$is_leaf) {
          n_collected <<- n_collected + 1
          pred_leaf_keys[n_collected] <<- current_key
          pred_leaf_probs[n_collected] <<- current_prob

          sum_probs <<- sum_probs + current_prob
          if(sum_probs >= 1) {
            collection_done <<- TRUE
          }

          return(NULL)
        }

        new_p <- current_node$get_prob_child_1(x)

        p0 <- current_prob * (1 - new_p)
        child0_key = paste(current_key, "0", sep="")
        collect_leaves(x, child0_key, p0)

        p1 <- current_prob * new_p
        child1_key = paste(current_key, "1", sep="")
        collect_leaves(x, child1_key, p1)

        return(NULL)
      }

      # Collect the leaves
      collect_leaves(x, "0", 1)

      # Truncate the vectors
      length(pred_leaf_keys) <- n_collected
      length(pred_leaf_probs) <- n_collected

      # Now compute the joint prediction
      if(return_std) {

        for(i in 1:n_collected) {

          leaf_key <- pred_leaf_keys[i]
          leaf_prob <- pred_leaf_probs[i]

          pred_with_err <- self$nodes[[leaf_key]]$wrapped_gp$predict(x, return_std=TRUE)
          y_pred <- pred_with_err[1]
          y_pred_std <- pred_with_err[2]

          # If the prediction uncertainty is 0 or negative, something has gone
          # wrong and we should assign a large uncertainty to this prediction
          if (y_pred_std <= 0) {
            message(msgprefix, "Got a negative uncertainty from the GP. Will assign a 100% uncertainty to this prediction.")
            y_pred_std <- abs(y_pred)
          }

          # Modify prediction with the empirical_error method?
          if (self$use_empirical_error) {
            y_pred_std <- y_pred_std * self$nodes[[leaf_key]]$error_scaler
          }

          # Enforce minimum node prediction uncertainty
          y_pred_std_limit <- max(self$min_abs_node_pred_err, abs(y_pred) * self$min_rel_node_pred_err)
          y_pred_std <- max(y_pred_std, y_pred_std_limit)

          message(msgprefix, "Prediction from node ", leaf_key, ": ", y_pred, " +- ", y_pred_std, ", leaf_prob: ", leaf_prob)

          # Add this leaf's contribution to the joint prediction
          joint_y_pred <- joint_y_pred + leaf_prob * y_pred
          joint_y_pred_var <- joint_y_pred_var + leaf_prob * (y_pred**2 + y_pred_std**2)
        }

        # Last part of Eq. (6)
        joint_y_pred_var <- joint_y_pred_var - joint_y_pred**2

        message(msgprefix, "Joint prediction from the tree: ", joint_y_pred, " +- ", sqrt(joint_y_pred_var))

      } else {

        for(i in 1:n_collected) {

          leaf_key <- pred_leaf_keys[i]
          leaf_prob <- pred_leaf_probs[i]

          y_pred <- self$nodes[[leaf_key]]$wrapped_gp$predict(x, return_std=FALSE)

          # Add this leaf's contribution to the joint prediction
          joint_y_pred <- joint_y_pred + leaf_prob * y_pred
        }

        message(msgprefix, "Joint prediction from the tree: ", joint_y_pred)
      }


      # Return result
      result <- c()

      if (return_std) {
        result <- c(joint_y_pred, sqrt(joint_y_pred_var))
      } else {
        result <- c(joint_y_pred)
      }

      return(result)
    }

  )
)
