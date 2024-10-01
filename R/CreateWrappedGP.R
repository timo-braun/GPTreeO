# The WrappedGP familiy of classes
#
# The Wrapped<PACKAGE>GP classes are wrapper classes, intended to contain a GP object specific to some
# GP package (e.g. DiceKriging) and provide a set of standard methods to allow GPNode and GPTree to interact
# with the underlying GP. We do this to separate the DLGP implementation (trees, nodes, overlaps, point caps, ...)
# from the details of the underlying GP implementation (choice of kernel, training procedure, ...)
#
# The following methods are expected by GPNode and/or GPTree:
#
#   - initialize((X, y, y_var, gp_control, init_covpars, retrain_buffer_length, add_buffer_in_prediction, estimate_covpars=TRUE, X_shared=NULL, y_shared=NULL, y_var_shared=NULL))
#   - update_init_covpars()
#   - get_lengthscales()
#   - get_X_data(include_shared = FALSE)
#   - get_y_data(include_shared = FALSE)
#   - get_y_var_data(include_shared = FALSE)
#   - get_cov_mat()
#   - update_add_y_var(max_cond_num)
#   - store_point(x, y, y_var, shared = FALSE, remove_shared = TRUE)
#   - delete_buffers()
#   - train(do_buffer_check = TRUE)
#   - predict(x, return_std=TRUE)
#   - delete_gp()
#
#   optional methods specifically used in WrappedDiceKrigingGP and WrappedmlegpGP (replace PACKAGE with DiceKriging or mlegp):
#   - create_PACKAGE_gp(X, y, y_var)
#   - call_PACKAGE_predict(x, use_gp = NULL)
#   
#   The following functions explicitly depend on the chosen GP package (excluding create_PACKAGE_gp and call_PACKAGE_predict)
#   - get_cov_mat()
#   - predict(x, return_std=TRUE)
#   - the list provided in gp_control
# 
# Also, the user needs to ensure that the GP parameters are retrieved correctly in the code. For example, the parameter "param" for a GP of class S3 is stored as gp$param, while for S4 class it is gp@param .



msgprefix <- paste(packageName(), ": ", sep = "")

#' Factory function called by GPNode to create the wrapper for a specified GP package
#'
#' @details A detailed list of expected functions from GPTree and GPNode can be found in the comments of this file. Currently, GPs from the \code{DiceKriging} package (\link{WrappedDiceKrigingGP}) and \code{mlegp} package (\link{WrappedmlegpGP}) are implemented. The user can create their own wrapper using \link{WrappedGP}.
#' @param wrapper A string specifying what GP implementation is used
#' @param X Input data matrix with x_dim columns and at maximum Nbar rows. Is used to create the first iteration of the local GP.
#' @param y Value of target variable at input point x; has to be a one-dimensional matrix or a vector; any further columns will be ignored
#' @param y_var Variance of the target variable; has to be a one-dimensional matrix or vector
#' @param gp_control A list of GP implementation-specific options, passed directly to the wrapped GP implementation
#' @param init_covpars Initial covariance parameters of the local GP
#' @param retrain_buffer_length Only retrain when the number of buffer points or collected points exceeds this value
#' @param add_buffer_in_prediction If TRUE, points in the data buffers are added to the GP before prediction. They are added into a temporarily created GP which contains the not yet included points. The GP in the node is not yet updated.
#' @return The wrapper of the chosen GP package, containing the respective GP and information on the shared points and those stored in the buffer.
#' @export
CreateWrappedGP <- function(wrapper, X, y, y_var, gp_control, init_covpars, retrain_buffer_length, add_buffer_in_prediction) {
  if (wrapper == "mlegp") {
    return(WrappedmlegpGP$new(X, y, y_var, gp_control, init_covpars, retrain_buffer_length, add_buffer_in_prediction, estimate_covpars=TRUE))
  } else if (wrapper == "DiceKriging") {
    return(WrappedDiceKrigingGP$new(X, y, y_var, gp_control, init_covpars, retrain_buffer_length, add_buffer_in_prediction, estimate_covpars=TRUE))
  } else if (wrapper == "dummy") {
    return(WrappedGP$new(X, y, y_var, gp_control, init_covpars, retrain_buffer_length, add_buffer_in_prediction, estimate_covpars=TRUE))
  } else {
    stop(msgprefix, "Received an unknown value for the 'wrapper' argument.")
  }
}
