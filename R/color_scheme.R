#' Color Scheme
#'
#' This function returns a predefined color scheme based on the input number.
#'
#' @param scheme_num An integer from 1 to 5 indicating the desired color scheme.
#'
#' @return A vector of color codes representing the selected color scheme.
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(your_package_name)
#'
#' # Get color scheme 3
#' colors <- color_scheme(3)
#' }
#' @export
color_scheme <- function(scheme_num=1) {
  # Define a list of color schemes
  color_schemes <- list(
    scheme_1 = c('#08306B', '#2171B5', '#6BAED6', '#9ECAE1'),
    scheme_2 = c('mediumpurple1', 'mediumpurple3', 'slateblue3', 'slateblue4'),
    scheme_3 = c('darkorchid4', 'darkorchid3', 'orchid3', 'orchid1'),
    scheme_4 = c('green', 'limegreen', 'forestgreen', 'darkgreen'),
    scheme_5 = c('chocolate1', 'coral2', 'firebrick3', 'firebrick4')
  )
  
  # Check if the scheme_num is valid
  if (scheme_num < 1 || scheme_num > length(color_schemes)) {
    stop("Invalid scheme number. Please choose a number from 1 to 5.")
  }
  
  # Return the selected color scheme
  return(color_schemes[[scheme_num]])
}
