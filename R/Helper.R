

#' Helper function to check the condition value
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#'
#' @return String of column for condition
#' @export
#'
get_condition_value <- function(se, condition){
  if(is.null(condition)){
    # check if condition in metadata of se
    if(is.null(S4Vectors::metadata(se)$condition)){
      # print error
      stop("No condition provided!")
    } else {
      condition <- S4Vectors::metadata(se)$condition
      message("Condition of SummarizedExperiment used!")
    }
  } else {
    # check if condition in colData
    if(! condition %in% colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))){
      stop(paste0("No column named ", condition, " in SummarizedExperiment!"))
    }
  }
  return(condition)
}


#' Helper function to get correct value for coloration of plots (color_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#'
#' @return String of column to color or NULL if no color should be applied
#' @export
#'
get_color_value <- function(se, color_by){
  if(is.null(color_by)){
    # check if condition in metadata of se
    if(is.null(S4Vectors::metadata(se)$condition)){
      color_by <- NULL
      message("No condition provided. Hence, no color bar added.")
    } else {
      color_by <- S4Vectors::metadata(se)$condition
      message("Condition of SummarizedExperiment used!")
    }
  } else if(color_by == "No"){
    color_by <- NULL
    message("No color bar added.")
  } else {
    #check if color_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!color_by %in% cols){
      stop("Color_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(color_by)
}



#' Helper function to get correct value for sample labeling of plots (label_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#'
#' @return String of column to label or NULL if no label should be applied
#' @export
#'
get_label_value <- function(se, label_by){
  show_sample_names <- TRUE
  if(is.null(label_by)){
    # check if condition in metadata of se
    if(is.null(S4Vectors::metadata(se)$label)){
      label_by <- "Column"
      show_sample_names <- FALSE
      message("No label provided. Hence, no labeling of samples.")
    } else {
      label_by <- S4Vectors::metadata(se)$label
      message("Label of SummarizedExperiment used!")
    }
  } else if (label_by == "No") {
    label_by <- "Column"
    show_sample_names <- FALSE
    message("No labeling of samples.")
  } else {
    # check if label_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!label_by %in% cols){
      stop("Label_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(list(show_sample_names, label_by))
}

#' Helper function to get correct value for shaping of plots (shape_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param shape_by String specifying the column to shape the samples (If NULL or "No", no shaping is done.)
#'
#' @return String of column to shape or NULL if no shaping should be done
#' @export
#'
get_shape_value <- function(se, shape_by){
  if(is.null(shape_by)){
      message("No shaping done.")
  } else if(shape_by == "No"){
    shape_by <- NULL
    message("No shaping done.")
  } else {
    # check if shape_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!shape_by %in% cols){
      stop("Shape_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(shape_by)
}


#' Helper function to get correct value for faceting of plots (facet_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param facet_by String specifying the column to facet the samples (If NULL or "No", no faceting is done.)
#'
#' @return String of column to facet or NULL if no faceting should be done
#' @export
#'
get_facet_value <- function(se, facet_by){
  if(is.null(facet_by)){
    message("No faceting done.")
  } else if(facet_by == "No"){
    facet_by <- NULL
    message("No faceting done.")
  } else {
    # check if facet_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!facet_by %in% cols){
      stop("Facet_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(facet_by)
}



#' Helper function to check whether all given assays are in SummarizedExperiment object
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#'
#' @return NULL if no methods in SummarizedExperiment object, else all available methods ready for visualization
#' @export
#'
check_input_assays <- function(se, ain) {
  if (is.null(ain)) {
    assays_se <- names(SummarizedExperiment::assays(se))
    message("All assays of the SummarizedExperiment will be used.")
    return(assays_se)
  } else {
    assays_se <- names(SummarizedExperiment::assays(se))
    not_existing_assays <- ain[!ain %in% assays_se]
    existing_assays <- ain[ain %in% assays_se]
    if (length(not_existing_assays) > 0) {
      # some methods of ain are not in SummarizedExperiment object
      # check if there are still some methods of ain in the SummarizedExperiment object
      if (length(existing_assays) > 0) {
        warning(
          paste0(
            paste0(not_existing_assays, collapse = ", "),
            " not in SummarizedExperiment object. Check with names(assays(se)) which assays can be used for visualization!"
          )
        )
        message(paste0(
          paste0(existing_assays, collapse = ", "),
          " used for visualization."
        ))
        return(existing_assays)
      } else {
        stop(
          paste0(
            paste0(not_existing_assays, collapse = ", "),
            " not in SummarizedExperiment object. Check with names(assays(se)) which assays can be used for visualization!"
          )
        )
        return(NULL)
      }
    } else {
      return(ain)
    }
  }
}
