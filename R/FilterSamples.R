################# ----------- Filter Samples Functions ----------- ###################

# Function to remove a specific sample from SE (by column and value)



# Function to remove reference samples



################# ----------- Plotting Functions ----------- ###################

#' Plot number of non-zero proteins per sample
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which data type should be used (default raw)
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#'
#' @return ggplot object
#' @export
#'
plot_nr_prot_samples <- function(se, ain="raw", color_by = NULL, label_by = NULL){
  # get color_by
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
  }

  # get label_by
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
  }

  # prepare data
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt$ID <- rownames(dt)
  dt <- data.table::melt(dt, id.vars = "ID", variable.name = "Column", value.name = "Intensity")
  dt$Intensity[is.na(dt$Intensity)] <- 0
  dt <- dt[dt$Intensity != 0,]
  dt <- dt %>% dplyr::count(Column)
  dt <- merge(dt, data.table::as.data.table(SummarizedExperiment::colData(se)), by="Column")
  dt <- data.table::as.data.table(dt)

  # order data
  if(!is.null(color_by)){
    dt[, color_by] <- factor(dt[,color_by], levels = unique(dt[,color_by]))
    dt <- dt[order(dt[,color_by]),]
    if(show_sample_names){
      dt[, label_by] <- factor(dt[,label_by], levels = unique(dt[,label_by]))
    } else {
      dt$Column <- factor(dt$Column, levels = dt$Column)
    }
  }

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector[54:73])

  # plot
  if(is.null(color_by)){
    p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(label_by), y=get("n"))) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      ggplot2::labs(x = "Samples", y ="Number of proteins \n (non-zero protein abundance)")
  } else {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(label_by), y=get("n"), fill=get(color_by))) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      ggplot2::scale_fill_manual(name = color_by, values = col_vector ) +
      ggplot2::labs(x = "Samples", y ="Number of proteins \n (non-zero protein abundance)")
  }
  if(!show_sample_names){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  }
  return(p)
}

# Plot Total Protein Intensity


# Outlier Detection Using POMA



# remove outliers
