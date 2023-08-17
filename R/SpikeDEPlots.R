
#' Get performance metrics of DE results of spike-in data set.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#'
#' @return data table with multiple performance metrics of the DE results
#' @export
#'
get_spiked_stats_DE <- function(se, de_res){
  spike_column <- S4Vectors::metadata(se)$spike_column
  spike_val <- S4Vectors::metadata(se)$spike_value
  stats <- de_res %>% dplyr::group_by(Assay, Comparison) %>%
    dplyr::summarise(TP = sum(Change %in% c("Up Regulated", "Down Regulated", "Significant Change") & get(spike_column) == spike_val, na.rm=TRUE),
                     FP = sum(Change %in% c("Up Regulated", "Down Regulated", "Significant Change") & get(spike_column) != spike_val, na.rm=TRUE),
                     FN = sum(Change == "No Change" & get(spike_column) == spike_val, na.rm=TRUE),
                     TN = sum(Change == "No Change" & get(spike_column) != spike_val, na.rm=TRUE),
                     Sensitivity = TP / (TP + FN),
                     Specificity = TN /(TN + FP),
                     Precision = TP / (TP + FP),
                     FPR = FP / (FP + TN),
                     F1Score = 2 * (Precision * Sensitivity) / (Precision + Sensitivity),
                     Accuracy = (TP + TN)/(TP + TN + FP + FN)) %>%
    data.table::as.data.table()
  return(stats)
}


#' Heatmap of performance metrics for spike-in data sets
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#' @param metrics vector of Strings specifying the metrics (must be colnames of stats)
#'
#' @return
#' @export
#'
#' @examples
plot_stats_spiked_heatmap <- function(stats, ain = NULL, comparisons = NULL, metrics = c("Accuracy", "Precision", "F1Score")){
  # check if measures in stats
  if(is.null(metrics)){
    message("All available metrics will be used for plotting")
    valid_metrics <- colnames(stats)[!colnames(stats) %in% c("TP", "FN", "TN", "FP")]
  } else {
    metrics <- c(metrics)
    valid_metrics <- metrics[metrics %in% colnames(stats)]
    non_valid_metrics <- metrics[!metrics %in% colnames(stats)]
    if(length(valid_metrics)== 0){
      stop("No valid metrics! Check colnames of stats.")
    } else if(length(non_valid_metrics)>0){
      warning(paste0("Not all metrics are valid!", metrics, " are not visualized!"))
    }
  }
  metrics <- valid_metrics
  dt <- stats[, colnames(stats) %in% c("Assay", "Comparison", metrics)]

  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]

  # prepare data
  melted_dt <- data.table::melt(dt, variable.name = "Metric", value.name = "Value", id.vars = c("Assay", "Comparison"))
  p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("Assay"), y = get("Comparison"), fill = get("Value"))) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::labs(x = "Normalization Method", y = "Comparison") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 2, barheight = 5, title = "Value"))
  # plot
  if(length(metrics) == 1){
   p <- p + ggplot2::scale_fill_distiller(name = metrics[1],palette = "YlOrRd", direction = 1, limits = c(0,1))
  } else {
    p <- p + ggplot2::facet_wrap(~Metric, ncol = 2) + ggplot2::scale_fill_distiller(name = "Measure",palette = "YlOrRd", direction = 1, limits = c(0,1))
  }
  return(p)
}


#' Barplot of true and false positives for specific comparisons and normalization methods
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object (barplot)
#' @export
#'
plot_TP_FP_spiked_bar <- function(stats, ain = NULL, comparisons = NULL){
  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- stats[, c("Assay", "Comparison", "TP", "FP")]
  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]

  dt$TP <- as.integer(dt$TP)
  dt$FP <- as.integer((-1) * as.integer(dt$FP))

  melted_dt <- data.table::melt(dt, variable.name = "Class", value.name = "Value", measure.vars = c("TP", "FP"))

  # plot
  p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x=get("Assay"), y=get("Value"), fill = get("Class"))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(y="False Positives | True Positives", x = "Normalization Method") +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "Class") +
    ggplot2::facet_wrap(~Comparison, scales="free_y") +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::coord_flip()
  return(p)
}

#' Boxplot of true and false positives for specific comparisons and normalization methods
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object (barplot)
#' @export
#'
plot_TP_FP_spiked_box <- function(stats, ain = NULL, comparisons = NULL){
  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- stats[, c("Assay", "Comparison", "TP", "FP")]
  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]

  dt$TP <- as.integer(dt$TP)
  dt$FP <- as.integer(dt$FP)

  melted_dt <- data.table::melt(dt, variable.name = "Class", value.name = "Value", measure.vars = c("TP", "FP"))

  # order methods
  tps <- melted_dt[melted_dt$Class == "TPs",]
  meds <- tps %>% dplyr::group_by(Assay) %>% dplyr::summarise(Median = median(Value, na.rm = TRUE)) %>% data.table::as.data.table()
  meds <- meds[order(meds$Median),]
  melted_dt_1 <- melted_dt
  melted_dt_1$Assay <- factor(melted_dt_1$Assay, levels = meds$Assay)

  # plot
  p <- ggplot2::ggplot(melted_dt_1, ggplot2::aes(x = get("Assay"), y = get("Value"), fill = get("Class"))) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "Class") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle =90, vjust =0.5)) +
    ggplot2::labs(x = "Normalization Method", y="Number of Proteins")
  return(p)
}
