

#' Volcano plots of DE results
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param facet_norm Boolean indicating whether to facet by normalization method (TRUE) or not (FALSE)
#' @param facet_comparison Boolean indicating whether to facet by comparison (TRUE) or not (FALSE). Only valid if facet_norm = FALSE.
#'
#' @return list of ggplot objects
#' @export
#'
plot_volcano_DE <- function(de_res, ain = NULL, comparisons = NULL, facet_norm = TRUE, facet_comparison = FALSE){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]
  de_res <- de_res[de_res$Assay %in% ain,]

  # plot
  p <- list()
  if(facet_norm){
    # facet by normalization method
    for(comp in comparisons){
      dt <- de_res[de_res$Comparison == comp,]
      if("Significant Change" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
      } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
      } else {
        color_values <- c("No Change" = "grey")
      }
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("logFC"), y=-log10(get("P.Value")))) +
        ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(ggplot2::aes(col=get("Change"))) +
        ggplot2::labs(y = expression(-log[10]~"P-value"), x= "logFC") +
        ggplot2::scale_color_manual(values = color_values, name = "Change") +
        ggplot2::facet_wrap(~Assay, scales = "free")
      p[[comp]] <- tmp
    }
  } else if (facet_comparison){
    # facet by comparison
    for(method in ain){
      dt <- de_res[de_res$Assay == ain,]
      if("Significant Change" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
      } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
      } else {
        color_values <- c("No Change" = "grey")
      }
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("logFC"), y=-log10(get("P.Value")))) +
        ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(ggplot2::aes(col=get("Change"))) +
        ggplot2::labs(y = expression(-log[10]~"P-value"), x= "logFC") +
        ggplot2::scale_color_manual(values = color_values, name = "Change") +
        ggplot2::facet_wrap(~Comparison, scales = "free")
      p[[method]] <- tmp
    }
  } else {
    # no facet at all --> individual plots for each method and each comparison
    for(comp in comparisons){
      for(method in ain){
        dt <- de_res[de_res$Assay == ain,]
        dt <- dt[dt$Comparison == comp,]
        if("Significant Change" %in% dt$Change){
          color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
        } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
          color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
        } else {
          color_values <- c("No Change" = "grey")
        }
        tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("logFC"), y=-log10(get("P.Value")))) +
          ggplot2::geom_vline(xintercept=0) +
          ggplot2::geom_point(ggplot2::aes(col=get("Change"))) +
          ggplot2::labs(y = expression(-log[10]~"P-value"), x= "logFC") +
          ggplot2::scale_color_manual(values = color_values, name = "Change")
        p[[paste0(comp, "_", method)]] <- tmp
      }
    }
  }
  return(p)
}



#TODO: sort bars by increasing number of DE proteins with a parameter

#' Overview plots of DE results
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#'
#' @return list of ggplot objects
#' @export
#'
plot_overview_DE_bar <- function(de_res, ain = NULL, comparisons = NULL){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  # get overview DE
  dt <- get_overview_DE(de_res)
  dt <- dt[dt$Assay %in% ain, ]
  melted_dt <- data.table::melt(dt, id.vars = c("Assay", "Comparison"), variable.name = "Change", value.name = "N")
  melted_dt$Assay <- factor(melted_dt$Assay, levels = sort(as.character(unique(melted_dt$Assay))))
  # plot
  p <- list()
  for(comp in comparisons){
    dt <- melted_dt[melted_dt$Comparison == comp,]
    if("Significant Change" %in% dt$Change){
      color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
    } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
      color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
    } else {
      color_values <- c("No Change" = "grey")
    }
    tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("N"), y = get("Assay"), fill = get("Change"), label = get("N"))) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_stack()) +
      ggplot2:: scale_fill_manual(values = color_values, name = "Change") +
      ggplot2::labs( y = "Normalization Method", x = "Number of DE Proteins")
    p[[comp]] <- tmp
  }
  return(p)
}

#' Upset plots of DE results of the different normalization methods
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param min_degree Minimal degree of an intersection for it to be included
#'
#' @return list of list of ComplexUpset objects and list of data tables with intersections (one for each comparison)
#' @export
#'
plot_upset_DE <- function(de_res, ain = NULL, comparisons = NULL, min_degree = 2) {
    # check parameters
    tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
    de_res <- tmp[[1]]
    ain <- tmp[[2]]
    comparisons <- tmp[[3]]
    dt <- de_res[de_res$Assay %in% ain, ]
    p <- list()
    t <- list()
    for (comp in comparisons) {
      dt <- dt[dt$Comparison == comp, ]
      # only extract significant changes
      dt <- dt[dt$Change %in% c("Up Regulated", "Down Regulated", "Significant Change"), ]
      if (nrow(dt) == 0) {
        warning(paste0("No significant changes for comparison ", comp, ": nothing to plot."))
      } else {
        # prepare data for ComplexUpset
        dt <- dt[, c("Assay", "Protein.IDs"), with = FALSE]
        dt <- unique(dt)
        dt <- data.table::dcast(dt, Assay ~ Protein.IDs)
        dt[is.na(dt)] <- 0
        assays <- dt$Assay
        dt$Assay <- NULL
        dt[dt != 0] <- 1
        dt <- as.data.frame(dt)
        row.names(dt) <- assays
        dt <- t(dt)
        dt <- dt == 1
        dt <- dt[names(sort(rowSums(dt))), ]
        dt <- as.data.frame(dt)
        # plot
        upset <- ComplexUpset::upset(
          dt,
          colnames(dt),
          name = "",
          set_sizes = ComplexUpset::upset_set_size(position = "right") + ggplot2::ylab("Set Size"),
          sort_sets = FALSE,
          keep_empty_groups = FALSE,
          sort_intersections = "descending",
          min_degree = min_degree,
          base_annotations = list("Intersection Size" =
                                    ComplexUpset::intersection_size(text = list(size = 3))),
          themes = ComplexUpset::upset_default_themes(text = ggplot2::element_text(size = 12))
        )
        p[[comp]] <- upset

        # prepare data table of intersections
        nr_methods <- rowSums(dt)
        t <- purrr::map2_df(dt, names(dt), ~  replace(.x, .x==TRUE, .y))
        t[t == FALSE] <- NA
        t <- as.data.frame(t)
        rownames(t) <- rownames(dt)
        res <- t %>% tidyr::unite(., col = "Assays", na.rm=TRUE, sep = ",")
        res$Nr <- nr_methods
        res <- res[order(-res$Nr),]
        res$Protein.IDs <- rownames(res)
        res <- res[, c("Protein.IDs", "Nr", "Assay")]
        colnames(res) <- c("Protein.IDs", "Number of Intersected Assays", "Assays")
        t[[comp]] <- res
      }
    }
    return(list("plots" = p, "tables" = t))
}

