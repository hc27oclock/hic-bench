#!/usr/bin/env Rscript
usage = "\
Rscript optimal_lambda.r OUTPUT-DIR INPUT-COMPARE-TSV MIN-IMPROVEMENT
"

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsignif))

stat_compare_means <- function(mapping = NULL, data = NULL,
                               method = NULL, paired = FALSE, method.args = list(exact=FALSE), ref.group = NULL,
                               comparisons = NULL, hide.ns = FALSE, label.sep = ", ",
                               label = NULL, label.x.npc = "left", label.y.npc = "top",
                               label.x = NULL, label.y = NULL, tip.length = 0.03,
                               symnum.args = list(),
                               geom = "text", position = "identity",  na.rm = FALSE, show.legend = NA,
                               inherit.aes = TRUE, ...) {
  .add_item <- function(.list, ...){
    pms <- list(...)
    for(pms.names in names(pms)){
      .list[[pms.names]] <- pms[[pms.names]]
    }
    .list
  }
  .method_info <- function(method){
    
    if(is.null(method))
      method = "wilcox.test"
    
    allowed.methods <- list(
      t = "t.test", t.test = "t.test", student = "t.test",
      wiloxon = "wilcox.test", wilcox = "wilcox.test", wilcox.test = "wilcox.test",
      anova = "anova", aov = "anova",
      kruskal = "kruskal.test", kruskal.test = "kruskal.test")
    
    method.names <- list(
      t.test = "T-test", wilcox.test = "Wilcoxon",
      anova = "Anova", kruskal.test = "Kruskal-Wallis")
    
    if(!(method %in% names(allowed.methods)))
      stop("Non-supported method specified. Allowed methods are one of: ",
           .collapse(allowed.methods, sep =", "))
    method <- allowed.methods[[method]]
    method.name <- method.names[[method]]
    
    list(method = method, name = method.name)
  }
  if(!is.null(comparisons)){
    
    method.info <- .method_info(method)
    method <- method.info$method
    method.args <- .add_item(method.args, paired = paired)
    if(method == "wilcox.test")
    pms <- list(...)
    size <- ifelse(is.null(pms$size), 0.2, pms$size)
    color <- ifelse(is.null(pms$color), "black", pms$color)
    
    map_signif_level <- FALSE
    if(is.null(label)) label <- "p.format"
    
    if(.is_p.signif_in_mapping(mapping) | (label %in% "p.signif"))
    {
      map_signif_level <- c("****"=0.0001, "***"=0.001, "**"=0.01,  "*"=0.05, "ns"=1)
      if(hide.ns) names(map_signif_level)[5] <- " "
    }
    
    step_increase <- ifelse(is.null(label.y), 0.12, 0)
    ggsignif::geom_signif(comparisons = comparisons, y_position = label.y,
                          test = method, test.args = method.args,
                          step_increase = step_increase, size = size, textsize = 2.5, color = color,
                          map_signif_level = map_signif_level, tip_length = tip.length, vjust=0.15)
  }
  
  else{
    mapping <- .update_mapping(mapping, label)
    layer(
      stat = StatCompareMeans, data = data, mapping = mapping, geom = geom,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(label.x.npc  = label.x.npc , label.y.npc  = label.y.npc,
                    label.x = label.x, label.y = label.y, label.sep = label.sep,
                    method = method, method.args = method.args,
                    paired = paired, ref.group = ref.group,
                    symnum.args = symnum.args,
                    hide.ns = hide.ns, na.rm = na.rm, ...)
    )
    
  }
}


StatCompareMeans<- ggproto("StatCompareMeans", Stat,
                           required_aes = c("x", "y"),
                           default_aes = aes(hjust = ..hjust.., vjust = ..vjust..),
                           
                           compute_panel = function(data, scales, method, method.args,
                                                    paired, ref.group, symnum.args,
                                                    hide.ns, label.x.npc, label.y.npc,
                                                    label.x, label.y, label.sep)
                           {
                             . <- x <- NULL
                             .is.multiple.grouping.vars <- !all(data$x == data$group)
                             
                             if(!is.null(ref.group)) {
                               if(ref.group != ".all.") ref.group <- scales$x$map(ref.group)
                             }
                             
                             # Guess the number of group to be compared
                             #::::::::::::::::::::::::::::::::::::::::::::::::::
                             if(.is.multiple.grouping.vars)
                               x.levels <- .levels(data$group)
                             else x.levels <- .levels(data$x)
                             two.groups <- length(x.levels) == 2 | !is.null(ref.group)
                             multi.groups <- length(x.levels) > 2
                             
                             # Guess the test to be performed
                             #::::::::::::::::::::::::::::::::::::::::::::::::::
                             if(two.groups & is.null(method))
                               method <- "wilcox.test"
                             else if(multi.groups & is.null(method))
                               method <- "kruskal.test"
                             
                             # Perform group comparisons
                             #::::::::::::::::::::::::::::::::::::::::::::::::::
                             method.args <- method.args %>%
                               .add_item(data = data, method = method,
                                         paired = paired, ref.group = ref.group,
                                         symnum.args = symnum.args)
                             
                             if(.is.multiple.grouping.vars){
                               method.args <- method.args %>%
                                 .add_item(formula = y ~ group, group.by = "x")
                               .test <- do.call(compare_means, method.args)
                             }
                             else{
                               method.args <- method.args %>%
                                 .add_item(formula = y ~ x)
                               .test <- do.call(compare_means, method.args)
                             }
                             
                             pvaltxt <- ifelse(.test$p < 2.2e-16, "p < 2.2e-16",
                                               paste("p =", signif(.test$p, 2)))
                             .test$label <- paste(.test$method, pvaltxt, sep =  label.sep)
                             
                             # Options for label positioning
                             #::::::::::::::::::::::::::::::::::::::::::::::::::
                             label.opts <- list(data = data, scales = scales,
                                                label.x.npc = label.x.npc, label.y.npc = label.y.npc,
                                                label.x = label.x, label.y = label.y,
                                                symnum.args = symnum.args, .by = "panel" )
                             
                             if(.is.multiple.grouping.vars){
                               
                               if(is.null(label.x) & length(label.x.npc) == 1)
                                 label.opts$label.x <- .test$x
                               
                               .label.pms <- label.opts %>%
                                 .add_item(group.ids = .test$x) %>%
                                 do.call(.label_params_by_group, .) # Returns a data frame with label: x, y, hjust, vjust
                               # .test <- dplyr::select(.test, -x)
                               .label.pms <- dplyr::select(.label.pms, -x)
                               
                             }
                             
                             else{
                               .label.pms <- label.opts %>%
                                 do.call(.label_params, .) %>% # Returns a data frame with label: x, y, hjust, vjust
                                 dplyr::mutate(hjust = 0.2)
                             }
                             if(!is.null(ref.group)){
                               group.ids <- as.numeric(.test$group2)
                               if(!is.null(label.y) & ref.group != ".all."){
                                 if(length(label.y) == length(group.ids))
                                   label.opts$label.y <- c(0, label.y)
                               }
                               .label.pms <- label.opts %>%
                                 .add_item(group.ids = group.ids) %>%
                                 do.call(.label_params_by_group, .)
                             }
                             
                             res <- cbind(.test, .label.pms)
                             if(!is.null(ref.group)){
                               # Set label x value to group names
                               other.group.index <- as.numeric(res$group2)
                               res$x <- scales$x$range$range[other.group.index ]
                               res <- res %>% dplyr::mutate(hjust = 0.5)
                             }
                             
                             if(hide.ns){
                               p.signif <- res$p.signif
                               p.format <- res$p.format
                               p.signif[p.signif == "ns"] <- " "
                               res$p.signif <- p.signif
                             }
                             res
                           }
                           
)


# Check if p.signif is in mapping
.is_p.signif_in_mapping <- function(mapping){
  
  res <- FALSE
  if(!is.null(mapping)){
    if(!is.null(mapping$label)){
      .label <- as.character(mapping$label)
      res <- "..p.signif.." %in% .label
    }
  }
  return(res)
}

# Update mapping with label
.update_mapping <- function (mapping, label){
  
  allowed.label <- list(
    "p.signif" = quote(..p.signif..),
    "..p.signif.." = quote(..p.signif..),
    "p.format" = quote(paste0("p = ",..p.format..)),
    "..p.format.." = quote(paste0("p = ",..p.format..)),
    "p" = quote(paste0("p = ",..p.format..)),
    "..p.." = quote(paste0("p = ",..p.format..))
  )
  
  if(!is.null(label)){
    if(!label %in% names(allowed.label) )
      stop("Allowed values for label are: ", .collapse(names(allowed.label) , sep = ", "))
  }
  if(!is.null(mapping) & is.character(label)){
    mapping$label <- allowed.label[[label]]
  }
  else if(is.character(label)){
    mapping <- aes()
    mapping$label <- allowed.label[[label]]
  }
  mapping
}

args <- commandArgs(trailingOnly = TRUE)
outdir=args[1]
inputs=args[2]
min_improvement=as.numeric(args[3])
paired=FALSE

compare=read.table(inputs,header=T)
for (k in compare$LAMBDA %>% sort %>% unique %>% -1){
  if(k!=0){
    p=wilcox.test((1.0+min_improvement)*compare[compare$LAMBDA==k,]$VALUE,compare[compare$LAMBDA==k+1,]$VALUE,alternative = 'less',paired=paired)$p.value
    cat(paste(unique(compare$SAMPLE.1),unique(compare$SAMPLE.2),k,p,sep="\t"))
    cat("\n")
  }
}

outplot = ggplot(compare,aes(x=as.factor(LAMBDA),y=VALUE,fill=as.factor(LAMBDA))) + geom_boxplot() + xlab(expression(lambda)) + ylab("Stratum adjusted correlation coefficient (SCC)") + labs(fill='Lambda') + scale_fill_manual(values=c("#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5")) + stat_compare_means(comparisons =list( c("1", "2"), c("2", "3"), c("3", "4"), c("4", "5"), c("5","6")), method = "wilcox.test", paired = F, method.args = list(alternative="less", exact=NULL)) + theme_bw() + theme(axis.text=element_text(size=8),axis.title=element_text(size=8),strip.text = element_text(size = 8), legend.position="none") + scale_x_discrete(labels=c("1" = "0", "2" = "0.2", "3" = "0.4", "4" = "0.6", "5" = "0.8", "6" = "1"))

ggsave(paste(outdir,"boxplots.pdf",sep='/'),plot=outplot, width=2.75, height=3.75)

