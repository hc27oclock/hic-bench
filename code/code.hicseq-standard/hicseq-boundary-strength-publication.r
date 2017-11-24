#!/usr/bin/env Rscript
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('tidyr'))
suppressPackageStartupMessages(library('stringr'))
suppressPackageStartupMessages(library('lubridate'))
suppressPackageStartupMessages(library('readr'))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('cowplot'))
suppressPackageStartupMessages(library('pheatmap'))
suppressPackageStartupMessages(library('tools'))
suppressPackageStartupMessages(library('ggpubr'))
suppressPackageStartupMessages(library('RColorBrewer'))

args <- commandArgs(trailingOnly = TRUE)


stat_compare_means <- function(mapping = NULL, data = NULL,
                               method = NULL, paired = FALSE, method.args = list(), ref.group = NULL,
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
      method.args$exact <- FALSE
    
    
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
                          step_increase = step_increase, size = size, textsize = 1.5, color = color,
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


for (file in strsplit(args[1],",") %>% unlist()){
ctcf=read_delim(file,delim = '\t',col_names = F)
colnames(ctcf) <- c("sample","Boundary Strength Category","Normalized Aggregate CTCF Level")
#ctcf = ctcf %>% mutate(`Normalized Aggregate CTCF Level`=as.numeric(`Normalized Aggregate CTCF Level`)) %>% mutate(`Boundary Strength Category`=gsub("b","",`Boundary Strength Category`))
ctcf = ctcf %>% mutate(`Normalized Aggregate CTCF Level`=as.numeric(`Normalized Aggregate CTCF Level`)) %>% mutate(`Boundary Strength Category`=gsub("b1","I",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b2","II",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b3","III",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b4","IV",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("b5","V",`Boundary Strength Category`))
figure=ctcf %>% 
  ggplot(aes(x=`Boundary Strength Category`, y=`Normalized Aggregate CTCF Level`, fill=`Boundary Strength Category`)) + 
  scale_fill_manual(values=brewer.pal(5, "Greens")) +
  geom_boxplot(width = 0.5,lwd=0.3) +
  geom_jitter(size=0.05) + 
  theme_bw() +
  theme(legend.position="none",axis.text=element_text(size=8),axis.title = element_blank(),plot.margin=grid::unit(c(0.01,0.01,0.01,0.01), "inch"), axis.ticks.length=unit(0.02, "inch")) +
  stat_compare_means(comparisons =list( c("V", "IV"), c("V", "III"), c("V", "II"), c("V", "I")),method = "wilcox.test", paired = T)
ggsave(paste0(file_path_sans_ext(file),".pdf"), plot = last_plot(), width = 2.25, height = 1.25)
figure_ctcf=ctcf %>% 
  ggplot(aes(x=`Boundary Strength Category`, y=`Normalized Aggregate CTCF Level`, fill=`Boundary Strength Category`)) + 
  scale_fill_manual(values=brewer.pal(5, "Greens")) +
  geom_boxplot(width = 0.5,lwd=0.3) +
  geom_jitter(size=0.05) + 
  theme_bw() +
  theme(legend.position="none",axis.text=element_text(size=8),axis.title = element_text(size=8),plot.margin=grid::unit(c(0.01,0.01,0.01,0.01), "inch"), axis.ticks.length=unit(0.02, "inch")) +
  stat_compare_means(comparisons =list( c("V", "IV"), c("V", "III"), c("V", "II"), c("V", "I")),method = "wilcox.test", paired = T)
  cwd=getwd()
save.image(file = paste0(file_path_sans_ext(file),".RData"))
}

se=read_delim(args[2],delim = '\t',col_names = F)
colnames(se) <- c("sample","Boundary Strength Category","Fraction of Super-enhancers")
#se = se %>% mutate(`Fraction of Super-enhancers`=as.numeric(`Fraction of Super-enhancers`)) %>% mutate(`Boundary Strength Category`=gsub("0","",`Boundary Strength Category`))
se = se %>% mutate(`Fraction of Super-enhancers`=as.numeric(`Fraction of Super-enhancers`)) %>% mutate(`Boundary Strength Category`=gsub("01","I",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("02","II",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("03","III",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("04","IV",`Boundary Strength Category`)) %>% mutate(`Boundary Strength Category`=gsub("05","V",`Boundary Strength Category`))
figure=se %>% 
  ggplot(aes(x=`Boundary Strength Category`,y=`Fraction of Super-enhancers`,fill=`Boundary Strength Category`)) + 
  scale_fill_manual(values=brewer.pal(5, "Purples")) +
  geom_boxplot(width = 0.5,lwd=0.3) +
  geom_jitter(size=0.05) + 
  theme_bw() +
  theme(legend.position="none",axis.text=element_text(size=8),axis.title = element_blank(),plot.margin=grid::unit(c(0.01,0.01,0.01,0.01), "inch"), axis.ticks.length=unit(0.02, "inch")) +
  stat_compare_means(comparisons =list( c("V", "IV"), c("V", "III"), c("V", "II"), c("V", "I")),method = "wilcox.test", paired = T)
ggsave(paste0(file_path_sans_ext(args[2]),".pdf"), plot = last_plot(), width = 2.25, height = 1.25)
save.image(file = paste0(file_path_sans_ext(args[2]),".RData"))
figure_se=se %>% 
  ggplot(aes(x=`Boundary Strength Category`,y=`Fraction of Super-enhancers`,fill=`Boundary Strength Category`)) + 
  scale_fill_manual(values=brewer.pal(5, "Purples")) +
  geom_boxplot(width = 0.5,lwd=0.3) +
  geom_jitter(size=0.05) + 
  theme_bw() +
  theme(legend.position="none",axis.text=element_text(size=8),axis.title = element_text(size=8),plot.margin=grid::unit(c(0.01,0.01,0.01,0.01), "inch"), axis.ticks.length=unit(0.02, "inch")) +
  stat_compare_means(comparisons =list( c("V", "IV"), c("V", "III"), c("V", "II"), c("V", "I")),method = "wilcox.test", paired = T)
  cwd=getwd()


df=read_delim(args[3],delim = '\t',col_names = T)
df = df %>% as.data.frame()
rownames(df)=df[,1]
df=df[,-1]
kmean_result=kmeans(df,10)
hm=pheatmap(df[order(kmean_result$cluster),],breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5), color=rev(brewer.pal(6, "RdYlBu")),legend_breaks=0:5,legend_labels=c("NA","I","II","III","IV","V"),cluster_rows = F,show_rownames = F,fontsize=8, width = 4.2935*2,height = 3.4705*2, treeheight_col=15,cellwidth=17)
legend_title=ggdraw() + draw_label("Boundary Strength Category", fontface='bold', angle=270, size = 8, hjust = 1, vjust = 0) + theme(plot.margin=grid::unit(c(0.01,0.01,0.01,0.01), "inch"))
myplot=plot_grid(hm$gtable,legend_title,nrow=1,rel_widths=c(24,1))
ggsave(paste0(file_path_sans_ext(args[3]),".pdf"),plot=myplot,height=4,width=6.25)
save.image(file = paste0(file_path_sans_ext(args[3]),".RData"))
