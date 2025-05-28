
#color.list <- yaml::yaml.load(color.yaml)

#args <- commandArgs(trailingOnly = FALSE)
#Bin <- dirname(sub("--file=", "", args[grep("--file=", args)]))
#color.yaml <- paste(Bin, "Colors.yml", sep = "/")

## if you cannot source this script correctly in nested scripts, 
## you could try like this:
## source("Colors.R", chdir = T)
#script.dir <- dirname(sys.frame(1)$ofile)
#color.yaml <- paste(script.dir, "Colors.yml", sep = "/")
#color.list <- yaml::yaml.load_file(color.yaml)

## you should source this script with the argu [chdir = T]
color.yaml <- "Colors.yml"
color.list <- yaml::yaml.load_file(color.yaml)

fetch_color <- function(n = 0, type = c(names(color.list), "random"), tag = NULL,
  is.extend = TRUE, verbose = FALSE) {
  
  type <- match.arg(type)

  if ( type == "random" ) {
    color.use <- fetch_random_color(n = n, usepalette = T)
    if ( verbose ) {
      message("color : ", type, "->", tag, "->", n)
	}
	return(color.use)
  }

  tag  <- match.arg(tag, names(color.list[[type]]))
  if (!is.numeric(n))
	stop("'n' must be numeric.")
  
  n.available <- as.numeric(names(color.list[[type]][[tag]]))
  n.select <- n.available[n.available >= n]
  if (length(n.select)){
	if(n == min(n.select)) {
		n.select <- as.character(min(n.select))
		color.use <- color.list[[type]][[tag]][[n.select]]
	} else {
		n.select <- as.character(min(n.select))
		color.use <- head(color.list[[type]][[tag]][[n.select]], n)
	}
  } else {
    n.select <- as.character(0)
    color.use <- color.list[[type]][[tag]][[n.select]]
    if (is.extend)
      color.use <- colorRampPalette(color.list[[type]][[tag]][[n.select]])(n)
  }
  if ( verbose ) {
    message("color : ", type, "->", tag, "->", n.select)
  }
  
  return(color.use)
}

fetch_random_color <- function(n = 1, usepalette = FALSE, hue = " ", luminosity = " ") {
  library(randomcoloR)
  if(usepalette == TRUE) {
    set.seed(1)
    color.use <- distinctColorPalette(k = n)
  } else {
    color.use <- randomColor(count = n, hue = hue, luminosity = luminosity)
  }
  return(color.use)
}

show_all_color <- function() {
  library(dplyr, warn.conflicts = F)
  dt <- reshape2::melt(color.list)
  dt <- dt %>% group_by(L1,L2,L3) %>% mutate(x = 1:n()) %>% arrange(L1, L2, as.numeric(L3), x)
  dt$L3 <- factor(dt$L3, levels = as.character(0:max(as.numeric(dt$L3))))
  
  cc <- unique(as.character(dt$value))
  names(cc) <- cc
  
  library(ggplot2, warn.conflicts = F)
  p <- list()
  for ( i in unique(dt$L1) ) {
    p[[i]] <-  ggplot(dt %>% filter(L1 == i), aes(x = x, y = L3)) + geom_tile(aes(fill = value), color = "grey") +
      facet_grid(L2 ~ ., scales = "free_y", switch = "y") +
      scale_fill_manual(values = cc) +
      scale_x_continuous(expand = expand_scale()) +
	  geom_text(aes(label = value)) + 
      ylab(i) + 
      theme_minimal() + 
      theme(legend.position = "none", strip.placement = "outside", panel.grid = element_blank(),
            axis.title.x = element_blank(), axis.text.x = element_blank())
  }

  ap <- cowplot::plot_grid(plotlist = p, ncol = 1, axis = "ltrb", align = "hv")
  cowplot::save_plot("Color.test.pdf", ap, base_width = 50, base_height = 100, limitsize = FALSE)

}

