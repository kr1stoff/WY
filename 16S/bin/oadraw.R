#!/usr/bin/env Rscript
options(warn = -1)
suppressPackageStartupMessages({
  library(funr)
  library(stringr)
  library(yaml)
  BIN = dirname(sys.script())
  f_software <- str_c(BIN, "../config/software.yml", sep = '/')
  software <- yaml.load_file(f_software)
  library(optparse)
  library(logging)
})
basicConfig()


option_list <- list(
  make_option(
    c("-t", "--type"),
    type = "character",
    default = "PCA",
    help = "defined which ordination to choose [default %default]"
  ),
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "otu table to plot [default %default]"
  ),
  make_option(
    c("-d", "--dist"),
    type = "character",
    default = NULL,
    help = "dist matrix to plot [default %default]"
  ),
  make_option(
    c("-c", "--constrain"),
    type = "character",
    default = NULL,
    help = "constrain table [default %default]"
  ),
  make_option(
    c("-f", "--family"),
    type = "character",
    default = NULL,
    help = "Set text font [default %default]"
  ),
  make_option(
    c("-s", "--size"),
    type = "numeric",
    default = "3",
    help = "size of points [default %default]"
  ),
  make_option(
    c("-a", "--alpha"),
    type = "numeric",
    default = "0.7",
    help = "alpha of points [default %default]"
  ),
  make_option(
    c("-g", "--groups"),
    type = "character",
    default = NULL,
    help = "treatment groups of samples, separated by tab, no spaces [default %default]"
  ),
  make_option(
    c("-e", "--ellipse"),
    action = "store_true",
    default = F,
    help = "draw ellipse for group points [default %default]"
  ),
  make_option(
    c("-x", "--scale"),
    action = "store_true",
    default = F,
    help = "Scale species to unit variance [default %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "output directory [default %default]"
  ),
  make_option(
    c("-l", "--loading"),
    type = "character",
    default = NULL,
    help = "whether or how to show species [default %default]"
  ),
  make_option(
    c("-p", "--plot"),
    type = "character",
    default = NULL,
    help = "general plot script [default %default]"
  ),
  make_option(
    c("-v", "--versus"),
    type = "character",
    default = NULL,
    help = "a&b or a&b&c [default %default]"
  )
)
opts <-
  parse_args(
    OptionParser(usage = "%prog [options]", option_list = option_list),
    positional_arguments = F
  )
logdebug(opts)

if (!dir.exists(opts$output))
  dir.create(opts$output, recursive = T)

log <- file(paste0(opts$output, "/analysis.log"))

datelog <- function(m) {
  sink(log, append = TRUE, split = TRUE)
  cat("[", date(), "] ", m, "\n", sep = "")
  sink()
}

error_test <- function(e, m) {
  if (e) {
    sink(
      paste(opts$output, "analysis.log", sep = "/"),
      append = TRUE,
      split = TRUE
    )
    cat("ERROR:\n\t", m, "\n\n")
    sink()
    datelog("Done with Error!")
    q()
  }
}

read_dist <- function(x) {
  tryCatch({
    loginfo("Reading distance matrix ...")
    dist.mt <-
      read.table(
        file = x,
        row = 1,
        sep = "\t",
        header = T,
        check.names = F
      )
    dist <- as.dist(dist.mt)
  },
    error = function(e) {
      error_test(T, conditionMessage(e))
    },
    finally = {
      return(dist)
    })
}

read_group <- function(x) {
  tryCatch({
    tb <-
      read.table(
        file = x,
        sep = "\t",
        header = T,
        check = F,
        comment = "",
        colClasses = "character",
        stringsAsFactors = F)
    tb <- unique(tb)
    rownames(tb) <- gsub('^ *| *$', '', tb[, 1])
    tb <- tb[, -1, drop = F]
    colnames(tb) <- sub('^ *| *$', '', colnames(tb))
    is_blank <- (colnames(tb) == "")
    if (sum(is_blank) > 0) {
      loginfo("Delete blank Column!")
      tb <- tb[, !is.blank]
    }
  },
    error = function(e) {
      error_test(T, conditionMessage(e))
    },
    finally = {
      return(tb)
    })
}


read_tb <- function(x) {
  tryCatch({
    tb <-
      read.table(
        file = x,
        sep = "\t",
        header = T,
        check = F,
        comment = "",
        stringsAsFactors = F
      )
    tb <- unique(tb)
    rownames(tb) <- gsub('^ *| *$', '', tb[, 1])
    tb <- tb[, -1, drop = F]
    colnames(tb) <- sub('^ *| *$', '', colnames(tb))
    is_blank <- (colnames(tb) == "")
    if (sum(is_blank) > 0) {
      loginfo("Delete blank Column!")
      tb <- tb[, !is.blank]
    }
  },
    error = function(e) {
      error_test(T, conditionMessage(e))
    },
    finally = {
      return(tb)
    })
}

loginfo("Pipe start now!")
oa <- list(type = opts$type)
knum <- 3

# read dist table
if (!is.null(opts$dist)) {
  oa$dist <- read_dist(opts$dist)
  loginfo(str(oa$dist))
}

# read otu table
if (!is.null(opts$input)) {
  loginfo("Reading otu table ...")
  oa$otu <- read_tb(opts$input)
  oa$otu <- oa$otu[rowSums(oa$otu) > 0,]
  loginfo(paste0("Samples num : ", ncol(oa$otu)))
  error_test(ncol(oa$otu) < 3, "Samples number of input file is less than 3!")
}

# read groups table
if (!is.null(opts$groups)) {
  loginfo("Reading groups file ...")
  oa$groups <- read_group(opts$groups)
  if (!is.null(opts$versus)) {
    comp <- unlist(strsplit(opts$versus, split = "&"))
    oa$groups <- oa$groups[oa$groups[, 1] %in% comp, , drop = F]
    loginfo(str(oa$groups))
  }
  tmp_gro <- oa$groups[, 1, drop = F]
  gro_file <- paste0(opts$output, '/group.txt')
  write.table(
    tmp_gro,
    file = gro_file,
    col.names = F,
    row.names = T,
    sep = "\t",
    quote = F
  )
}

# check samples
if (!is.null(oa$groups)) {
  groups_name <- row.names(oa$groups)
  loginfo("Check samples ...")
  if (!is.null(oa$otu)) {
    loginfo('Choose samples in table by groups!\n')
    oa$otu <- oa$otu[, row.names(oa$groups)]
    loginfo(str(oa$otu, list.len = 5))
  } else if (!is.null(oa$dist)) {
    loginfo('Choose samples in dist by groups!\n')
    oa$dist <- as.matrix(oa$dist)
    oa$dist <- oa$dist[groups_name, groups_name]
    oa$dist <- as.dist(oa$dist)
  }
}

# read constrain table
if (!is.null(opts$constrain)) {
  loginfo("Reading constrain file ...")
  oa$constrain <- read_tb(opts$constrain)
  loginfo(str(oa$constrain))
}

# do ordination analysis
if (oa$type == "PCA") {
  library(gmodels)
  error_test(is.null(oa$otu), "Input file was not defined!")
  loginfo(paste("Do ordination analysis of", oa$type, "..."))
  pca <-
    fast.prcomp(t(oa$otu),
                retx = T,
                scale = opts$scale,
                center = T)
  pca <- summary(pca)
  oa$sites <-
    as.data.frame(pca$x, check.names = F, stringsAsFactors = F)
  oa$eig <-
    data.frame(PC = colnames(pca$importance), t(pca$importance))
  oa$labs <- pca$importance[2, 1:3]
  oa$labs <-
    paste(names(oa$labs), sprintf("(%.2f%%)", oa$labs * 100), sep = "")

} else if (oa$type == "PCoA" | oa$type == "NMDS") {
  library(vegan, quietly = T)

  if (!is.null(oa$dist)) {
    #error_test(T, "Dist matrix was not defined!")
  } else if (!is.null(oa$otu)) {
    loginfo("Analysis distance of Bray-Curtis ...")
    oa$dist <- vegdist(t(oa$otu), method = "bray")
  } else {
    error_test(T, "Input file was not defined!")
  }

  loginfo(paste("Do ordination analysis of", oa$type, "..."))

  knum <- ifelse(attr(oa$dist, "Size") > 3, 3, 2)
  ifelse(
    oa$type == "PCoA",
    mds <-
      cmdscale(oa$dist, k = knum, eig = T),
    mds <- metaMDS(oa$dist, k = knum)
  )
  oa$sites <-
    as.data.frame(mds$points,
                  check.names = F,
                  stringsAsFactors = F)

  if (oa$type == "PCoA") {
    colnames(oa$sites) <- paste("PCO", 1:ncol(oa$sites), sep = "")
    ev <- summary(eigenvals(mds))
    oa$eig <- ev
    if (is.list(ev))
      oa$eig <- ev$importance
    axis <- paste0('PCO', 1:ncol(oa$eig))
    oa$eig <- data.frame(Axis = axis, t(oa$eig)[, -3])
    oa$labs <-
      sprintf("%s(%.2f%%)", colnames(oa$sites)[1:knum], oa$eig[1:knum, 3] * 100)
  } else {
    oa$labs <- colnames(oa$sites)[1:knum]
    oa$text <- sprintf("stress = %.3f", mds$stress)
    sink(paste0(opts$output, '/label.json'))
    cat(rjson::toJSON(list(stress = oa$text)))
    sink()
  }

} else if (oa$type == "RDA" | oa$type == "CCA") {
  library(vegan, quietly = T)

  ## check data
  error_test(is.null(oa$otu), "Input file was not defined!")
  error_test(is.null(oa$constrain), "Constrain file was not defined!")

  # get intersection of samples
  intercect <- intersect(row.names(oa$constrain), colnames(oa$otu))
  error_test(length(intercect) == 0, "No intercection!")
  cat("intercection (", length(intercect), ") :", intercect, "\n")
  oa$constrain <- oa$constrain[intercect, , drop = F]
  oa$otu <- oa$otu[, intercect]
  oa$otu <- oa$otu[rowSums(oa$otu) > 0,]
  oa$otu <- oa$otu[order(rowSums(oa$otu), decreasing = T),]

  # do DCA before CCA/RDA
  loginfo("Do ordination analysis of DCA ...")
  sink(
    paste(opts$output, "DCA.summary.txt", sep = "/"),
    append = F,
    split = TRUE
  )
  loginfo(decorana(oa$otu))
  sink()

  # do CCA/RDA
  loginfo(paste("Do ordination analysis of", oa$type, "..."))
  ifelse(
    oa$type == "RDA",
    cca <- rda(t(oa$otu) ~ ., oa$constrain, scale = opts$scale),
    cca <-
      cca(t(oa$otu) ~ ., oa$constrain, scale = opts$scale)
  )
  oa$mod <- cca
  # do anvoa for cca
  anova_out <- anova(cca, by = 'term')
  write.table(
    anova_out,
    paste0(opts$output, "/", 'anova.xls'),
    sep = "\t",
    col.names = NA,
    quote = F
  )

  # do envfit for cca
  fit <- envfit(cca, oa$constrain, perm = 999)
  envfit_out <-
    data.frame(
      fit$vectors$arrows,
      r2 = fit$vectors$r,
      'Pr(>r)' = fit$vectors$pvals,
      check.names = F
    )
  envfit_out <- round(envfit_out, 4)
  write.table(
    envfit_out,
    paste0(opts$output, "/", 'envfit_result.xls'),
    sep = "\t",
    col.names = NA,
    quote = F
  )

  # output cca result
  cca.sum <- summary(cca)
  # print(unclass(cca.sum)); q()
  oa$sites <-
    as.data.frame(cca.sum$sites,
                  check.names = F,
                  stringsAsFactors = F)
  oa$biplot <-
    as.data.frame(cca.sum$biplot,
                  check.names = F,
                  stringsAsFactors = F)
  oa$species <-
    as.data.frame(cca.sum$species,
                  check.names = F,
                  stringsAsFactors = F)
  if (nrow(oa$biplot) == 1)
    row.names(oa$biplot) <- colnames(oa$constrain)
  oa$eig <-
    data.frame(
      Axis = colnames(cca.sum$concont$importance),
      t(cca.sum$concont$importance),
      check.names = F
    )
  oa$labs <-
    sprintf("%s(%.2f%%)", colnames(oa$sites)[1:3], oa$eig[1:3, 3] * 100)
} else {
  error_test(T, "Ordination type was not true!")
}

## print head of sites
loginfo("sites :\n")
loginfo(str(oa$sites[, 1:knum]))

## get ellipse function
get_ell <- function(data) {
  library(ellipse)
  dat_ell <- data.frame()
  for (g in levels(data$Group)) {
    q1 <- data[data$Group == g, 1]
    q2 <- data[data$Group == g, 2]
    ell <-
      ellipse(cor(q1, q2),
              scale = c(sd(q1), sd(q2)),
              centre = c(mean(q1), mean(q2)))
    dat_ell <- rbind(dat_ell, data.frame(ell, Group = g))
  }
  return(dat_ell)
}

## pca_plot function
#library(plotly, quietly = T)
library(ggplot2)
library(ggrepel)

pca_plot <- function(oa, x, y) {
  #shape_value=c(15:17, 19, 18, 10:14, 9:0, 20:25)
  shape_value = c(15:17, 18, 10:14, 9:0, 20:25)
  samples <- row.names(oa$sites)
  sites <- oa$sites[, c(x, y)]
  colnames(sites) <- c("x", "y")

  if (!is.null(oa$groups)) {
    groups <- oa$groups
    samples <- row.names(groups)
    groups_name <- colnames(groups)
    if (ncol(groups) == 1) {
      colnames(groups) <- "Group"
      groups$Group <-
        factor(groups$Group, levels = unique(groups$Group))
    } else {
      colnames(groups) <- c("Group1", "Group2")
      groups$Group1 <-
        factor(groups$Group1, levels = unique(groups$Group1))
      groups$Group2 <-
        factor(groups$Group2, levels = unique(groups$Group2))
    }
    # Sample must be defined for plotly
    sites <-
      data.frame(
        sites[samples,],
        Sample = samples,
        groups,
        check.names = F,
        stringsAsFactors = F
      )
    color_num <- length(unique(groups[, ncol(groups)]))
  } else {
    color_num <- length(samples)
  }
  # print(str(sites))
  # print(samples)
  # cat("\n")

  pca.plot <- ggplot(data = sites, aes(x = x, y = y))
  if (is.null(oa$groups)) {
    pca.plot <- pca.plot +
      geom_point(alpha = opts$alpha,
                 size = opts$size,
                 aes(color = samples))
  } else if (ncol(groups) == 1) {
    pca.plot <- pca.plot +
      geom_point(alpha = opts$alpha,
                 size = opts$size,
                 #aes(shape = Group, color = Group, group = Sample)) +
                 aes(color = Group, group = Sample)) +
      #scale_shape_manual(values = shape_value) +
      labs(shape = groups_name, color = groups_name)
    if (opts$ellipse) {
      dat_ell <- get_ell(sites)
      #print(head(dat_ell));q()
      pca.plot <-
        pca.plot + geom_path(data = dat_ell, aes(x, y, color = Group))
      #pca.plot <- pca.plot + stat_ellipse(aes(color = Group))
    }

  } else {
    pca.plot <- pca.plot +
      geom_point(
        alpha = opts$alpha,
        size = opts$size,
        aes(
          shape = Group1,
          color = Group2,
          group = Sample
        )
      ) +
      scale_shape_manual(values = shape_value) +
      labs(shape = groups_name[1], color = groups_name[2])
  }

  pca.plot <-
    pca.plot +
      labs(x = oa$labs[x],
           y = oa$labs[y],
           title = oa$type) +
      geom_hline(yintercept = 0,
                 linetype = 4,
                 color = "grey") +
      geom_vline(xintercept = 0,
                 linetype = 4,
                 color = "grey") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 10))

  if (!is.null(oa$text))
    pca.plot <- pca.plot + annotate(
      "text",
      label = oa$text,
      size = 5,
      x = max(sites$x) * .85,
      y = max(sites$y)
    )

  #if(color_num < 10) pca.plot <- pca.plot + scale_color_brewer(palette = "Set1")
  set1_col <- RColorBrewer::brewer.pal(9, 'Set1')[c(1:5, 7:9)]
  if (color_num < 9)
    pca.plot <- pca.plot + scale_color_manual(values = set1_col)

  if (!is.null(oa$biplot)) {
    bip <- scores(oa$mod, c(x, y), "bp", "species")
    mul <- ordiArrowMul(bip, 2)
    if (mul == Inf) {
      mul <- 2
    }
    cat('mul:', mul, '\n')
    bip.scl <- bip * mul

    plot(oa$mod)
    bip.lab <- ordiArrowTextXY(bip.scl, rescale = FALSE)

    pca.plot <- pca.plot +
      geom_segment(
        data = data.frame(bip.scl),
        aes(
          x = 0,
          y = 0,
          xend = bip.scl[, 1],
          yend = bip.scl[, 2]
        ),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "red"
      ) +
      geom_text(
        data = data.frame(bip.lab),
        aes(
          label = row.names(bip.lab),
          x = bip.lab[, 1],
          y = bip.lab[, 2]
        ),
        color = "red",
        parse = T
      )

    # show species
    if (!is.null(opts$loading)) {
      num <- nrow(oa$species)
      if (opts$loading == 'top10') {
        num <- ifelse(nrow(oa$species) < 10, nrow(oa$species), 10)
      }
      sp_show <- oa$species[1:num,]
      pca.plot <-
        pca.plot + geom_text(
          data = sp_show,
          color = "blue",
          aes(
            label = row.names(sp_show),
            x = sp_show[, x],
            y = sp_show[, y]
          )
        )
    }

  }

  if (!is.null(opts$family)) {
    cat('family:', opts$family, '\n')
    pca.plot <-
      pca.plot + theme(text = element_text(family = opts$family))
    pca.text <-
      pca.plot + geom_text_repel(aes(label = samples),
                                 size = opts$size,
                                 family = opts$family)
  } else {
    #pca.text <- pca.plot + geom_text(aes(label = samples), size = opts$size, vjust = -1)
    pca.text <-
      pca.plot + geom_text_repel(aes(label = samples), size = opts$size)
  }

  out_pre <- paste0(opts$output, "/", oa$type, x, "-", y)
  ggsave(paste(out_pre, ".png", sep = ""),
         pca.text,
         width = 7,
         height = 6)
  ggsave(paste(out_pre, ".pdf", sep = ""),
         pca.text,
         width = 7,
         height = 6)
  ggsave(
    paste(out_pre, ".nonlabs.png", sep = ""),
    pca.plot,
    width = 7,
    height = 6
  )
  ggsave(
    paste(out_pre, ".nonlabs.pdf", sep = ""),
    pca.plot,
    width = 7,
    height = 6
  )
  #ply <- ggplotly(pca.plot)
  #htmlwidgets::saveWidget(ply, file = paste(out_pre, ".html", sep = ""), selfcontained = F,
  #                        libdir = paste(opts$output, "/lib", sep = ""))
}


## output statistic
sites_file <-
  paste(opts$output, "/", oa$type, "_sites.xls", sep = "")
sites <-
  data.frame(
    Samples = row.names(oa$sites),
    oa$sites,
    stringsAsFactors = F,
    check.names = F
  )
write.table(
  sites,
  file = sites_file,
  col = T,
  row = F,
  sep = "\t",
  quote = F
)

if (!is.null(oa$eig)) {
  eig_file <- paste(opts$output, "/", oa$type, "_eig.xls", sep = "")
  write.table(
    oa$eig,
    file = eig_file,
    col.names = T,
    row.names = F,
    sep = "\t",
    quote = F
  )
}
if (!is.null(oa$biplot)) {
  file <- paste0(opts$output, "/", oa$type)
  fits <-
    data.frame(
      Constrain = row.names(oa$biplot),
      oa$biplot,
      stringsAsFactors = F,
      check.names = F
    )
  write.table(
    fits,
    file = paste0(file, "_biplot.xls"),
    row = F,
    sep = "\t",
    quote = F
  )
  sps <-
    data.frame(
      Species = row.names(oa$species),
      oa$species,
      stringsAsFactors = F,
      check.names = F
    )
  write.table(
    sps,
    file = paste0(file, "_species.xls"),
    row = F,
    sep = "\t",
    quote = F
  )
}

## plot
if (is.null(opts$plot)) {
  logdebug("Plot ...")
  pca_plot(oa, 1, 2)
} else {
  tmp_file <-
    paste(opts$output, "/", oa$type, ".data4fig.xls", sep = "")
  tmp_data <- sites
  if (ncol(sites) > 4) {
    tmp_data <- sites[1:4]
  }
  colnames(tmp_data)[-1] <- oa$labs
  write.table(
    tmp_data,
    file = tmp_file,
    col.names = T,
    row.names = F,
    sep = "\t",
    quote = F
  )
  tmp_group <-
    ifelse(!is.null(opts$groups), paste('-group', gro_file), '')
  outprefix <- paste0(opts$output, "/", oa$type, '1-2')
  title <-
    ifelse(!is.null(oa$text),
           paste0(oa$type, '(', oa$text, ')'),
           oa$type)
  title <- sprintf('\'%s\'', title)
  cmd <-
    paste(
      'perl',
      opts$plot,
      'pca -skipstat T -file',
      tmp_file,
      tmp_group,
      '-title',
      title,
      '-outprefix',
      outprefix,
      '-label T'
    )
  cat('[CMD]', cmd, '\n')
  system(cmd)
  outprefix <- paste0(outprefix, '.nonlabs')
  cmd <-
    paste(
      'perl',
      opts$plot,
      'pca -skipstat T -file',
      tmp_file,
      tmp_group,
      '-title',
      title,
      '-outprefix',
      outprefix,
      '-label F'
    )
  system(cmd)
}
loginfo("All Done!")
