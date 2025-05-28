args <- commandArgs(trailingOnly = FALSE)
Bin <- dirname(sub("--file=", "", args[grep("--file=", args)]))
source(paste(Bin, "Colors.R", sep = "/"), chdir = T)

args  <- commandArgs(TRUE)
color_type <- args[1]
color_tag  <- args[2]
color_num  <- args[3]

color_num <- as.numeric(color_num)
color.return <- fetch_color(n = color_num, type = color_type, tag = color_tag)
cat(paste(unlist(color.return), collapse=","))

