#-----------------------------------------------------------------------------
## source this script with chdir = T : source("all.theme.r", chdir = T)
.my_now_dir <- getwd()
ttf_dir     <- paste0(.my_now_dir, "/../../fonts/tff/msttcore")
arial_ttf   <- paste0(.my_now_dir, "/../../fonts/tff/msttcore/arial.ttf")
ggsave_R    <- paste0(.my_now_dir, "/ggsave.R")

line_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	options(scipen = -1)
	mytheme <- theme_bw() + 
	theme(
	    panel.grid = element_blank(),
		panel.border = element_rect(color = "#000000", size = 0.8),
    
		axis.text  = element_text(color = "#000000", size = 11),
	    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		axis.text.y = element_text(hjust = 1, vjust = 0.5),
	    axis.title = element_text(color = "#000000", size = 14, face = "plain"),
		axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
	    axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
		axis.ticks = element_line(color = "#000000", size = 0.5),
	    axis.ticks.length = unit(0.1, 'cm'),

		legend.title = element_blank(),
	    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
		# font_use = "Arial"
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

	if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
		library(Cairo)
		CairoFonts(
			regular = "Arial:style=Regular",
			bold = "Arial:style=Bold",
			italic = "Arial:style=Italic",
			bolditalic = "Arial:style=Bold Italic,BoldItalic",
		)
		source(ggsave_R)
	}

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
		showtext_auto(enable = TRUE)
		font.add( 'Arial', regular = paste0(ttf_dir, "/arial.ttf"), 
		bold = paste0(ttf_dir, "/arialbd.ttf"), italic = paste0(ttf_dir, "/ariali.ttf"), 
			bolditalic = paste0(ttf_dir, "/arialbi_ttf"))
		mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme
}


bar_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	options(scipen = -1)
	mytheme <- theme_bw() + 
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(color = "#000000", size = 0.8),
    
	    axis.text  = element_text(color = "#000000", size = 11),
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	    axis.text.y = element_text(hjust = 1, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 14, face = "plain"),
	    axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
	    axis.ticks = element_line(color = "#000000", size = 0.5),
		axis.ticks.length = unit(0.1, 'cm'),

	    legend.title = element_blank(),
		plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

	if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
		library(Cairo)
		CairoFonts(
			regular = "Arial:style=Regular",
			bold = "Arial:style=Bold",
			italic = "Arial:style=Italic",
			bolditalic = "Arial:style=Bold Italic,BoldItalic"
		)
		source(ggsave_R)
	}

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}


pie_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	mytheme <- theme_void() +
	theme(
	    legend.title = element_blank(),
		legend.key.height = unit(5, "mm"),
	    legend.key.width = unit(5, "mm"),

		strip.text = element_text(size = 16, face = "plain", hjust = 0.5),
	    plot.title = element_text(size = 20, face = "plain", hjust = 0.5),
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    } 

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}


hist_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	mytheme <- theme_bw() + 
	theme(
	    panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
	    panel.border = element_rect(color = "#000000", size = 0.8),

	    axis.text  = element_text(color = "#000000", size = 11),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
	    axis.text.y = element_text(hjust = 1, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 14, face = "plain"),
	    axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
	    axis.ticks = element_line(color = "#000000", size = 0.5),
		axis.ticks.length = unit(0.1, 'cm'),

	    legend.title = element_blank(),
		strip.text = element_blank(),
	    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    }     

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}

heatmap_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	mytheme <- theme_void() +
	theme(
		axis.text  = element_text(color = "#000000", size = 11),
	    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
		axis.text.y = element_text(hjust = 0, vjust = 0.5),

	    legend.title = element_blank(),
		legend.justification = c(0.5, 0.5),
		legend.margin = margin(5, 5, 5, 5),
		legend.text = element_text(hjust = 0.5, vjust = 0.5, margin = margin(2, 0, 0, 0)),

	    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    } 

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}


heatmap_theme_ggcor <- function(font_use = "Arial") {

	library(ggplot2)
	library(ggcor)
	mytheme <- theme_cor() +
	theme(
		axis.text = element_text(color = "#000000", size = 11),
	    legend.title = element_blank(),
		legend.text = element_text(hjust = 0.5, vjust = 0.5, margin = margin(2, 0, 0, 0)),
	    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    }     

	if( FALSE & "showtext" %in% installed.packages() ) {
	    library(showtext)
		showtext_auto(enable = TRUE)
	    font_add("Arial", regular = arial_ttf)
		mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}


dot_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	mytheme <- theme_bw() +
	theme(
		panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_rect(color = "#000000", size = 0.8),

		axis.text = element_text(color = "#000000", size = 11),
	    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
	    axis.title = element_text(color = "#000000", size = 14, face = "plain"),
		axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
	    axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
		axis.ticks = element_line(color = "#000000", size = 0.5),
	    axis.ticks.length = unit(0.1, 'cm'),

		legend.title = element_blank(),
	    legend.text = element_text(size = 11),
		plot.title = element_text(color = "#000000", size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    } 

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}


dot_theme_manhattan <- function( font_use = "Arial" ) {

	library(ggplot2)
	mytheme <- theme_bw() +
	theme(
		panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		panel.border = element_rect(color = "#000000", size = 0.8),

	    axis.text = element_text(color = "#000000", size = 11),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
	    axis.text.y = element_text(hjust = 1, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 14, face = "bold"),
	    axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
	    axis.ticks = element_line(color = "#000000", size = 0.5),
		axis.ticks.length = unit(0.1, 'cm'),

	    legend.title = element_blank(),
		legend.justification = "center",
	    plot.title = element_text(color = "#000000", size = 16, face = "bold", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

	mytheme

}


box_theme_default <- function(font_use = "Arial") {

	library(ggplot2)
	mytheme <- theme_bw() +
	theme(
		panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		panel.border = element_rect(color = "#000000", size = 0.8),

	    axis.text = element_text(color = "#000000", size = 11),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
	    axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 14, face = "plain"),
	    axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
	    axis.ticks = element_line(color = "#000000", size = 0.5),
		axis.ticks.length = unit(0.11, 'cm'),

	    legend.title = element_blank(),
		legend.justification = "center",
	    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    } 

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}

box_theme_splitviolin <- function(font_use = "Arial") {

	library(ggplot2)
	mytheme <- theme_bw() +
	theme(
	    panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
	    panel.border = element_rect(color = "#000000", size = 0.8),

	    axis.text = element_text(color = "#000000", size = 11),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
	    axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 14, face = "plain"),
	    axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
	    axis.ticks = element_line(color = "#000000", size = 0.5),
		axis.ticks.length = unit(0.11, 'cm'),

	    legend.title = element_blank(),
		legend.justification = "center",
	    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(5,5,5,5), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    } 

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		#font_add("Arial", regular = arial_ttf)
	    #mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}

