#!/usr/bin/env Rscript

##加载R包，初始传递命令、变量等
library(getopt)

spec <- matrix(c(
	'input', 'i', 1, 'character',
	'output', 'o', 1, 'character',
	'help', 'h', 0, 'logical'
	), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print_usage <- function(spec = NULL) {
	getopt(spec, usage = TRUE)
	cat('
用法示例: 
1) Rscript depth_GC_plot.r -i [depth_base.stat] -o [output_file_name]
2) Rscript depth_GC_plot.r --input [depth_base.stat] --output [output_file_name]

-i --input	character	按基因组滑窗统计的 reads 覆盖深度及碱基含量结果文件 depth_base.stat [required]
-o --output	character	输出图片名称，输出 png 格式图片 [required]
-h --help	logical		帮助信息
\n')
	q('no')
}
if (!is.null(opt$help)) print_usage(spec)

library(ggplot2)
library(grid)


#读取文件
depth_base <- read.delim(opt$input)

#统计平均 GC 含量（百分比）
GC_mean <- round(mean(depth_base$GC), 2)

#统计平均测序深度（百分比）
depth_mean <- round(mean(depth_base$depth), 2)
depth_base <- subset(depth_base, depth <= 3 * depth_mean)

#depth 深度、GC 含量散点密度图
depth_GC <- ggplot(depth_base, aes(GC, depth)) +
	geom_point(color = 'gray', alpha = 0.6, pch = 19, size = 0.5) +
	geom_vline(xintercept = GC_mean, color = 'red', lty = 2, lwd = 0.5) + 
	geom_hline(yintercept = depth_mean, color = 'red', lty = 2, lwd = 0.5) +
	stat_density2d(aes(fill = ..density.., alpha = ..density..), geom = 'tile', contour = FALSE, n = 500) +
	scale_fill_gradientn(colors = c('transparent', 'gray', 'yellow', 'red')) +
	theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
	labs(x = 'GC%', y = 'Depth') +
	theme(axis.text = element_text(size = 10)) +
	theme(axis.title = element_text(size = 12)) +
	theme(legend.position = 'none')

#depth 深度直方密度图
depth_hist <- ggplot(depth_base, aes(depth)) +
	geom_histogram(binwidth = (max(depth_base$depth) - min(depth_base$depth))/100, fill = 'gray', color = 'gray40', size = 0.1) +
	geom_rug(color = 'gray', alpha = 0.6) +
	theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
	theme(axis.line = element_line(color = 'black', size = 0.3), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
	labs(x = '', y = 'Numbers') +
	coord_flip() +
	geom_vline(xintercept = depth_mean, color = 'red', lty = 2, lwd = 0.5)

#GC 含量直方密度图
GC_hist <- ggplot(depth_base, aes(GC)) +
	geom_histogram(binwidth = (max(depth_base$GC) - min(depth_base$GC))/100, fill = 'gray', color = 'gray40', size = 0.1) +
	geom_rug(color = 'gray', alpha = 0.6) +
	theme(panel.grid.major = element_line(color = 'gray', linetype = 2, size = 0.25), panel.background = element_rect(color = 'black', fill = 'transparent')) +
	theme(axis.line = element_line(color = 'black', size = 0.3), axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
	labs(x = '', y = 'Numbers') +
	geom_vline(xintercept = GC_mean, color = 'red', lty = 2, lwd = 0.5)

#组合图片并输出
pdf(paste(opt$output, 'pdf', sep = '.'), width = 8, height = 8)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(3, 3)))
	print(depth_GC, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 1:2))
	print(GC_hist, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
	print(depth_hist, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 3))
dev.off()

#png(paste(opt$output, 'png', sep = '.'), width = 4000, height = 4000, res = 600, units = 'px')
#	grid.newpage()
#	pushViewport(viewport(layout = grid.layout(3, 3)))
#	print(depth_GC, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 1:2))
#	print(GC_hist, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
#	print(depth_hist, vp = viewport(layout.pos.row = 2:3, layout.pos.col = 3))
#dev.off()
