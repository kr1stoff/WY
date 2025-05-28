## hongshimiao 's ggplot.r setttings

## 图片宽度
FigWidth = function(sample.name, width = 8) {
    sample.num = length(sample.name)
    if(sample.num >= 800) {
        width = 30
    } else if(sample.num >= 500) {
        width = 25
    } else if(sample.num >= 100) {
        width = switch(as.integer(sample.num / 100), 11, 14, 17, 20)
    } else {
        if(sample.num >= 50) {
            width = 10
        } else if(sample.num >= 20) {
            width = 9
        }
    }

    width
}

## 字体大小
AxisTextSize <- function(value, width = 8) {
    name <- as.character(unique(value))
    name.len <- sum(nchar(name))
    name.meanlen <- name.len / length(name)

    if(name.meanlen <= 10) {
        size = 12
    } else if(name.meanlen <= 20) {
        size = 11
    } else {
        size = 10
    }

    name.num <- length(name)
    if(name.num >= 30) {
        size = width * 0.6 * 72 / name.num
    }

    size
}

## 柱形图 添加文字大小
GeomTextSize = function(value, sample.num = 10, width = 8, fontsize = 4, legend = T) {
    name = as.character(value)
    name.len = sum(nchar(name))
    name.maxlen = max(nchar(name))
    name.meanlen = name.len / length(name)

    name.font.len = fontsize * 1 / 72 * name.len * 1
    name.max.font.len = fontsize * 1 / 72 * name.maxlen * 1
    draw.width = ifelse(legend, width - 3, width - 2)

    if(name.max.font.len >= draw.width / sample.num * 0.8) {
        fontsize = draw.width * 72 / sample.num / name.maxlen * 0.5
    } else if(name.max.font.len >= draw.width / (sample.num + 2) * 0.8) {    # recheck
		fontsize = draw.width * 72 / sample.num / name.maxlen * 0.5
	}

    fontsize
}

## 根据名字长度设置角度 axis.text
## xlab text names 
xlab_angle = function(sample_name, width = 7, fontsize = 14) { 
    name = as.character(unique(sample_name))
    name_len = sum(nchar(name))
    name_mean = name_len / length(name)

    name_font_len = fontsize * 1/72 * name_len * 0.6  ## size * pt2in * len * 大小写比例,长度不一

    cut = c(0.6, 2) * width  ## 文字 估算占 图片的0.6  放松一些, 0.7 名字会连在一起
    
	if(name_font_len < cut[1]) {
        angle = 0   ## 水平
        hjust = 0.5 ## 左右 居中
        vjust = 0   ## 貌似没用
    } else if(name_font_len < cut[2]) {
        if(name_mean > 15){ ## 名字过长，角度太大不好看
            angle = 45
            hjust = 1
            vjust = 1
        } else {
            angle = 45
            hjust = 1  ## 左右, 最右
            vjust = 1  ## 上下, 最上
        }
    } else {
        angle = 45  ## 竖直
        hjust = 1   ## 上下, 最上
        vjust = 1   ## 左右, 中间
    }

    theme(axis.text.x = element_text(angle=angle, hjust = hjust, vjust = vjust))
}

SaveFig = function(fig = fig, outpfx = outpfx, width = 7, height = 7) {
	ggsave(file = paste(outpfx,".png",sep=""), fig, width=width, height=height, limitsize = FALSE)
	ggsave(file = paste(outpfx,".pdf",sep=""), fig, width=width, height=height, limitsize = FALSE)
}

