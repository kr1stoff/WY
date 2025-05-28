library("officer")
library("flextable")
library("magrittr")
library("this.path")

# library('docxtractr')
# convert_to_pdf(outfile, gsub('.docx', '.pdf', outfile))

# 表格字体字号格式
table_font_property <- function(ft) {
  ft <- fontsize(ft, size = 9, part = "header")
  ft <- fontsize(ft, size = 9, part = "body")
  ft <- font(ft, fontname = "微软雅黑", part = "header")
  ft <- font(ft, fontname = "微软雅黑", part = "body")
  ft <- bold(ft, bold = T, part = "header") # 表头加粗
  return(ft)
}

# 未发现表格
weifaxian <- function(doc, key_word, describe = "未发现") {
  x <- describe
  ft <- flextable(as.data.frame(x))
  ft <- table_font_property(ft) # 表格字体属性
  ft <- align(ft, align = "left", part = "all") # 居中
  ft <- delete_part(ft, part = "header") # 删除表头
  # 设置单元格的高度和宽度
  ft <- width(ft, j = 1, width = 6.9)
  ft <- height(ft, i = 1, height = .5)
  ft <- height(ft, height = .2, part = "header")
  # 设置边框和线条颜色
  big_border <- fp_border(color = "#81D2E6", width = 1)
  border_v <- fp_border(color = "#81D2E6", width = 1)
  border_h <- fp_border(color = "#81D2E6", width = 1)
  ft <- border_remove(x = ft)
  ft <- border_outer(ft, part = "all", border = big_border)
  ft <- border_inner_h(ft, part = "all", border = border_h)
  ft <- border_inner_v(ft, part = "all", border = border_v)
  cursor_reach(doc, keyword = key_word) %>% body_add_flextable(value = ft) # 在对应变量的下方插入表格
  cursor_reach(doc, keyword = key_word) %>% body_remove() # 删除表格标签
}

# 表格宽度调整 宽度和大概6.5左右合适
table_width <- function(ft, key_word) {
  if (key_word == "Table_FenXing") {
    ft <- width(ft, j = c(1, 2), width = 0.9)
    ft <- width(ft, j = c(3, 5), width = 2)
    ft <- width(ft, j = 4, width = 0.7)
  } else if (key_word == "Table_FenXingBatch") {
    ft <- width(ft, j = 1, width = 0.7)
    ft <- width(ft, j = c(2, 3), width = 0.8)
    ft <- width(ft, j = c(4, 6), width = 1.8)
    ft <- width(ft, j = 5, width = 0.6)
  } else if (key_word == "Table_BianYi") {
    ft <- width(ft, j = 1, width = 1)
    ft <- width(ft, j = c(2:4, 8, 9), width = 0.6)
    ft <- width(ft, j = 5, width = 0.8)
    ft <- width(ft, j = c(6, 7), width = 0.9)
  } else if (key_word == "Table_ShuJuTongJi") {
    ft <- width(ft, j = c(1:3), width = 1.1)
    ft <- width(ft, j = 4, width = 1.2)
    ft <- width(ft, j = c(5, 6), width = 1)
  } else if (key_word == "Table_BiDuiJieGuo") {
    ft <- width(ft, width = 0.9) # 7*0.9=6.3
  } else {
    ft <- set_table_properties(ft, layout = "autofit")
  }
  return(ft)
}

# 普通表格
normal_table <- function(doc, data, key_word) {
  ft <- flextable(data)
  ft <- table_font_property(ft) # 表格字体属性
  ft <- align(ft, align = "center", part = "all") # 居中
  # 配置表头的颜色
  ft <- bg(ft, bg = "#80D1E6", part = "header") # 008CD6
  # 设置单元格的高度和宽度
  ft <- table_width(ft, key_word)
  ft <- height(ft, i = 1, height = .5)
  ft <- height(ft, height = .2, part = "header")
  # 设置边框和线条颜色
  big_border <- fp_border(color = "#C6D9F0")
  border_v <- fp_border(color = "#C6D9F0")
  border_h <- fp_border(color = "#C6D9F0")
  ft <- border_remove(x = ft)
  ft <- border_outer(ft, part = "all", border = big_border)
  ft <- border_inner_h(ft, part = "all", border = border_h)
  ft <- border_inner_v(ft, part = "all", border = border_v)
  cursor_reach(doc, keyword = key_word) %>% body_add_flextable(value = ft) # 在对应变量的下方插入表格
  cursor_reach(doc, keyword = key_word) %>% body_remove() # 删除表格标签
}

# 使用哪个表格,'正常或未发现'表格
assign_table <- function(doc, tab, key_word, describe = "未发现") {
  dataframe <- read.table(tab, sep = "\t", header = T, stringsAsFactors = F, encoding = "UTF-8", check.names = F, colClasses = "character")
  if (nrow(dataframe) < 1) {
    weifaxian(doc, key_word, describe)
  } else {
    normal_table(doc, dataframe, key_word)
  }
}

# 当前样本Pangolin结果
single_pangolin_table <- function(doc, dir_upload_sample, sample_number) {
  dataframe <- read.table(paste0(dir_upload_sample, "/source/lineage_report_trans.xls"),
    sep = "\t", header = T, stringsAsFactors = F, encoding = "UTF-8", colClasses = "character"
  )
  normal_table(doc, dataframe[dataframe$样本名 == sample_number, 2:ncol(dataframe)], key_word = "Table_FenXing")
}

# 图像
# 230117 FASTQ2不存在的图片删除标签
import_image <- function(doc, key_word, path_png, width = 6, height = 3) {
  if (file.exists(path_png)) {
    cursor_reach(doc, keyword = key_word) %>% body_add_img(path_png, pos = "before", width = width, height = height, style = "center")
  }
  cursor_reach(doc, keyword = key_word) %>% body_remove()
}

# 更新文本类型的标签
update_text <- function(doc, sample_info) {
  df_text <- read.table(sample_info, sep = "\t", encoding = "UTF-8", header = F, stringsAsFactors = F, row.names = 1, colClasses = "character") %>%
    t() %>%
    as.data.frame()
  attach(df_text)
  body_replace_all_text(doc, "DanWei", DanWei, only_at_cursor = F, fixed = T)
  body_replace_all_text(doc, "KeShi", KeShi, only_at_cursor = F, fixed = T)
  body_replace_all_text(doc, "SampleNumber", SampleNumber, only_at_cursor = F, fixed = T)
  body_replace_all_text(doc, "ReportNumber", ReportNumber, only_at_cursor = F, fixed = T)
  body_replace_all_text(doc, "SongJianDate", SongJianDate, only_at_cursor = F, fixed = T)
  body_replace_all_text(doc, "SampleType", SampleType, only_at_cursor = F, fixed = T)
  detach(df_text)
}

# 更新表格类型的标签
update_table <- function(doc, dir_upload_sample, sample_number) {
  assign_table(doc, paste0(dir_upload_sample, "/1.qc/", sample_number, ".basic.stat.txt"), key_word = "Table_ShuJuTongJi")
  assign_table(doc, paste0(dir_upload_sample, "/2.map/", sample_number, ".bam_stats.txt"), key_word = "Table_BiDuiJieGuo")
  single_pangolin_table(doc, dir_upload_sample, sample_number)
  assign_table(doc, paste0(dir_upload_sample, "/source/lineage_report_trans.xls"), key_word = "Table_FenXingBatch")
  assign_table(doc, paste0(dir_upload_sample, "/3.variant/", sample_number, ".trans.tsv"), key_word = "Table_BianYi")
}

# 更新图片类型的标签
update_image <- function(doc, dir_upload_sample, sample_number) {
  import_image(doc, "Image_GuoLvQian_Fastqc1", paste0(dir_upload_sample, "/1.qc/", sample_number, ".1_before_per_base_quality.png"), 4, 2)
  import_image(doc, "Image_GuoLvQian_Fastqc2", paste0(dir_upload_sample, "/1.qc/", sample_number, ".2_before_per_base_quality.png"), 4, 2)
  import_image(doc, "Image_GuoLvHou_Fastqc1", paste0(dir_upload_sample, "/1.qc/", sample_number, ".1_after_per_base_quality.png"), 4, 2)
  import_image(doc, "Image_GuoLvHou_Fastqc2", paste0(dir_upload_sample, "/1.qc/", sample_number, ".2_after_per_base_quality.png"), 4, 2)
  import_image(doc, "Image_FuGaiDu", paste0(dir_upload_sample, "/2.map/", sample_number, ".genome_coverage_depth.png"))
  import_image(doc, "Image_1000FuGaiDu", paste0(dir_upload_sample, "/2.map/", sample_number, ".genome_coverage_depth_ylim1000.png"))
  import_image(doc, "Image_SingleTree", paste0(dir_upload_sample, "/5.phylogenetic/rectangular.png"), 5.8, 2.9)
  import_image(doc, "Image_SNPTree", paste0(dir_upload_sample, "/source/SNPTree/rectangular.png"), 5.8, 2.9)
  import_image(doc, "Image_MSATree", paste0(dir_upload_sample, "/source/MSATree/rectangular.png"), 5.8, 2.9)
}

which_template <- function(company) {
  if (company == "KP") {
    template <- gsub("bin", "templates/新冠模板_凯普.docx", this.path::this.dir())
  } else {
    template <- gsub("bin", "templates/新冠模板_WY.docx", this.path::this.dir())
  }
  return(template)
}

# 命令行参数
get_args <- function() {
  args <- commandArgs(trailingOnly = T)
  print_usage <- function() {
    cat("PROGRAM <样本upload目录> <样本信息表> <模板公司>\n")
    quit(save = "no", status = 1)
  }
  if (length(args) != 3) {
    cat("Need 3 Arguments!\n")
    print_usage()
  } else if (!file.exists(args[1]) | !file.exists(args[2])) {
    cat("<depth file> or <summary_numbers_file> not exists!\n")
    print_usage()
  }
  return(args)
}

# 配置
args <- get_args()
dir_upload_sample <- args[1]
sample_info <- args[2]
company <- args[3]
template <- which_template(company)
sample_number <- basename(dir_upload_sample)
report_date <- gsub("-", "", as.character(Sys.Date()))
outfile <- paste0(dir_upload_sample, "/", sample_number, "-新冠病毒全基因组测序报告-", report_date, ".docx")
# 执行
doc <- read_docx(template)
update_text(doc, sample_info)
update_table(doc, dir_upload_sample, sample_number)
update_image(doc, dir_upload_sample, sample_number)
print(doc, outfile) # print输出word
