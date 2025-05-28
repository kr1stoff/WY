package com.kmlab.cli;

import java.util.Optional;

import org.apache.commons.cli.Options;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class OptionsParser {
    private static final Logger logger = LogManager.getLogger(OptionsParser.class);
    private static final Options options = new Options();
    private static final CommandLineParser parser = new DefaultParser();

    static {
        options.addOption(Option.builder("g")
                .longOpt("genomes")
                .hasArg(true)
                .required(true)
                .argName("dir")
                .desc("(必选) 输入基因组目录, 包含一个物种下的多个基因组序列.")
                .build());

        options.addOption(Option.builder("r")
                .longOpt("reference-accession")
                .hasArg(true)
                .required(false)
                .argName("string")
                .desc("(可选) 参考基因组登录号. 例: GCF_008632635.1")
                .build());

        options.addOption(Option.builder("o")
                .longOpt("result-dir")
                .hasArg(true)
                .required(false)
                .argName("dir")
                .desc("(可选) 输出结果路径, 默认: pgi_results.")
                .build());

        options.addOption(Option.builder("t")
                .longOpt("threads")
                .hasArg(true)
                .required(false)
                .argName("int")
                .desc("(可选) 线程数, 最小值8. 默认: 当前系统最大线程数/2.")
                .build());

        options.addOption(Option.builder("k")
                .longOpt("kingdom")
                .hasArg(true)
                .required(false)
                .argName("string")
                .desc("(可选) 物种所属界, 细菌,真菌,病毒. {Bacteria(默认), Fungi, Virus}")
                .build());

        options.addOption(Option.builder()
                .longOpt("pe")
                .required(false)
                .argName("boolean")
                .desc("(可选) 模拟 PE FASTQ 数据用于下游分析. 默认: SE.")
                .build());

        options.addOption(Option.builder()
                .longOpt("inclusiveness-filter")
                .required(false)
                .argName("boolean")
                .desc("(可选) 整合基因组 Contig 是否进行包容性过滤.")
                .build());

        options.addOption(Option.builder()
                .longOpt("assembly-software")
                .hasArg(true)
                .required(false)
                .argName("string")
                .desc("(可选) 组装软件选择. {SPAdes(默认), SOAPdenovo2, MEGAHIT}.")
                .build());

        options.addOption(Option.builder()
                .longOpt("remove-clipping-reads")
                .required(false)
                .argName("boolean")
                .desc("(可选) 删除软剪切序列, 默认: 保留.")
                .build());

        options.addOption(Option.builder()
                .longOpt("remove-plasmid-genome")
                .required(false)
                .argName("boolean")
                .desc("(可选) 基因组筛选, 是否删除包含质粒的基因组.")
                .build());

        options.addOption(Option.builder()
                .longOpt("remove-low-mapped")
                .required(false)
                .argName("boolean")
                .desc("(可选) 基因组筛选, 删除比对率低的基因组.")
                .build());

        options.addOption(Option.builder()
                .longOpt("no-genomes-screening")
                .required(false)
                .argName("boolean")
                .desc("(可选) 不做基因组筛选, 直接生成整合基因组. 默认先做基因组筛选信息表.")
                .build());

        options.addOption("h", "help", false, "显示帮助信息.");
    }

    /**
     * 解析命令行参数。
     *
     * @param args 命令行参数数组。
     * @return 如果解析成功，返回包含解析结果的Optional；如果解析失败，返回空的Optional。
     */
    public static Optional<CommandLine> parse(String[] args) {
        try {
            // 尝试解析命令行参数
            CommandLine comParse = parser.parse(options, args);
            return Optional.of(comParse);
        } catch (MissingOptionException e) {
            // 缺少必要选项时记录错误并返回空Optional
            logger.warn("缺少必要选项");
            return Optional.empty();
        } catch (ParseException e) {
            // 解析异常时打印堆栈信息并返回空Optional
            e.printStackTrace();
            return Optional.empty();
        }
    }

    /**
     * 打印帮助信息。
     * 该方法用于输出关于病原基因组整合工具的帮助信息，包括工具的使用说明和可用选项等。
     * 该方法不接受参数，也不返回任何值。
     */
    public static void printHelp() {
        // 创建HelpFormatter实例用于格式化帮助信息
        HelpFormatter formatter = new HelpFormatter();
        // 打印帮助信息
        formatter.printHelp("pgi", "病原基因组整合工具.\n", options, null);
    }
}
