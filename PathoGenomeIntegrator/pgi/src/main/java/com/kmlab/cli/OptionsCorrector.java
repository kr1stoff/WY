package com.kmlab.cli;

import java.util.Arrays;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.FileOperator;
import com.kmlab.util.StringProcessor;
import com.kmlab.util.SystemInfo;

public class OptionsCorrector {
    private static final Logger logger = LogManager.getLogger(OptionsCorrector.class);

    /**
     * 获取参考基因组的绝对路径。
     * 此方法首先检查是否通过命令行提供了参考基因组的访问编号（accession），然后验证accession的格式。
     * 如果格式正确，将尝试在指定的目录中找到对应的基因组文件，并返回其绝对路径。
     *
     * @param commandLine 包含命令行参数的对象，预期包含"r"和"g"两个选项。
     *                    "r"选项代表参考基因组的accession号，"g"选项代表基因组文件所在的目录。
     * @return 返回参考基因组文件的绝对路径。如果未提供参考基因组accession或accession格式不正确，则返回空字符串。
     */
    public static String getReferenceAbsolutePath(CommandLine commandLine) {
        // 定义用于检查accession格式的正则表达式
        final Pattern CHECK_PATTERN = Pattern.compile("^(GCF|GCA)_[0-9]+$");
        String rawReferenceAccession = commandLine.getOptionValue("r");
        String dirGenomes = commandLine.getOptionValue("g");

        // 检查是否提供了参考基因组的accession
        if (rawReferenceAccession == null) {
            logger.info("未提供 --reference-accession 参数. 默认使用使 checkM 计算最合适的参考基因组.");
            return "";
        } else if (!CHECK_PATTERN.matcher(rawReferenceAccession).matches()) {
            // 检查accession的格式
            logger.info("--reference-accession 参数格式错误，请使用 GCF/GCA_XXXXX 格式. 默认使用使 checkM 计算最合适的参考基因组.");
            return "";
        }

        // 在指定目录中查找参考基因组文件
        String[] extensions = { "fna", "fa", "fasta" };
        String strReferenceAbsolutePath = FileOperator.findFileByContentAndExtensions(dirGenomes,
                rawReferenceAccession.toUpperCase(), extensions);
        logger.info("参考基因组绝对路径 " + strReferenceAbsolutePath);
        return strReferenceAbsolutePath;
    }

    /**
     * 校正过滤长度选项的值。
     * 该方法用于校验通过命令行传递的过滤长度参数，并确保其值符合特定的条件。
     * 如果传入的值无效或不在指定的范围内（300到2000之间），则会使用默认值N90。
     *
     * @param commandLine 命令行对象，包含用户提供的过滤长度选项。
     * @return 经校正后的过滤长度选项值，要么是用户输入的有效值，要么是默认值N90。
     */
    public static String correctFilterLengthOption(CommandLine commandLine) {
        String rawFilterLength = commandLine.getOptionValue("filter-length", "N90").toUpperCase();
        String regex = "(?i)N90|\\d{3,4}";
        String filterLength = rawFilterLength;
        if (!rawFilterLength.matches(regex)) {
            logger.warn("--filter-length 无效参数: " + rawFilterLength + ", 默认使用N90.");
        } else if (rawFilterLength.equalsIgnoreCase("N90")) {
        } else {
            int length = Integer.parseInt(rawFilterLength);
            if (!(length >= 300 && length <= 2000)) {
                logger.warn("无效的过滤长度 (length >= 300 && length <= 2000): " + rawFilterLength + ", 改用默认值N90");
                filterLength = "N90";
            } else {
                filterLength = rawFilterLength;
            }
        }

        return filterLength;
    }

    /**
     * 校正组装软件的名称。
     * 该方法从命令行中获取组装软件的参数，如果参数无效，则默认使用 SPAdes。
     *
     * @param commandLine 包含组装软件参数的命令行对象。
     * @return 经校正后的有效组装软件名称。
     */
    public static String correctAssemblySoftware(CommandLine commandLine) {
        String rawAssemblySoftware = commandLine.getOptionValue("assembly-software", "spades").toLowerCase();
        String[] validAssemblySoftware = { "spades", "megahit", "soapdenovo2" };
        String assemblySoftware = rawAssemblySoftware;
        if (!Arrays.asList(validAssemblySoftware).contains(rawAssemblySoftware)) {
            logger.warn("无效的组装软件名: " + rawAssemblySoftware + ", 默认使用 SPAdes.");
            assemblySoftware = "spades";
        }

        return assemblySoftware;
    }

    /**
     * 从命令行参数中获取删除低比对率基因组的阈值。
     *
     * @param commandLine 命令行对象，包含各种运行时参数。
     * @return 返回一个浮点数，表示删除低比对率基因组的比例阈值。默认值为0.1f。
     */
    public static float getRemoveLowMappedThreshold(CommandLine commandLine) {
        String rawStringTholdRemoveLowMapped = commandLine.getOptionValue("remove-low-mapped-threshold", "0.1f");
        float tholdRemoveLowMapped = (float) 0.1;
        float rawTholdRemoveLowMapped = Float.parseFloat(rawStringTholdRemoveLowMapped);

        if (rawTholdRemoveLowMapped > 1) {
            logger.warn("删除比对率低的基因组比例阈值参数需要小于1: " + rawStringTholdRemoveLowMapped);
        } else {
            tholdRemoveLowMapped = rawTholdRemoveLowMapped;
        }

        logger.warn("无效删除比对率低的基因组比例阈值参数: " + rawStringTholdRemoveLowMapped);

        return tholdRemoveLowMapped;
    }

    /**
     * 获取线程数。
     * 此方法根据提供的命令行参数来获取线程数。首先尝试从命令行中读取"threads"参数的值，
     * 如果未提供该参数或者提供的值无效（不是数字或者小于8），则会使用系统当前最大线程数的一半作为默认值。
     *
     * @param commandLine 命令行对象，用于获取"threads"参数的值。
     * @return 返回有效线程数的字符串表示。要么是用户提供的合法数值，要么是系统默认计算的值。
     */
    public static String getThreads(CommandLine commandLine) {
        String rawStringThreadNumber = commandLine.getOptionValue("threads");
        String defaultStringThreadNumber = SystemInfo.getThreadNumber() / 2 + "";

        if (rawStringThreadNumber == null) {
            logger.info("未提供 --threads 参数. 默认使用当前系统最大线程数/2");
            return defaultStringThreadNumber;
        } else if (!rawStringThreadNumber.matches("\\d+")) {
            logger.warn("无效的线程数参数: " + rawStringThreadNumber + ", 默认使用当前系统最大线程数/2");
            return defaultStringThreadNumber;
        } else if (Integer.parseInt(rawStringThreadNumber) < 8) {
            logger.warn("输入线程数小于8: " + rawStringThreadNumber + ", 默认使用当前系统最大线程数/2");
            return defaultStringThreadNumber;
        }
        return rawStringThreadNumber;
    }

    /**
     * 根据命令行提供的参数来校正并返回生物界的名称。
     * 参数 "kingdom" 用于指定生物界，如果不提供则默认为 "bacteria"。
     * 有效的生物界名称包括 "bacteria", "fungi", "virus"。
     * 如果提供的生物界名称无效，则记录警告信息并使用默认值 "bacteria"。
     *
     * @param commandLine 命令行对象，包含用户指定的生物界参数。
     * @return 校正后的生物界名称。
     */
    public static String correctKingdom(CommandLine commandLine) {
        String rawKingdom = commandLine.getOptionValue("kingdom", "bacteria").toLowerCase();
        String[] validKingdom = { "bacteria", "fungi", "virus" };
        String kingdom = rawKingdom;

        if (!Arrays.asList(validKingdom).contains(rawKingdom)) {
            logger.warn("无效的界名: " + rawKingdom + ", 默认使用 Bacteria.");
            kingdom = "bacteria";
        }

        // 在这里就把 kingdom 转成 title 格式 (Bactreia)
        return StringProcessor.toTitle(kingdom);
    }
}
