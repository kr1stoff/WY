package com.kmlab.cli;

import java.util.Optional;
import java.util.Map;
import java.util.HashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.commons.cli.CommandLine;

public class OptionsToHashMap {
    private static final Logger logger = LogManager.getLogger(OptionsToHashMap.class);

    /**
     * 解析命令行参数并将其填充到一个Map中。
     *
     * @param args 命令行参数数组。
     * @return 如果解析成功，返回一个包含所有选项和其值的Map；如果解析失败或缺少必要参数，返回一个空的Optional。
     */
    public static Optional<Map<String, Object>> map(String[] args) {
        Map<String, Object> mapOptions = new HashMap<>();

        try {
            Optional<CommandLine> commandLine = OptionsParser.parse(args);
            if (!commandLine.isPresent()) {
                OptionsParser.printHelp();
                return Optional.empty();
            }
            fillMapOptions(commandLine.get(), mapOptions);
            OptionsValidator.validateGenomeDirAndRefFile((String) mapOptions.get("dirGenomes"));
            mapOptions.forEach((key, value) -> logger.info(key + ": " + value));
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("获取命令行参数是发生错误", e);
        }
        return Optional.ofNullable(mapOptions);
    }

    /**
     * 将命令行参数填充到选项映射中。
     *
     * @param commandLine 命令行对象，包含已解析的命令行参数。
     * @param mapOptions  用于存储选项名称和其对应值的映射。
     */
    private static void fillMapOptions(CommandLine commandLine, Map<String, Object> mapOptions) {
        mapOptions.put("dirGenomes", commandLine.getOptionValue("g"));
        mapOptions.put("dirResult", commandLine.getOptionValue("o", "pgi_results"));
        mapOptions.put("fileReference", OptionsCorrector.getReferenceAbsolutePath(commandLine));
        mapOptions.put("assemblySoftware", OptionsCorrector.correctAssemblySoftware(commandLine));
        mapOptions.put("sequenceMode", commandLine.hasOption("pe") ? "PE" : "SE");
        mapOptions.put("boolInclusivenessFilter", commandLine.hasOption("inclusiveness-filter"));
        mapOptions.put("boolRemoveClippingReads", commandLine.hasOption("remove-clipping-reads"));
        mapOptions.put("boolRemovePlasmidGenome", commandLine.hasOption("remove-plasmid-genome"));
        mapOptions.put("boolRemoveLowMapped", commandLine.hasOption("remove-low-mapped"));
        mapOptions.put("boolNoGenomesScreening", commandLine.hasOption("no-genomes-screening"));
        mapOptions.put("threads", OptionsCorrector.getThreads(commandLine));
        mapOptions.put("kingdom", OptionsCorrector.correctKingdom(commandLine));
    }
}
