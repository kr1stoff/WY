package com.kmlab.module;

import java.nio.file.Paths;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.ScriptGenerator;

public class AlignmentInformationStats {
    private static final Logger logger = LogManager.getLogger(AlignmentInformationStats.class);

    private final String RESULT_DIRECTORY;

    public AlignmentInformationStats(
            String resultDirectory) {
        this.RESULT_DIRECTORY = resultDirectory;
    }

    public void stats() {
        logger.info("统计比对结果");
        String resourceScriptFile = "/script/alignment_information_stats.py";
        String targetScriptPath = Paths.get(RESULT_DIRECTORY, ".script/bin", "alignment_information_stats.py")
                .toString();

        ScriptGenerator.generateScript(resourceScriptFile, targetScriptPath);
        BioinformaticToolsExecutor.alignmentStats(RESULT_DIRECTORY);
    }
}
