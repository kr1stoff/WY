package com.kmlab.module;

import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.ScriptGenerator;

public class FakeReadGenerator {
    private static final Logger logger = LogManager.getLogger(FakeReadGenerator.class);

    private final String RESULT_DIRECTORY;
    private final Map<String, String> PARAMETERS_PROPERTIY_MAP;
    private final int PARALLEL_NUMBER;
    private final List<String> PRIMARY_ACCESSION_NUMBERS;

    public FakeReadGenerator(
            String resultDirectory,
            Map<String, String> pararmetersPropertyMap,
            int parallelNumber,
            List<String> primaryAccessionNumbers) {
        this.RESULT_DIRECTORY = resultDirectory;
        this.PARAMETERS_PROPERTIY_MAP = pararmetersPropertyMap;
        this.PARALLEL_NUMBER = parallelNumber;
        this.PRIMARY_ACCESSION_NUMBERS = primaryAccessionNumbers;
    }

    public void generateFakeReads() {
        logger.info("批量生成 Fake Fastq");
        String genomeDirectory = Paths.get(RESULT_DIRECTORY, ".genomes").toString();
        String logDirectory = Paths.get(RESULT_DIRECTORY, ".log").toString();
        String resourceScriptFile = "/script/snippy_fake_reads.pl";
        String targetScriptPath = Paths.get(RESULT_DIRECTORY, ".script/bin", "snippy_fake_reads.pl").toString();

        ScriptGenerator.generateScript(resourceScriptFile, targetScriptPath);
        BioinformaticToolsExecutor.snippyFakeReads(PRIMARY_ACCESSION_NUMBERS, genomeDirectory, PARAMETERS_PROPERTIY_MAP,
                RESULT_DIRECTORY, logDirectory, PARALLEL_NUMBER);
    }
}
