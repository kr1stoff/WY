package com.kmlab.module;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import com.kmlab.util.FileOperator;

public class ReferenceSequenceProcessor {
    private static final Logger logger = LogManager.getLogger(ReferenceSequenceProcessor.class);

    private final String REFERENCE_FILE;
    private final String KINGDOM;
    private final int THREAD_NUMBER;
    private final int PARALLEL_NUMBER;
    private final int SINGLE_TASK_THREAD;
    private final List<String> PRIMARY_ACCESSION_NUMBERS;
    private final String RESULT_DIRECTORY;
    private final Map<String, String> PARAMETERS_PROPERTIY_MAP;
    private final String LOG_DIRECTORY;

    public ReferenceSequenceProcessor(
            String referenceFile,
            String resultDirectory,
            int threadNumber,
            String kingdom,
            int parallelNumber,
            List<String> primaryAccessionNumbers,
            Map<String, String> pararmetersPropertyMap) {
        this.REFERENCE_FILE = referenceFile;
        this.RESULT_DIRECTORY = resultDirectory;
        this.KINGDOM = kingdom;
        this.THREAD_NUMBER = threadNumber;
        this.PARALLEL_NUMBER = parallelNumber;
        this.PRIMARY_ACCESSION_NUMBERS = primaryAccessionNumbers;
        this.PARAMETERS_PROPERTIY_MAP = pararmetersPropertyMap;
        this.LOG_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".log").toString();
        this.SINGLE_TASK_THREAD = THREAD_NUMBER / PARALLEL_NUMBER;
    }

    /**
     * 分配参考序列访问码。
     * 本方法用于根据提供的参考基因组信息，生成或选择一个参考序列，并将其访问码分配给后续流程使用。
     * 如果未提供参考基因组文件，则通过一系列生物信息学工具（checkm、prokka和seqtk）对基因组进行质量评估和注释，
     * 并基于评估结果生成一个参考序列访问码。如果已经提供了参考基因组文件，则直接使用该文件的访问码。
     * 最后，会为使用的参考序列建立一个符号链接，供后续流程使用。
     */
    public void assignReferenceAccession() {

        Path linkReferencePath = Paths.get(RESULT_DIRECTORY, ".reference/ref.fna");
        String genomeDirectory = Paths.get(RESULT_DIRECTORY, ".genomes").toString();

        if (REFERENCE_FILE == "") {
            // 完整度,污染度,N碱基,假基因排序
            logger.info("未提供参考基因组登录号, 计算并指定参考序列登录号");

            BioinformaticToolsExecutor.checkmGenomes(RESULT_DIRECTORY, LOG_DIRECTORY, THREAD_NUMBER, genomeDirectory);
            BioinformaticToolsExecutor.prokkaGenomes(RESULT_DIRECTORY, PRIMARY_ACCESSION_NUMBERS, genomeDirectory,
                    SINGLE_TASK_THREAD, KINGDOM, LOG_DIRECTORY, PARALLEL_NUMBER);
            BioinformaticToolsExecutor.seqtkCompGenomes(RESULT_DIRECTORY, PARALLEL_NUMBER, PRIMARY_ACCESSION_NUMBERS,
                    genomeDirectory, LOG_DIRECTORY);

            GenomeQualityCalculator genomeQualityCalculator = new GenomeQualityCalculator(
                    RESULT_DIRECTORY, PARAMETERS_PROPERTIY_MAP);
            String referenceAccession = genomeQualityCalculator.calculateGenomeQuality();

            Path sourcePath = Paths.get(RESULT_DIRECTORY, ".genomes", referenceAccession + ".fna");
            FileOperator.linkFile(sourcePath, linkReferencePath);
        } else {
            Path sourcePath = Paths.get(REFERENCE_FILE);
            FileOperator.linkFile(sourcePath, linkReferencePath);
        }

        String logFile = Paths.get(LOG_DIRECTORY, "bwa_index_reference.log").toString();
        BioinformaticToolsExecutor.bwaIndex(linkReferencePath.toString(), logFile);
    }
}
