package com.kmlab.module;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import com.kmlab.util.YamlWriter;

public class SPAdesGenomeAssembler {
    private static final Logger logger = LogManager.getLogger(SPAdesGenomeAssembler.class);

    private final List<String> PRIMARY_ACCESSION_NUMBERS;
    private final String RESULT_DIRECTORY;
    private final int THREAD_NUMBER;

    public SPAdesGenomeAssembler(
            List<String> primaryAccessionNumbers,
            String resultDirectory,
            int threadNumber) {
        this.PRIMARY_ACCESSION_NUMBERS = primaryAccessionNumbers;
        this.RESULT_DIRECTORY = resultDirectory;
        this.THREAD_NUMBER = threadNumber;
    }

    public void runSPAdes() {
        String datasetsYaml = createSingleReadsDatasets();
        BioinformaticToolsExecutor.runSPAdesGenomeAssembler(RESULT_DIRECTORY, datasetsYaml, THREAD_NUMBER);
    }

    public String createSingleReadsDatasets() {
        logger.info("生成 SPAdes 输入 datasets.yaml 文件");
        String datasetsYaml = Paths.get(RESULT_DIRECTORY, "unmapped_reads_assembly/datasets.yaml").toString();
        List<Map<String, Object>> datasetsList = new ArrayList<>();
        Map<String, Object> datasetsMap = new HashMap<>();
        datasetsMap.put("type", "single");
        List<String> singleReads = new ArrayList<>();

        for (String accessionNumber : PRIMARY_ACCESSION_NUMBERS) {
            String fastqPath = Paths
                    .get(RESULT_DIRECTORY, "align_fastq_to_reference", accessionNumber,
                            accessionNumber + ".unmapped.fq")
                    .toString();
            File fastqFile = new File(fastqPath);
            if (fastqFile.length() == 0) {
                logger.debug(fastqPath + " 是空文件");
            } else {
                singleReads.add(fastqPath);
            }
        }

        datasetsMap.put("single reads", singleReads);
        datasetsList.add(datasetsMap);

        YamlWriter.write(datasetsList, datasetsYaml);

        return datasetsYaml;
    }
}
