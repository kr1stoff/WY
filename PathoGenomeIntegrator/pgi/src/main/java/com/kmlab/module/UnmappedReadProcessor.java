package com.kmlab.module;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.FileOperator;
import com.kmlab.util.TXTWriter;

public class UnmappedReadProcessor {
    private static final Logger logger = LogManager.getLogger(UnmappedReadProcessor.class);

    private final List<String> PRIMARY_ACCESSION_NUMBERS;
    private final String RESULT_DIRECTORY;
    private final int PARALLEL_NUMBER;
    private final int SINGLE_TASK_THREAD;

    public UnmappedReadProcessor(
            List<String> primaryAccessionNumbers,
            String resultDirectory,
            int parallelNumber,
            int threadNumber) {
        this.PRIMARY_ACCESSION_NUMBERS = primaryAccessionNumbers;
        this.RESULT_DIRECTORY = resultDirectory;
        this.PARALLEL_NUMBER = parallelNumber;
        this.SINGLE_TASK_THREAD = threadNumber / parallelNumber;
    }

    public void getUnmappedRead() {
        logger.info("获取 Unmmaped Reads");
        String shellScriptsDirectory = Paths.get(RESULT_DIRECTORY, ".script/unmapped_reads").toString();
        FileOperator.mkdir(shellScriptsDirectory);

        for (String accessionNumber : PRIMARY_ACCESSION_NUMBERS) {
            String unmappedReadsDirectory = Paths.get(RESULT_DIRECTORY, "align_fastq_to_reference", accessionNumber)
                    .toString();
            String bamFile = Paths.get(unmappedReadsDirectory, accessionNumber + ".bam").toString();
            String unmappedBamFile = Paths.get(unmappedReadsDirectory, accessionNumber + ".unmapped.bam").toString();
            String clippingBamFile = Paths.get(unmappedReadsDirectory, accessionNumber + ".clipping.bam").toString();
            String mergedBamFile = Paths.get(unmappedReadsDirectory, accessionNumber + ".unmapped_clipping.bam")
                    .toString();
            String unmappedReadFASTQ = Paths.get(unmappedReadsDirectory, accessionNumber + ".unmapped.fq").toString();
            String shellScript = Paths.get(shellScriptsDirectory, accessionNumber + ".sh").toString();

            List<String> commands = new ArrayList<>();
            commands.add(String.format("samtools view -@ %s -f 4 -b %s -o %s",
                    SINGLE_TASK_THREAD, bamFile, unmappedBamFile));
            commands.add("");
            commands.add(String.format("samtools view -@ %s -H %s > %s",
                    SINGLE_TASK_THREAD, bamFile, clippingBamFile));
            commands.add("");
            commands.add(String.format("samtools view -@ %s %s | awk '$6~/S/' >> %s",
                    SINGLE_TASK_THREAD, bamFile, clippingBamFile));
            commands.add("");
            commands.add(String.format("samtools merge -f -o %s %s %s",
                    mergedBamFile, unmappedBamFile, clippingBamFile));
            commands.add("");
            commands.add(String.format("samtools bam2fq %s > %s", mergedBamFile, unmappedReadFASTQ));

            TXTWriter.write(shellScript, commands);
        }

        BioinformaticToolsExecutor.getUnmmapedReadsFastq(PRIMARY_ACCESSION_NUMBERS, RESULT_DIRECTORY, PARALLEL_NUMBER);
    }
}
