package com.kmlab.module;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.FileOperator;
import com.kmlab.util.TXTWriter;

public class FastqReferenceAligner {
    private static final Logger logger = LogManager.getLogger(FastqReferenceAligner.class);

    private final List<String> PRIMARY_ACCESSION_NUMBERS;
    private final String RESULT_DIRECTORY;
    private final int PARALLEL_NUMBER;
    private final int SINGLE_TASK_THREAD;

    public FastqReferenceAligner(
            List<String> primaryAccessionNumbers,
            String resultDirectory,
            int parallelNumber,
            int threadNumber) {
        this.PRIMARY_ACCESSION_NUMBERS = primaryAccessionNumbers;
        this.RESULT_DIRECTORY = resultDirectory;
        this.PARALLEL_NUMBER = parallelNumber;
        this.SINGLE_TASK_THREAD = threadNumber / parallelNumber;
    }

    public void alignFastqToReference() {
        logger.info("比对 Fake Reads 到参考基因组");

        String referencePath = Paths.get(RESULT_DIRECTORY, ".reference/ref.fna").toString();
        String fakeReadsDirectory = Paths.get(RESULT_DIRECTORY, "fake_reads").toString();
        String shellScriptsDirectory = Paths.get(RESULT_DIRECTORY, ".script/align_fastq_to_reference").toString();

        FileOperator.mkdir(shellScriptsDirectory);

        for (String accessionNumber : PRIMARY_ACCESSION_NUMBERS) {
            String alignFastqToReferenceDirectory = Paths
                    .get(RESULT_DIRECTORY, "align_fastq_to_reference", accessionNumber)
                    .toString();
            String fakeFastq = Paths.get(fakeReadsDirectory, accessionNumber + ".fq").toString();
            String sortedBam = Paths.get(alignFastqToReferenceDirectory, accessionNumber + ".bam").toString();
            String shellScript = Paths.get(shellScriptsDirectory, accessionNumber + ".sh").toString();

            FileOperator.mkdir(alignFastqToReferenceDirectory);

            List<String> commands = new ArrayList<>();
            commands.add(String.format(
                    "bwa mem -t %s -Y -M %s %s | samtools view  -@ %s -bS - | samtools sort -@ %s -o %s",
                    SINGLE_TASK_THREAD, referencePath, fakeFastq,
                    SINGLE_TASK_THREAD, SINGLE_TASK_THREAD, sortedBam));
            commands.add(""); // 加个空行, 方便看一点
            commands.add(String.format("samtools stats %s | grep ^SN | cut -f 2- > %s.stats",
                    sortedBam, sortedBam));

            TXTWriter.write(shellScript, commands);
        }

        BioinformaticToolsExecutor.bwaFastqToDatabaseX(PRIMARY_ACCESSION_NUMBERS, RESULT_DIRECTORY, PARALLEL_NUMBER,
                "reference");
    }
}
