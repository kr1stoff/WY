package com.kmlab.module;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.FileOperator;
import com.kmlab.util.TXTWriter;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class FastqIntegratedGenomeAligner {
    private static final Logger logger = LogManager.getLogger(FastqIntegratedGenomeAligner.class);

    private final List<String> PRIMARY_ACCESSION_NUMBERS;
    private final String RESULT_DIRECTORY;
    private final int PARALLEL_NUMBER;
    private final int SINGLE_TASK_THREAD;

    public FastqIntegratedGenomeAligner(
            List<String> primaryAccessionNumbers,
            String resultDirectory,
            int parallelNumber,
            int threadNumber) {
        this.PRIMARY_ACCESSION_NUMBERS = primaryAccessionNumbers;
        this.RESULT_DIRECTORY = resultDirectory;
        this.PARALLEL_NUMBER = parallelNumber;
        this.SINGLE_TASK_THREAD = threadNumber / parallelNumber;
    }

    public void alignFastqToIntegratedGenome() {
        logger.info("比对 Fake Reads 到整合基因组");
        String fakeReadsDirectory = Paths.get(RESULT_DIRECTORY, "fake_reads").toString();
        String shellScriptsDirectory = Paths.get(RESULT_DIRECTORY, ".script/align_fastq_to_integrated").toString();
        String integrateGenomeFasta = Paths.get(RESULT_DIRECTORY, "unmapped_reads_assembly/integrate_genome.fna")
                .toString();

        FileOperator.mkdir(shellScriptsDirectory);

        for (String accessionNumber : PRIMARY_ACCESSION_NUMBERS) {
            String alignToIntegratedDirectory = Paths
                    .get(RESULT_DIRECTORY, "align_fastq_to_integrated", accessionNumber).toString();
            String fakeFastq = Paths.get(fakeReadsDirectory, accessionNumber + ".fq").toString();
            String sortedBam = Paths.get(alignToIntegratedDirectory, accessionNumber + ".bam").toString();
            String shellScript = Paths.get(shellScriptsDirectory, accessionNumber + ".sh").toString();

            FileOperator.mkdir(alignToIntegratedDirectory);

            List<String> commands = new ArrayList<>();
            commands.add(String.format(
                    "bwa mem -t %s -Y -M %s %s | samtools view -@ %s -bS - | samtools sort -@ 8 -o %s -",
                    SINGLE_TASK_THREAD, integrateGenomeFasta, fakeFastq, SINGLE_TASK_THREAD, sortedBam));
            commands.add("");
            commands.add(String.format("samtools stats %s | grep ^SN | cut -f 2- > %s.stats", sortedBam, sortedBam));

            TXTWriter.write(shellScript, commands);
        }

        BioinformaticToolsExecutor.bwaFastqToDatabaseX(PRIMARY_ACCESSION_NUMBERS, RESULT_DIRECTORY, PARALLEL_NUMBER,
                "integrated");
    }
}
