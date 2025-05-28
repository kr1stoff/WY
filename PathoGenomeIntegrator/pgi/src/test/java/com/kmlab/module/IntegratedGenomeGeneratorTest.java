package com.kmlab.module;

import java.nio.file.Paths;

import org.junit.Test;

public class IntegratedGenomeGeneratorTest {
    @Test
    public void test() {
        final String RESULT_DIRECTORY = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426";
        IntegratedGenomeGenerator integratedGenomeGenerator = new IntegratedGenomeGenerator(RESULT_DIRECTORY);
        String integrateGenomeFasta = integratedGenomeGenerator.concatenateIntegratedGenome();

        String logFile = Paths.get(RESULT_DIRECTORY, ".log/bwa_index_integrated_genome.log").toString();
        BioinformaticToolsExecutor.bwaIndex(integrateGenomeFasta, logFile);
    }
}
