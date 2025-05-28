package com.kmlab.module;

import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

public class FastqIntegratedGenomeAlignerTest {
    @Test
    public void test() {
        final String RESULT_DIRECTORY = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426";
        final String GENOMES_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".genomes").toString();

        StartDataPreparer startDataPreparer = new StartDataPreparer(GENOMES_DIRECTORY, RESULT_DIRECTORY);
        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();

        FastqIntegratedGenomeAligner fastqIntegratedGenomeAligner = new FastqIntegratedGenomeAligner(
                primaryAccessionNumbers, RESULT_DIRECTORY, 8, 128);
        fastqIntegratedGenomeAligner.alignFastqToIntegratedGenome();
    }
}
