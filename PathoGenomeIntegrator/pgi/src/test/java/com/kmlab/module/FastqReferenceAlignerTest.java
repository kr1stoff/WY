package com.kmlab.module;

import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

public class FastqReferenceAlignerTest {
    @Test
    public void main() {
        final String RESULT_DIRECTORY = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426";
        final String GENOMES_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".genomes").toString();
        final int PARALLEL_NUMBER = 8;
        final int THREAD_NUMBER = 128;

        StartDataPreparer startDataPreparer = new StartDataPreparer(GENOMES_DIRECTORY, RESULT_DIRECTORY);
        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();

        FastqReferenceAligner fastqReferenceAligner = new FastqReferenceAligner(primaryAccessionNumbers,
                RESULT_DIRECTORY, PARALLEL_NUMBER, THREAD_NUMBER);
        fastqReferenceAligner.alignFastqToReference();
    }
}
