package com.kmlab.module;

import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

public class SPAdesGenomeAssemblerTest {
    @Test
    public void main() {
        final String RESULT_DIRECTORY = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426";
        final String GENOMES_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".genomes").toString();
        final int THREAD_NUMBER = 128;

        StartDataPreparer startDataPreparer = new StartDataPreparer(GENOMES_DIRECTORY, RESULT_DIRECTORY);
        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();

        SPAdesGenomeAssembler spadesGenomeAssembler = new SPAdesGenomeAssembler(primaryAccessionNumbers,
                RESULT_DIRECTORY, THREAD_NUMBER);
        spadesGenomeAssembler.runSPAdes();
    }
}
