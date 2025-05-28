package com.kmlab.module;

import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.junit.Test;

import com.kmlab.util.PropertyLoader;

public class FakeReadGeneratorTest {
    @Test
    public void generateFakeReadsTest() {
        final String RESULT_DIRECTORY = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426";
        final String GENOMES_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".genomes").toString();
        final Map<String, String> PARAMETERS_PROPERTIY_MAP = PropertyLoader.loadParametersProperties();
        final int PARALLEL_NUMBER = 8;

        StartDataPreparer startDataPreparer = new StartDataPreparer(GENOMES_DIRECTORY, RESULT_DIRECTORY);
        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();

        FakeReadGenerator fakeReadGenerator = new FakeReadGenerator(RESULT_DIRECTORY, PARAMETERS_PROPERTIY_MAP,
                PARALLEL_NUMBER, primaryAccessionNumbers);

        fakeReadGenerator.generateFakeReads();
    }
}
