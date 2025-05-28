package com.kmlab.module;

import java.util.List;
import java.util.Map;

import org.junit.Test;

import com.kmlab.util.PropertyLoader;

public class ReferenceSequenceProcessorTest {
    @Test
    public void testAssgin() {
        Map<String, String> PARAMETERS_PROPERTIY_MAP = PropertyLoader.loadParametersProperties();

        StartDataPreparer startDataPreparer = new StartDataPreparer(
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes",
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426");

        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();

        ReferenceSequenceProcessor referenceSequenceProcessor = new ReferenceSequenceProcessor(
                "", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426", 128, "bacteria", 8,
                primaryAccessionNumbers, PARAMETERS_PROPERTIY_MAP);

        referenceSequenceProcessor.assignReferenceAccession();
    }
}
