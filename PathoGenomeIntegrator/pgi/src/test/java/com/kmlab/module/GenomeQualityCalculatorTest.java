package com.kmlab.module;

import java.util.Map;

import org.junit.Test;

import com.kmlab.util.PropertyLoader;

public class GenomeQualityCalculatorTest {
    @Test
    public void testGenomeQualityCalculator() {
        Map<String, String> PARAMETERS_PROPERTIY_MAP = PropertyLoader.loadParametersProperties();

        GenomeQualityCalculator genomeQualityCalculator = new GenomeQualityCalculator(
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426", PARAMETERS_PROPERTIY_MAP);

        String referenceAccession = genomeQualityCalculator.calculateGenomeQuality();
        System.out.println(referenceAccession);
    }
}
