package com.kmlab.module;

import org.junit.Test;

public class StartDataPreparerTest {
    @Test
    public void main() {
        StartDataPreparer startDataPreparer = new StartDataPreparer(
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes",
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426");

        startDataPreparer.linkGenomeSequences();
        // startDataPreparer.generatePrimaryAccessionFile();
    }
}
