package com.kmlab.cli;

import org.junit.Test;

public class OptionsToHashMapTest {

    @Test
    public void testMap() {
        String[] args = {
                "-k", "virus",
                "-g", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes",
                "-o", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426",
                "--no-genomes-screening"
        };
        OptionsToHashMap.map(args);
    }
}
