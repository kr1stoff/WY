package com.kmlab;

import org.junit.Test;

/**
 * AppTest
 */
public class AppTest {
    @Test
    public void testCommand() {
        String[] args = { "-h" };
        // String[] args = {
        //         "-g", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes",
        //         "-o", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240520",
        //         "--no-genomes-screening" };
        App.main(args);
    }
}
