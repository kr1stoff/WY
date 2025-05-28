package com.kmlab.module;

import java.util.List;

import org.junit.Test;

/**
 * ResultDirectoriesMakerTest
 */
public class ResultDirectoriesMakerTest {

    @Test
    public void main() {
        StartDataPreparer startDataPreparer = new StartDataPreparer(
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes",
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426");

        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();

        ResultDirectoriesMaker resultDirectoriesMaker = new ResultDirectoriesMaker(
                "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426");
        resultDirectoriesMaker.makeResultDirectories();
        resultDirectoriesMaker.makeLogDirectories(primaryAccessionNumbers);
        resultDirectoriesMaker.makeScriptDirectories();
    }
}
