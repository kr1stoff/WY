package com.kmlab.module;

import org.junit.Test;

public class AlignmentInformationStatsTest {
    @Test
    public void test() {
        String resultDirectory = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240426";
        AlignmentInformationStats alignmentInformationStats = new AlignmentInformationStats(resultDirectory);
        alignmentInformationStats.stats();
    }
}
