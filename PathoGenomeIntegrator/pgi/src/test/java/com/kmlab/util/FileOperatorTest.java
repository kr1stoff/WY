package com.kmlab.util;

import org.junit.Test;

public class FileOperatorTest {
    @Test
    public void testFindFileByContentAndExtensions() {
        String[] extensions = { "fna", "fa", "fasta" };
        String dirGenomes = "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes";
        String targetContent = "GCF_008632635.1";
        String res = FileOperator.findFileByContentAndExtensions(dirGenomes, targetContent, extensions);
        if (res == "") {
            System.out.println("\"\"");
        } else {
            System.out.println(res);
        }
    }
}
