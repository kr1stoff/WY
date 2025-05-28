package com.kmlab.pipe;

import java.util.Map;
import java.util.Optional;

import org.junit.Test;

import com.kmlab.cli.OptionsToHashMap;

public class NoGenomesScreeningPipeTest {
    @Test
    public void test() {
        String[] args = {
                "-g", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/data/txid470_small_genomes",
                "-o", "/sdbb/bioinfor/mengxf/TASKS/WY24012501/result/pgi_results_240520",
                "--no-genomes-screening"
        };
        Optional<Map<String, Object>> mapOptions = OptionsToHashMap.map(args);

        NoGenomesScreeningPipe noGenomesScreeningPipe = new NoGenomesScreeningPipe(mapOptions);
        noGenomesScreeningPipe.run();
    }
}
