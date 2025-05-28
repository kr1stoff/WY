package com.kmlab;

import java.util.Optional;
import java.util.Map;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import com.kmlab.cli.OptionsToHashMap;
import com.kmlab.pipe.NoGenomesScreeningPipe;

public class App {

    private static final Logger logger = LogManager.getLogger(App.class);

    public static void main(String[] args) {
        logger.info("开始分析");
        Optional<Map<String, Object>> mapOptions = OptionsToHashMap.map(args);

        if (mapOptions.isPresent()) {
            if ((boolean) mapOptions.get().get("boolNoGenomesScreening")) {
                NoGenomesScreeningPipe noGenomesScreeningPipe = new NoGenomesScreeningPipe(mapOptions);
                noGenomesScreeningPipe.run();
            }
        }
    }
}
