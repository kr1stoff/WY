package com.kmlab.util;

import java.io.FileWriter;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.yaml.snakeyaml.DumperOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.representer.Representer;

public class YamlWriter {
    private static final Logger logger = LogManager.getLogger(YamlWriter.class);

    /**
     * 将对象内容写入到指定的YAML文件中。
     *
     * @param content  要写入的对象。这个对象会被转换成YAML格式。
     * @param filePath 要写入的文件路径。如果文件不存在，将会被创建。
     */
    public static void write(Object content, String filePath) {
        logger.info("写 YAML 文件: " + filePath);
        DumperOptions options = new DumperOptions();
        options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK);

        Yaml yaml = new Yaml(new Representer(options));

        try (FileWriter writer = new FileWriter(filePath)) {
            yaml.dump(content, writer);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
