package com.kmlab.cli;

import java.nio.file.Paths;
import java.nio.file.Files;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class OptionsValidator {
    private static final Logger logger = LogManager.getLogger(OptionsValidator.class);

    /**
     * 验证基因组目录的有效性。
     *
     * @param dirGeno 基因组目录的路径。必须是一个存在的目录。
     * @throws RuntimeException 如果基因组目录不存在，抛出此异常。
     */
    public static void validateGenomeDirAndRefFile(String dirGeno) {
        logger.info("验证基因组目录的有效性");

        if (!Files.exists(Paths.get(dirGeno))) {
            throw new RuntimeException("基因组目录不存在: " + dirGeno);
        }
    }
}
