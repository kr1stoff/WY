package com.kmlab.util;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ScriptGenerator {
    private static final Logger logger = LogManager.getLogger(RankingCalculator.class);

    /**
     * 从资源文件路径生成脚本到目标文件路径。
     *
     * @param resourceFilePath 资源文件路径，相对于类路径的路径
     * @param targetFilePath   目标文件路径，完整文件路径
     * @throws RuntimeException 如果资源文件不存在或复制文件时发生IO异常
     */
    public static void generateScript(String resourceFilePath, String targetFilePath) {
        logger.info("复制 resource 下脚本到指定位置");

        InputStream scriptStream = ScriptGenerator.class.getResourceAsStream(resourceFilePath);
        Path scriptPath = Paths.get(targetFilePath);

        if (scriptStream == null) {
            throw new RuntimeException("没找到 " + resourceFilePath + " 文件");
        }

        try {
            Files.copy(scriptStream, scriptPath, StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            throw new RuntimeException("复制 " + targetFilePath + " 文件失败", e);
        }
    }
}
