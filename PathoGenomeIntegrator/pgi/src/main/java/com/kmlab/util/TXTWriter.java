package com.kmlab.util;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TXTWriter {
    private static final Logger logger = LogManager.getLogger(TXTWriter.class);

    /**
     * 将字符串列表写入指定的文件。
     *
     * @param filePath 要写入的文件的路径。
     * @param lines    要写入文件的字符串列表。
     *                 注意：此方法不会在每行末尾添加换行符，每个字符串参数将作为单独的行写入文件。
     *
     * @return 无返回值。
     */
    public static void write(String filePath, List<String> lines) {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filePath))) {
            for (String line : lines) {
                writer.println(line);
            }
        } catch (IOException e) {
            logger.warn("写入文件失败: " + filePath);
        }
    }
}
