package com.kmlab.util;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TSVReader {
    private static final Logger logger = LogManager.getLogger(TSVReader.class);

    /**
     * 解析TSV文件为CSVParser对象。
     *
     * @param inFile 指定要解析的TSV文件的路径。
     * @return 返回一个配置为使用制表符作为分隔符，并且会裁剪字段的CSVParser对象。
     */
    public static CSVParser parseTSV(String inFile) {
        logger.info("读 TSV 文件: " + inFile);
        try {
            Reader reader = new FileReader(inFile);
            CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT.builder()
                    .setDelimiter('\t')
                    .setTrim(true)
                    .build());
            return csvParser;
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("读取TSV文件失败: " + e.getMessage());
        }
    }

    /**
     * 解析TSV文件的函数。
     *
     * @param inFile    指定要解析的TSV文件的路径。
     * @param hasHeader 指示TSV文件是否包含头部行。如果为true，则解析时会将第一行作为列名。
     * @return 返回一个CSVParser对象，用于进一步解析TSV文件的内容。
     */
    public static CSVParser parseTSV(String inFile, boolean hasHeader) {
        logger.info("读 TSV 文件: " + inFile);
        try {
            Reader reader = new FileReader(inFile);
            CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT.builder()
                    .setDelimiter('\t')
                    .setHeader()
                    .setSkipHeaderRecord(hasHeader)
                    .setIgnoreHeaderCase(true)
                    .setTrim(true)
                    .build());
            return csvParser;
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("读取TSV文件失败: " + e.getMessage());
        }
    }
}
