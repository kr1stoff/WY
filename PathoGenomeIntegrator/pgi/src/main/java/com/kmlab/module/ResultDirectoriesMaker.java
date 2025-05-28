package com.kmlab.module;

import java.nio.file.Paths;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.FileOperator;

public class ResultDirectoriesMaker {
    private static final Logger logger = LogManager.getLogger(ResultDirectoriesMaker.class);

    private static final String[] BASIC_DIRS = { "fake_reads", "align_fastq_to_reference", "unmapped_reads_assembly",
            "align_fastq_to_integrated", "alignment_stats", ".script", ".reference", ".genomes", ".log" };

    private final String RESULT_DIRECTORY;

    public ResultDirectoriesMaker(String resultDir) {
        this.RESULT_DIRECTORY = resultDir;
    }

    /**
     * 创建结果目录及其子目录结构。
     * 该方法首先创建一个指定的结果目录，然后在该目录下创建一系列的子目录，用于组织和存储不同的结果文件。
     */
    public void makeResultDirectories() {
        logger.info("创建结果目录: " +String.join(" ", BASIC_DIRS));
        FileOperator.mkdir(RESULT_DIRECTORY);
        FileOperator.mkdirs(RESULT_DIRECTORY, BASIC_DIRS);
    }

    /**
     * 创建日志目录，并为每个提供的主要访问号建立子目录。
     *
     * @param primaryAccessionNumbers 包含主要访问号的列表，每个访问号将对应一个日志子目录。
     */
    public void makeLogDirectories(List<String> primaryAccessionNumbers) {
        logger.info("创建日志目录下基因组子目录");

        for (String accesion : primaryAccessionNumbers) {
            FileOperator.mkdir(Paths.get(RESULT_DIRECTORY, ".log", accesion).toString());
        }
    }

    public void makeScriptDirectories() {
        logger.info("创建脚本目录下子目录");

        final String[] scriptDirectory = { "bin" };

        for (String dir : scriptDirectory) {
            FileOperator.mkdir(Paths.get(RESULT_DIRECTORY, ".script", dir).toString());
        }
    }
}
