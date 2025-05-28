package com.kmlab.module;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.FileOperator;
import com.kmlab.util.TXTWriter;

public class StartDataPreparer {
    private static final Logger logger = LogManager.getLogger(StartDataPreparer.class);

    private static final String[] EXTENSIONS = { ".fna", ".fa", ".fasta" };
    private final String GENOMES_DIRECTORY;
    private final String RESULT_DIRECTORY;

    public StartDataPreparer(String dirGenomes, String dirResult) {
        this.GENOMES_DIRECTORY = dirGenomes;
        this.RESULT_DIRECTORY = dirResult;
    }

    /**
     * 查找基因组文件的方法。
     * 该方法不接受任何参数，它会在指定的基因组目录中搜索所有具有特定扩展名的文件。
     *
     * @return 返回一个包含找到的所有基因组文件的列表。列表中的每个元素都是一个File对象。
     */
    private List<File> findGenomeFiles() {
        // 在基因组目录中，根据指定的扩展名查找所有文件
        List<File> genomeFiles = FileOperator.findFilesWithExtensions(GENOMES_DIRECTORY, EXTENSIONS);
        if (genomeFiles.isEmpty()) {
            throw new RuntimeException("没有发现基因组文件. 基因组目录:" + GENOMES_DIRECTORY +
                    " 文件后缀: " + String.join(",", EXTENSIONS));
        }
        return genomeFiles;
    }

    /**
     * 从文件名中提取访问序号（accession）。
     * 此方法假设文件名包含一个点（.），该点之前的部分即为访问序号。
     *
     * @param file 表示要处理的文件对象。
     * @return 返回提取出的访问序号字符串。
     */
    private String getAccession(File file) {
        // 获取文件的完整名称，包括文件扩展名
        String fileNameWithExt = file.getName();
        // 查找文件名中点的位置，为分割文件名和扩展名做准备
        int dotIndex = fileNameWithExt.indexOf('.');
        // 截取点之前的部分作为访问序号
        String accession = fileNameWithExt.substring(0, dotIndex);
        return accession;
    }

    /**
     * 软连接基因组序列文件的方法。
     * 此方法遍历找到的所有基因组文件，并为每个文件创建一个软链接到指定的结果目录中。
     */
    public void linkGenomeSequences() {
        logger.info("软连接基因组序列文件");

        for (File file : findGenomeFiles()) {
            Path linkTargetPath = Paths.get(RESULT_DIRECTORY + "/.genomes/" + getAccession(file) + ".fna");
            FileOperator.linkFile(file.toPath(), linkTargetPath);
        }
    }

    /**
     * 获取原始基因组目录中所有文件的登录号列表。
     * 此方法不接受参数。
     *
     * @return 返回一个包含所有文件主要访问编号（登录号）的字符串列表。
     */
    public List<String> getPrimaryAccessionNumbers() {
        logger.info("获取原始基因组目录中登录号");

        List<String> primaryAccessionNumbers = new ArrayList<>();

        for (File file : findGenomeFiles()) {
            primaryAccessionNumbers.add(getAccession(file));
        }

        return primaryAccessionNumbers;
    }

    /**
     * 生成原始Accession列表文件。
     * 此方法不接受参数，也不返回任何值。
     * 获取主要Accession号码列表, 将这些Accession号码写入到一个文本文件中。
     */
    public void generatePrimaryAccessionFile() {
        logger.info("生成原始Accession列表文件");
        String accessionNumbersFile = Paths.get(RESULT_DIRECTORY, ".genomes/primary_accession_numbers.txt").toString();
        logger.debug("accessionNumbersFile: " + accessionNumbersFile);
        List<String> primaryAccessionNumbers = getPrimaryAccessionNumbers();
        TXTWriter.write(accessionNumbersFile, primaryAccessionNumbers);
    }
}
