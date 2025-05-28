package com.kmlab.module;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.CondaEnvExecutor;
import com.kmlab.util.FileOperator;

public class BioinformaticToolsExecutor {
    private static final Logger logger = LogManager.getLogger(BioinformaticToolsExecutor.class);

    /**
     * 使用BWA工具创建索引
     *
     * @param fastaPath Fasta文件的路径，BWA将为这个文件创建索引
     * @param logFile   用于记录执行过程和结果的日志文件路径
     */
    public static void bwaIndex(String fastaPath, String logFile) {
        logger.info("创建 BWA 索引");
        String logFileFullPath = Paths.get(logFile).toString();
        String command = "bwa index " + fastaPath;
        CondaEnvExecutor.executeCommand(command, "bio", logFileFullPath);
    }

    /**
     * 使用CheckM工具对基因组进行质量评估。
     *
     * @param resultDirectory 结果目录的路径，CheckM的输出将保存在此目录的子目录".reference/checkm_out"中。
     * @param logDirectory    日志目录的路径，用于存储运行CheckM时产生的日志文件"find_reference_checkm.log"。
     * @param threadNumber    使用的线程数量，用于加速CheckM的执行。
     * @param genomeDirectory 基因组目录的路径，其中包含待检查的基因组文件。
     */
    public static void checkmGenomes(
            String resultDirectory,
            String logDirectory,
            int threadNumber,
            String genomeDirectory) {
        logger.info("运行 CheckM");

        String CHECKM_OUTPUT_DIRECTORY = Paths.get(resultDirectory, ".reference/checkm_out").toString();
        String FILE_LOG = Paths.get(logDirectory, "find_reference_checkm.log").toString();

        String command = String.format("checkm lineage_wf -t %s --tab_table -f %s/checkm.tsv -x fna %s %s",
                threadNumber, CHECKM_OUTPUT_DIRECTORY, genomeDirectory, CHECKM_OUTPUT_DIRECTORY);

        FileOperator.mkdir(CHECKM_OUTPUT_DIRECTORY);
        CondaEnvExecutor.executeCommand(command, "checkm", FILE_LOG);
    }

    // TODO 病毒使用checkV

    /**
     * 使用Prokka对基因组进行注释。
     *
     * @param resultDirectory        结果目录的路径，Prokka输出将存储在此目录的子目录".reference/prokka_out"中。
     * @param primaryAccessionNumber 基因组的主要访问号列表，每个基因组将根据此访问号进行标识和命名。
     * @param genomeDirectory        基因组文件目录的路径，其中每个基因组文件应以".fna"格式存储。
     * @param singleTaskThread       单个Prokka任务使用的线程数。
     * @param kingdom                为基因组指定的生物界，例如"bacteria"。
     * @param logDirectory           日志文件存储的目录路径，每个Prokka任务的日志将存储在此目录中。
     * @param parallelNumber         并行执行的Prokka任务数。
     */
    public static void prokkaGenomes(
            String resultDirectory,
            List<String> primaryAccessionNumber,
            String genomeDirectory,
            int singleTaskThread,
            String kingdom,
            String logDirectory,
            int parallelNumber) {
        logger.info("并行运行 prokka");

        String PROKKA_OUTPUT_DIRECTORY = Paths.get(resultDirectory, ".reference/prokka_out").toString();
        List<String[]> listCommandCondaLog = new ArrayList<>();

        FileOperator.mkdir(PROKKA_OUTPUT_DIRECTORY);

        for (String accessionNumber : primaryAccessionNumber) {
            String singleProkkaOut = Paths.get(PROKKA_OUTPUT_DIRECTORY, accessionNumber).toString();
            String fa = Paths.get(genomeDirectory, accessionNumber + ".fna").toString();

            String prokkaCommand = String.format(
                    "prokka %s --prefix %s --cpus %s --outdir %s --kingdom %s --force --addgenes",
                    fa, accessionNumber, singleTaskThread, singleProkkaOut, kingdom);
            String prokkaLog = Paths.get(logDirectory, accessionNumber, "find_reference_prokka.log").toString();
            String[] prokkaCommandCondaLog = { prokkaCommand, "assembly", prokkaLog };
            listCommandCondaLog.add(prokkaCommandCondaLog);
        }

        CondaEnvExecutor.parallelExecuteCommand(listCommandCondaLog, parallelNumber);
    }

    /**
     * 使用seqtk工具对给定的基因组进行压缩并行化处理。
     *
     * @param resultDirectory         结果目录的路径，处理后的压缩文件将保存在此目录的子目录".reference/seqtk_comp_out"中。
     * @param parallelNumber          并行处理的数量，决定同时执行的seqtk命令数。
     * @param primaryAccessionNumbers 基因组的主要访问号列表，每个访问号对应的基因组将被处理。
     * @param genomeDirectory         基因组文件所在的目录，每个基因组文件应以访问号+.fna命名。
     * @param logDirectory            日志文件保存的目录，每个处理过程的日志将保存在此目录下。
     */
    public static void seqtkCompGenomes(
            String resultDirectory,
            int parallelNumber,
            List<String> primaryAccessionNumbers,
            String genomeDirectory,
            String logDirectory) {
        logger.info("并行运行 seqtk comp");

        List<String[]> listCommandCondaLog = new ArrayList<>();
        String COMP_OUTPUT_DIRECTORY = Paths.get(resultDirectory, ".reference/seqtk_comp_out").toString();

        FileOperator.mkdir(COMP_OUTPUT_DIRECTORY);

        for (String accessionNumber : primaryAccessionNumbers) {
            String fa = Paths.get(genomeDirectory, accessionNumber + ".fna").toString();
            String command = String.format("seqtk comp %s > %s/%s.fna.comp",
                    fa, COMP_OUTPUT_DIRECTORY, accessionNumber);
            String logFile = Paths.get(logDirectory, accessionNumber, "find_reference_seqtk_comp.log").toString();
            String[] commandCondaLog = { command, "bio", logFile };
            listCommandCondaLog.add(commandCondaLog);
        }

        CondaEnvExecutor.parallelExecuteCommand(listCommandCondaLog, parallelNumber);
    }

    /**
     * 使用并行方式生成假的阅读数（fake reads）。
     *
     * @param primaryAccessionNumbers 主要的访问编号列表，用于指定需要生成假阅读数的基因组。
     * @param genomeDirectory         基因组文件所在的目录。
     * @param parametersPropertyMap   包含各种参数的属性映射，例如fake.reads.minimum.coverage。
     * @param resultDirectory         生成的假阅读数将保存在这个目录下。
     * @param logDirectory            日志文件将保存在这个目录下。
     * @param parallelNumber          并行执行命令的数量。
     */
    public static void snippyFakeReads(
            List<String> primaryAccessionNumbers,
            String genomeDirectory,
            Map<String, String> parametersPropertyMap,
            String resultDirectory,
            String logDirectory,
            int parallelNumber) {
        logger.info("并行运行生成 Fake Reads");
        List<String[]> listCommandCondaLog = new ArrayList<>();

        for (String accessionNumber : primaryAccessionNumbers) {
            String fa = Paths.get(genomeDirectory, accessionNumber + ".fna").toString();
            String outputFastq1 = Paths.get(resultDirectory, "fake_reads", accessionNumber + ".fq").toString();
            String minimumCoverage = parametersPropertyMap.get("fake.reads.minimum.coverage");
            String command = String.format(
                    "perl %s/.script/bin/snippy_fake_reads.pl --infa %s --outfq1 %s --mincov %s",
                    resultDirectory, fa, outputFastq1, minimumCoverage);
            String logFile = Paths.get(logDirectory, accessionNumber, "fake_reads.log").toString();
            String[] commandCondaLog = { command, "perl5.32", logFile };
            listCommandCondaLog.add(commandCondaLog);
        }

        CondaEnvExecutor.parallelExecuteCommand(listCommandCondaLog, parallelNumber);
    }

    /**
     * 将多个FASTQ文件对齐到参考基因组，使用BWA工具，并行执行。
     *
     * @param primaryAccessionNumbers 主要访问号列表，代表需要对齐的FASTQ样本。
     * @param resultDirectory         结果目录路径，所有结果和日志将保存在此目录下。
     * @param parallelNumber          并行执行的任务数，用于控制同时运行的BWA任务数量。
     */
    public static void bwaFastqToDatabaseX(
            List<String> primaryAccessionNumbers,
            String resultDirectory,
            int parallelNumber,
            String database) {
        logger.info("并行运行 BWA FASTQ to Reference");

        String baseName;
        if (database.equals("integrated")) {
            baseName = "align_fastq_to_integrated";
        } else if (database.equals("reference")) {
            baseName = "align_fastq_to_reference";
        } else {
            throw new RuntimeException("database must be integrated or reference");
        }

        List<String[]> listCommandCondaLog = new ArrayList<>();

        String shellScriptDirectory = Paths.get(resultDirectory, ".script", baseName).toString();
        String logDirectory = Paths.get(resultDirectory, ".log").toString();

        for (String accessionNumber : primaryAccessionNumbers) {
            String shellScript = Paths.get(shellScriptDirectory, accessionNumber + ".sh").toString();
            String command = "bash " + shellScript;
            String logFile = Paths.get(logDirectory, accessionNumber, baseName + ".log").toString();
            String[] commandCondaLog = { command, "bio", logFile };
            listCommandCondaLog.add(commandCondaLog);
        }

        CondaEnvExecutor.parallelExecuteCommand(listCommandCondaLog, parallelNumber);
    }

    /**
     * 并行获取未映射的Reads的FASTQ文件。
     *
     * @param primaryAccessionNumbers 主要访问编号列表，代表需要处理的样本。
     * @param resultDirectory         结果目录的路径，用于存储处理后的FASTQ文件和脚本。
     * @param parallelNumber          并行执行的任务数，用于控制同时执行的脚本数量。
     */
    public static void getUnmmapedReadsFastq(
            List<String> primaryAccessionNumbers,
            String resultDirectory,
            int parallelNumber) {
        logger.info("并行运行获取 Unmapped Reads FASTQ");
        List<String[]> listCommandCondaLog = new ArrayList<>();
        String shellScriptDirectory = Paths.get(resultDirectory, ".script/unmapped_reads").toString();
        String logDirectory = Paths.get(resultDirectory, ".log").toString();

        for (String accessionNumber : primaryAccessionNumbers) {
            String shellScript = Paths.get(shellScriptDirectory, accessionNumber + ".sh").toString();
            String command = "bash " + shellScript;
            String logFile = Paths.get(logDirectory, accessionNumber, "get_unmapped_reads.log").toString();
            String[] commandCondaLog = { command, "bio", logFile };
            listCommandCondaLog.add(commandCondaLog);
        }

        CondaEnvExecutor.parallelExecuteCommand(listCommandCondaLog, parallelNumber);
    }

    /**
     * 使用SPAdes基因组组装器进行基因组组装。
     *
     * @param resultDirectory 结果目录的路径，SPAdes的输出将存储在此目录下。
     * @param datasetsYaml    包含输入数据集信息的YAML文件路径。
     * @param threadNumber    用于组装过程的线程数量。
     */
    public static void runSPAdesGenomeAssembler(
            String resultDirectory,
            String datasetsYaml,
            int threadNumber) {
        logger.info("运行 SPAdes 组装");

        String spadesOutputDirectory = Paths.get(resultDirectory, "unmapped_reads_assembly/spades_out").toString();
        String logFile = Paths.get(resultDirectory, ".log/SPAdes.log").toString();

        String command = String.format("spades.py -t %s --phred-offset 33 -o %s --dataset %s",
                threadNumber, spadesOutputDirectory, datasetsYaml);

        CondaEnvExecutor.executeCommand(command, "assembly", logFile);
    }

    /**
     * 在给定的结果目录中运行统计脚本，以生成对齐信息的统计。
     *
     * @param resultDirectory 用于存储脚本结果的目录路径。
     */
    public static void alignmentStats(String resultDirectory) {
        logger.info("运行统计脚本 alignment_information_stats.py");
        String logFile = Paths.get(resultDirectory, ".log/alignment_stats.log").toString();

        String command = String.format(
                "cd %s && python .script/bin/alignment_information_stats.py -a .genomes/primary_accession_numbers.txt -A . -o alignment_stats",
                resultDirectory);

        CondaEnvExecutor.executeCommand(command, "python3.8", logFile);
    }
}
