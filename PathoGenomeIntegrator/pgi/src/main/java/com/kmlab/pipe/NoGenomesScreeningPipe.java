package com.kmlab.pipe;

import java.util.Optional;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.module.AlignmentInformationStats;
import com.kmlab.module.BioinformaticToolsExecutor;
import com.kmlab.module.FakeReadGenerator;
import com.kmlab.module.FastqIntegratedGenomeAligner;
import com.kmlab.module.FastqReferenceAligner;
import com.kmlab.module.IntegratedGenomeGenerator;
import com.kmlab.module.ReferenceSequenceProcessor;
import com.kmlab.module.ResultDirectoriesMaker;
import com.kmlab.module.SPAdesGenomeAssembler;
import com.kmlab.module.StartDataPreparer;
import com.kmlab.module.UnmappedReadProcessor;
import com.kmlab.util.PropertyLoader;

public class NoGenomesScreeningPipe {
    private static final Logger logger = LogManager.getLogger(NoGenomesScreeningPipe.class);

    private final String RESULT_DIRECTORY;
    private final String GENOMES_DIRECTORY;
    private final String REFERENCE_FILE;
    private final int THREAD_NUMBER;
    private final String KINGDOM;
    private final int PARALLEL_NUMBER;
    private final Map<String, String> PARAMETERS_PROPERTIY_MAP;

    public NoGenomesScreeningPipe(Optional<Map<String, Object>> mapOptions) {
        this.RESULT_DIRECTORY = mapOptions.get().get("dirResult").toString();
        this.GENOMES_DIRECTORY = mapOptions.get().get("dirGenomes").toString();
        this.REFERENCE_FILE = mapOptions.get().get("fileReference").toString();
        this.THREAD_NUMBER = Integer.parseInt(mapOptions.get().get("threads").toString());
        this.KINGDOM = mapOptions.get().get("kingdom").toString();
        Map<String, String> basicProperty = PropertyLoader.loadBasicProperties();
        this.PARALLEL_NUMBER = Integer.parseInt(basicProperty.get("parallel.number"));
        this.PARAMETERS_PROPERTIY_MAP = PropertyLoader.loadParametersProperties();
    }

    public void run() {
        logger.info("不做任何基因组筛选构建整合基因组");
        // 创建结果目录
        ResultDirectoriesMaker resultDirectoriesMaker = new ResultDirectoriesMaker(RESULT_DIRECTORY);
        resultDirectoriesMaker.makeResultDirectories();

        // 1.软连接菌株基因组
        // 2.获取原始基因组目录登录号列表
        StartDataPreparer startDataPreparer = new StartDataPreparer(GENOMES_DIRECTORY, RESULT_DIRECTORY);
        List<String> primaryAccessionNumbers = startDataPreparer.getPrimaryAccessionNumbers();
        startDataPreparer.generatePrimaryAccessionFile();
        startDataPreparer.linkGenomeSequences();

        // 创建日志,脚本目录
        resultDirectoriesMaker.makeLogDirectories(primaryAccessionNumbers);
        resultDirectoriesMaker.makeScriptDirectories();

        // 指定参考基因组并建索引
        ReferenceSequenceProcessor referenceSequenceProcessor = new ReferenceSequenceProcessor(
                REFERENCE_FILE, RESULT_DIRECTORY, THREAD_NUMBER, KINGDOM, PARALLEL_NUMBER,
                primaryAccessionNumbers, PARAMETERS_PROPERTIY_MAP);
        referenceSequenceProcessor.assignReferenceAccession();

        // 生成 Fake Reads
        FakeReadGenerator fakeReadGenerator = new FakeReadGenerator(RESULT_DIRECTORY,
                PARAMETERS_PROPERTIY_MAP,
                PARALLEL_NUMBER, primaryAccessionNumbers);
        fakeReadGenerator.generateFakeReads();

        // 比对 Fake Reads 到参考基因组
        FastqReferenceAligner fastqReferenceAligner = new FastqReferenceAligner(primaryAccessionNumbers,
                RESULT_DIRECTORY, PARALLEL_NUMBER, THREAD_NUMBER);
        fastqReferenceAligner.alignFastqToReference();

        // 获取 Unmapped Reads
        UnmappedReadProcessor unmappedReadProcessor = new UnmappedReadProcessor(
                primaryAccessionNumbers,
                RESULT_DIRECTORY,
                PARALLEL_NUMBER,
                THREAD_NUMBER);
        unmappedReadProcessor.getUnmappedRead();

        // SPAdes 组装
        SPAdesGenomeAssembler spadesGenomeAssembler = new SPAdesGenomeAssembler(primaryAccessionNumbers,
                RESULT_DIRECTORY, THREAD_NUMBER);
        spadesGenomeAssembler.runSPAdes();

        // 生成整合基因组
        IntegratedGenomeGenerator integratedGenomeGenerator = new IntegratedGenomeGenerator(RESULT_DIRECTORY);
        String integrateGenomeFasta = integratedGenomeGenerator.concatenateIntegratedGenome();
        String logFile = Paths.get(RESULT_DIRECTORY, ".log/bwa_index_integrated_genome.log").toString();
        BioinformaticToolsExecutor.bwaIndex(integrateGenomeFasta, logFile);

        // 比对 Fake Reads 到整合基因组
        FastqIntegratedGenomeAligner fastqIntegratedGenomeAligner = new FastqIntegratedGenomeAligner(
                primaryAccessionNumbers, RESULT_DIRECTORY, PARALLEL_NUMBER, THREAD_NUMBER);
        fastqIntegratedGenomeAligner.alignFastqToIntegratedGenome();

        // 生成比对统计文件
        AlignmentInformationStats alignmentInformationStats = new AlignmentInformationStats(RESULT_DIRECTORY);
        alignmentInformationStats.stats();
    }
}
