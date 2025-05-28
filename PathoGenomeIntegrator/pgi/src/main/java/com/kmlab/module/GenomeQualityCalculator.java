package com.kmlab.module;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.kmlab.util.RankingCalculator;
import com.kmlab.util.TSVReader;
import com.kmlab.util.TXTWriter;

public class GenomeQualityCalculator {
    private static final Logger logger = LogManager.getLogger(GenomeQualityCalculator.class);

    private final String RESULT_DIRECTORY;
    private final String CHECKM_OUT_DIRECTORY;
    private final String PROKKA_OUTPUT_DIRECTORY;
    private final String COMP_OUTPUT_DIRECTORY;
    private final String GENOME_QUANLITY_FILE;
    private final String CHECKM_FILE;
    private final Map<String, String> PARAMETERS_PROPERTIY_MAP;

    public GenomeQualityCalculator(String resultDirectory, Map<String, String> pararmetersPropertyMap) {
        this.RESULT_DIRECTORY = resultDirectory;
        this.PARAMETERS_PROPERTIY_MAP = pararmetersPropertyMap;
        this.CHECKM_OUT_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".reference/checkm_out").toString();
        this.PROKKA_OUTPUT_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".reference/prokka_out").toString();
        this.COMP_OUTPUT_DIRECTORY = Paths.get(RESULT_DIRECTORY, ".reference/seqtk_comp_out").toString();
        this.GENOME_QUANLITY_FILE = Paths.get(RESULT_DIRECTORY, ".reference/genome_quality.tsv").toString();
        this.CHECKM_FILE = Paths.get(CHECKM_OUT_DIRECTORY, "checkm.tsv").toString();
    }

    /**
     * 计算基因组质量。该方法首先会汇总基因组信息，然后根据特定标准选出质量最高的基因组。
     *
     * @return 返回质量最高的基因组的访问码（Accession）。
     */
    public String calculateGenomeQuality() {
        logger.info("汇总计算基因组质量");
        Map<String, Map<String, Object>> accessionGenomeQualityMap = getAccessionGenomeQualityMap();
        List<String> genomeCountRank1Accessions = getGenomeCountRank1Accessions(accessionGenomeQualityMap);
        List<String> filteredByContigCountAccessionsv = filterAccessionsByContigCount(
                genomeCountRank1Accessions, accessionGenomeQualityMap);
        String referenceAccession = getReferenceAccession(accessionGenomeQualityMap, filteredByContigCountAccessionsv);

        return referenceAccession;
    }

    /**
     * 根据基因组的contig数量过滤接入序列列表。
     * 该方法遍历给定的接入序列列表，并使用提供的基因组质量映射来检查每个接入序列的contig数量。
     * 如果contig数量小于10，则认为该接入序列符合条件，并将其添加到过滤后的列表中。
     *
     * @param accessions                接入序列列表，是一组需要进行筛选的序列标识符。
     * @param accessionGenomeQualityMap 基因组质量映射，映射中包含了每个接入序列及其对应的基因组质量信息，其中包含了contig数量。
     * @return filteredByContigCountAccessions 过滤后的接入序列列表，仅包含contig数量小于10的序列。
     */
    public List<String> filterAccessionsByContigCount(
            List<String> accessions,
            Map<String, Map<String, Object>> accessionGenomeQualityMap) {
        int minimumContigCount = Integer.parseInt(PARAMETERS_PROPERTIY_MAP.get("minimun.contig.count"));

        List<String> filteredByContigCountAccessions = new ArrayList<>();

        for (String accession : accessions) {
            Map<String, Object> genomeQualityMap = accessionGenomeQualityMap.get(accession);
            int contigCount = (int) genomeQualityMap.get("contigCount");

            if (contigCount < minimumContigCount) {
                filteredByContigCountAccessions.add(accession);
            }
        }

        return filteredByContigCountAccessions;
    }

    /**
     * 从给定的accessions映射中，根据总分数确定参考accession。
     * 首先，此方法将根据提供的accession列表过滤映射，仅保留那些在列表中的accession的品质信息。
     * 然后，计算所有这些accession的“总分数”（totalScore）并找出最高的分数。
     * 最后，返回具有最高总分数的accession作为参考accession。
     *
     * @param accessionGenomeQualityMap 一个映射，其中包含accession到其基因组品质信息的映射。
     *                                  每个品质信息映射包含各种指标，例如“总分数”。
     * @param accessions                一个列表，包含根据基因数量排名的accession。
     *                                  该方法将仅考虑此列表中的accession。
     * @return 返回具有最高总分数的accession。如果没有提供accession或所有accession的总分数相同，则返回空字符串。
     */
    public String getReferenceAccession(
            Map<String, Map<String, Object>> accessionGenomeQualityMap,
            List<String> accessions) {
        Map<String, Map<String, Object>> genomeCountRank1AccessionGenomeQualityMap = new HashMap<>();

        for (String accession : accessions) {
            Map<String, Object> genomeQualityMap = accessionGenomeQualityMap.get(accession);
            genomeCountRank1AccessionGenomeQualityMap.put(accession, genomeQualityMap);
        }

        String referenceAccession = "";
        List<Float> totalScores = collectGenomeQualityMetricValues(genomeCountRank1AccessionGenomeQualityMap,
                "totalScore", Float.class);
        Collections.sort(totalScores);
        float maxScore = totalScores.get(totalScores.size() - 1);

        for (String accession : genomeCountRank1AccessionGenomeQualityMap.keySet()) {
            Map<String, Object> genomeQuality = genomeCountRank1AccessionGenomeQualityMap.get(accession);
            if ((float) genomeQuality.get("totalScore") == maxScore) {
                referenceAccession = accession;
            }
        }

        return referenceAccession;
    }

    /**
     * 获取基因组计数排名为1的accession列表。
     *
     * @param accessionGenomeQualityMap 包含accession和其基因组质量指标的映射的映射。其中，外部键是accession，内部映射的键是指标名称（如基因组计数），值是指标的值。
     * @return 一个包含所有基因组计数排名为1(lineage级别更低，属相对于科，更精准的分类)的accession的列表。
     */
    public List<String> getGenomeCountRank1Accessions(Map<String, Map<String, Object>> accessionGenomeQualityMap) {
        List<Integer> genomeCounts = collectGenomeQualityMetricValues(
                accessionGenomeQualityMap, "checkmGenomeCount", Integer.class);
        Map<Integer, Integer> genomeCountRank = RankingCalculator.intMapRank(genomeCounts);
        List<String> genomeCountRank1Accessions = new ArrayList<String>();

        for (String accession : accessionGenomeQualityMap.keySet()) {
            Map<String, Object> genomeQuality = accessionGenomeQualityMap.get(accession);
            int checkmGenomeCount = (int) genomeQuality.get("checkmGenomeCount");

            if (genomeCountRank.get(checkmGenomeCount) == 1) {
                genomeCountRank1Accessions.add(accession);
            }
        }

        return genomeCountRank1Accessions;
    }

    /**
     * 根据提供的checkm文件，解析并生成基因组质量映射表。
     *
     * @return 返回一个映射表，键为基因组的访问序号，值为对应的基因组质量对象，包含了完整性、污染度、假基因占比、N碱基占比和总分等信息。
     */
    public Map<String, Map<String, Object>> getAccessionGenomeQualityMap() {
        List<String> outputLines = new ArrayList<String>();
        outputLines.add(String.join("\t", "accessionNumber", "checkmGenomeCount",
                "completeness", "contamination", "pseudogenePercent", "nPercent", "contigCount", "totalScore"));

        Map<String, Map<String, Object>> accessionGenomeQualityMap = new HashMap<>();
        CSVParser csvParser = TSVReader.parseTSV(CHECKM_FILE, true);

        for (CSVRecord csvRecord : csvParser) {
            String accessionNumber = csvRecord.get("Bin Id");
            int checkmGenomeCount = Integer.parseInt(csvRecord.get("# genomes"));
            float completeness = Float.parseFloat(csvRecord.get("Completeness"));
            float contamination = Float.parseFloat(csvRecord.get("Contamination"));
            String prokkaTSV = Paths.get(PROKKA_OUTPUT_DIRECTORY, accessionNumber, accessionNumber + ".tsv").toString();
            float pseudogenePercent = getPseudogenePercent(prokkaTSV);
            String compFile = Paths.get(COMP_OUTPUT_DIRECTORY, accessionNumber + ".fna.comp").toString();
            float nPercent = getNBasePercent(compFile);
            int contigCount = getContigCount(compFile);
            float totalScore = completeness + (100 - contamination) + (100 - pseudogenePercent)
                    + (100 - nPercent);

            Map<String, Object> genomeQuality = new HashMap<>();
            genomeQuality.put("accessionNumber", accessionNumber);
            genomeQuality.put("checkmGenomeCount", checkmGenomeCount);
            genomeQuality.put("contigCount", contigCount);
            genomeQuality.put("totalScore", totalScore);

            outputLines.add(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    accessionNumber, checkmGenomeCount, completeness, contamination,
                    pseudogenePercent, nPercent, contigCount, totalScore));
            accessionGenomeQualityMap.put(accessionNumber, genomeQuality);
        }

        TXTWriter.write(GENOME_QUANLITY_FILE, outputLines);
        return accessionGenomeQualityMap;
    }

    /**
     * 收集基因组质量指标值。
     * 该方法遍历提供的访问控制序列表和其对应的基因组质量映射，从中提取特定指标的值，并将这些值转换为给定类型后存储在一个列表中。
     *
     * @param accessionGenomeQualityMap 一个映射，其键为访问控制序号，值为另一个映射，该映射的键为指标名称，值为指标的数值或其他对象。
     * @param metricKey                 要提取的指标的键名。
     * @param valueType                 指标值期望的类型，必须是Number或其子类。
     * @return 一个包含所有提取的指标值的列表，这些值已转换为指定的类型T。
     */
    public static <T extends Number> List<T> collectGenomeQualityMetricValues(
            Map<String, Map<String, Object>> accessionGenomeQualityMap,
            String metricKey, Class<T> valueType) {
        List<T> metricValues = new ArrayList<>();

        for (String accession : accessionGenomeQualityMap.keySet()) {
            Map<String, Object> genomeQualityMap = accessionGenomeQualityMap.get(accession);
            Number metricValue = (Number) genomeQualityMap.get(metricKey);
            if (metricValue != null) {
                metricValues.add(valueType.cast(metricValue));
            }
        }

        return metricValues;
    }

    /**
     * 计算给定文件中所有序列的N基平均百分比
     *
     * @param compFile 需要分析的文件路径，假设文件格式为TSV
     * @return N基在所有序列总长度中的平均百分比
     */
    public static float getNBasePercent(String compFile) {
        int totalLength = 0;
        int totalNbaseCount = 0;
        CSVParser csvParser = TSVReader.parseTSV(compFile);

        for (CSVRecord csvRecord : csvParser) {
            int seqLength = Integer.parseInt(csvRecord.get(1));
            int nBaseCount = Integer.parseInt(csvRecord.get(8));
            totalLength += seqLength;
            totalNbaseCount += nBaseCount;
        }

        return (float) totalNbaseCount / totalLength * 100;
    }

    /**
     * 获取给定组件文件中 contig 的数量。
     *
     * @param compFile 组件文件的路径，该文件被预期为TSV格式。
     * @return 返回文件中的 contig 数量。
     */
    public static int getContigCount(String compFile) {
        CSVParser csvParser = TSVReader.parseTSV(compFile);
        return csvParser.getRecords().size();
    }

    /**
     * 计算伪基因在基因组中的百分比
     *
     * @param prokkaTSV Prokka格式的TSV文件内容，用于分析基因和CDS的数量
     * @return 伪基因占基因总数的百分比
     */
    public static float getPseudogenePercent(String prokkaTSV) {
        int geneCount = 0;
        int CDSCount = 0;
        CSVParser csvParser = TSVReader.parseTSV(prokkaTSV, true);

        for (CSVRecord record : csvParser) {
            if (record.get("ftype").equals("gene")) {
                geneCount++;
            } else if (record.get("ftype").equals("CDS")) {
                CDSCount++;
            }
        }
        return ((float) (geneCount - CDSCount) / geneCount) * 100;
    }
}
