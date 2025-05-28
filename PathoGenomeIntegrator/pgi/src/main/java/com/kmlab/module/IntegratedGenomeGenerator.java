package com.kmlab.module;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;

public class IntegratedGenomeGenerator {
    private static final Logger logger = LogManager.getLogger(FastqIntegratedGenomeAligner.class);
    private final String RESULT_DIRECTORY;

    public IntegratedGenomeGenerator(String resultDirectory) {
        RESULT_DIRECTORY = resultDirectory;
    }

    /**
     * 将多个DNA序列文件合并为一个集成的基因组序列文件。
     * 此方法首先读取指定的参考基因组和组装的基因组序列，然后将它们合并到一个单一的FASTA文件中。
     *
     * @return 综合基因组的文件路径。
     */
    public String concatenateIntegratedGenome() {
        logger.info("合并参考和组装基因组生成整合基因组");
        String integrateGenomeFasta = Paths.get(RESULT_DIRECTORY, "unmapped_reads_assembly/integrate_genome.fna")
                .toString();
        String[] inputFiles = {
                Paths.get(RESULT_DIRECTORY, ".reference/ref.fna").toString(),
                Paths.get(RESULT_DIRECTORY, "unmapped_reads_assembly/spades_out/scaffolds.fasta").toString(),
        };
        Collection<DNASequence> allSequences = new ArrayList<>();

        for (String inputFile : inputFiles) {
            Map<String, DNASequence> sequences = readAmbiguityDNAFileToMap(inputFile);
            allSequences.addAll(sequences.values());
        }

        writeNucleotideSequenceToFile(new File(integrateGenomeFasta), allSequences);

        return integrateGenomeFasta;
    }

    /**
     * 将核酸序列集合写入指定的文件中。
     *
     * @param inputFile 指定的文件对象，核酸序列将被写入这个文件。
     * @param sequences 核酸序列集合，需要写入文件的多个DNASequence对象。
     * @throws RuntimeException 如果写入核酸序列过程中发生异常，则抛出运行时异常。
     */
    public void writeNucleotideSequenceToFile(File inputFile, Collection<DNASequence> sequences) {
        try {
            FastaWriterHelper.writeNucleotideSequence(inputFile, sequences);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("写入核酸序列失败", e);
        }
    }

    /**
     * 从指定的FASTA文件中读取DNA序列，并将它们存储在一个Map中，其中键为序列ID，值为对应的DNASequence对象。
     * 这个方法特别处理了含有不确定核苷酸（ambiguity nucleotides）的DNA序列。
     *
     * @param inputFile 必须是合法的FASTA文件路径，文件中包含了需要读取的DNA序列。
     * @return 一个Map，键为序列ID，值为对应的DNASequence对象，包含了从输入文件中读取的所有DNA序列。
     * @throws RuntimeException 如果读取文件或者处理序列时发生错误，会抛出此异常。
     */
    public Map<String, DNASequence> readAmbiguityDNAFileToMap(String inputFile) {
        try {
            Map<String, DNASequence> outputSequencesMap = new HashMap<>();
            AmbiguityDNACompoundSet ambiguityDNACompoundSet = AmbiguityDNACompoundSet.getDNACompoundSet();
            Map<String, DNASequence> sequences = FastaReaderHelper.readFastaDNASequence(new File(inputFile));

            for (String sequenceID : sequences.keySet()) {
                // TODO 长度过滤在这里 dnaSequence.getLength()
                String sequenceString = sequences.get(sequenceID).getSequenceAsString();
                DNASequence dnaSequence = new DNASequence(sequenceString, ambiguityDNACompoundSet);
                dnaSequence.setOriginalHeader(sequenceID);
                // logger.info(String.format("sequenceID: %s | dnaSequenceLength: %s",
                //         sequenceID, dnaSequence.getLength()));
                outputSequencesMap.put(sequenceID, dnaSequence);
            }

            return outputSequencesMap;
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("合并参考和 Unmmaped 组装序列失败", e);
        } catch (CompoundNotFoundException e) {
            e.printStackTrace();
            throw new RuntimeException("未找到核苷酸化合物对象", e);
        }
    }
}
