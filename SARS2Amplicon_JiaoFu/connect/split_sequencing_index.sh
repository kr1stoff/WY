## bcl2fastq
## 20210903-DIGUANGZHOU001-MN01481-PE150test
export RUN_FOLDER="/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/20210903-DIGUANGZHOU001-MN01481-PE150test/210902_MN01481_0035_A000H3KCTJ"
export FQ_OUTDIR="/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/bcl_to_fastq/20210903-DIGUANGZHOU001-MN01481-PE150test"
/usr/local/bin/bcl2fastq \
    --barcode-mismatches 0 \
    --runfolder-dir ${RUN_FOLDER} \
    --input-dir ${RUN_FOLDER}/Data/Intensities/BaseCalls \
    --output-dir ${FQ_OUTDIR}/fastq \
    --sample-sheet ${FQ_OUTDIR}/test1210_20210903-DIGUANGZHOU001-MN01481-PE150test.csv \
    --no-lane-splitting \
    > ${FQ_OUTDIR}/bcl2fastq.o \
    2> ${FQ_OUTDIR}/bcl2fastq.e

## 20210908-DIGUANGZHOU1032-TPNB500564-YF
export RUN_FOLDER="/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/20210908-DIGUANGZHOU1032-TPNB500564-YF/210908_TPNB500564_0018_AHCV3LBGXK"
export FQ_OUTDIR="/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/bcl_to_fastq/20210908-DIGUANGZHOU1032-TPNB500564-YF"
/usr/local/bin/bcl2fastq \
    --barcode-mismatches 0 \
    --runfolder-dir ${RUN_FOLDER} \
    --input-dir ${RUN_FOLDER}/Data/Intensities/BaseCalls \
    --output-dir ${FQ_OUTDIR}/fastq \
    --sample-sheet ${FQ_OUTDIR}/20210908-DIGUANGZHOU1032-TPNB500564-YF_sample_sheet.csv \
    --no-lane-splitting \
    > ${FQ_OUTDIR}/bcl2fastq.o \
    2> ${FQ_OUTDIR}/bcl2fastq.e

## liuzu split barcode
/usr/local/bin/bcl2fastq --barcode-mismatches 0 \
    --no-lane-splitting \
    --sample-sheet $SAMPLE_SHEET \
    -R $SEQ_DIR \
    -i $SEQ_DIR/Data/Intensities/BaseCalls \
    --output-dir $FQ_OUTDIR \
    > $FQ_OUTDIR/bcl2fastq.o \
    2> $FQ_OUTDIR/bcl2fastq.e