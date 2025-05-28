include { VFDB_BLAST;VFDB_PARSE;VFDB_FASTA } from '../modules/vfdb.nf'
include { RGI;RGI_PARSE1;RGI_PARSE2 } from '../modules/card.nf'
include { RESFINDER;RESFINDER_PARSE } from '../modules/resfinder.nf'
include { EMAPPER;EMAPPER_PARSE } from '../modules/eggnog.nf'
include { CAZY;CAZY_PARSE1;CAZY_PARSE2 } from '../modules/cazy.nf'
include { SWISSPROT } from '../modules/swissprot.nf'
include { PLASMID_BLAST;PLASMID_FILTER } from '../modules/plasmid.nf'


workflow ANNOTATION {
  take:
  ch_dir_pred
  ch_depth_pred

  main:
  VFDB_BLAST(ch_dir_pred)
  VFDB_PARSE(VFDB_BLAST.out, ch_depth_pred)
  VFDB_FASTA(VFDB_PARSE.out, ch_dir_pred)
  RGI(ch_dir_pred)
  RGI_PARSE1(ch_dir_pred, RGI.out)
  RGI_PARSE2(ch_depth_pred, RGI_PARSE1.out[0])
  RESFINDER(ch_dir_pred) | RESFINDER_PARSE
  EMAPPER(ch_dir_pred) | EMAPPER_PARSE
  CAZY(ch_dir_pred)
  CAZY_PARSE1(CAZY.out)
  CAZY_PARSE2(CAZY_PARSE1.out[0])
  SWISSPROT(ch_dir_pred)
  PLASMID_BLAST(ch_dir_pred) | PLASMID_FILTER
}



include { KEGG;KEGG_ANNO_STATS } from '../modules/kegg.nf'


workflow VIR_ANNOTATION {
  take:
  ch_dir_pred
  ch_depth_pred

  main:
  SWISSPROT(ch_dir_pred)
  KEGG(ch_dir_pred)
  KEGG_ANNO_STATS(KEGG.out)
}