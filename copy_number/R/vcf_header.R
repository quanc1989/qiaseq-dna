## LINES TO PUT INTO VCF HEADER
vcf.header <- function (chromosomes) {
out <- '##fileformat=VCFv4.1
##INFO=<ID=CHROM,Number=1,Type=String,Description="Chromosome name">
##INFO=<ID=POS,Number=1,Type=Integer,Description="Position in chromosome">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=0,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=“Imprecise structural variation”>
##INFO=<ID=END,Number=1,Type=Integer,Description=“End position of the variant described in this record”>
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##ALT=<ID=CNV,Description="Copy number variable region">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=“Copy number genotype for imprecise events”>
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=“Copy number genotype quality for imprecise events”>
##FORMAT=<ID=CNL,Number=.,Type=Float,Description=“Copy number genotype likelihood for imprecise events”>'

#  for (chr in chromosomes) {
#    
#  }
  return(out)
}
