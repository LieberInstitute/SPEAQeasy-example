#  Install 'getopt' if needed, without requiring manual work from the user
if (!requireNamespace('getopt')) {
    install.packages('getopt')
}

library('getopt')

spec = matrix(c(
    'dir', 'd', 1, 'character', 'directory containing FASTQ reads'
), byrow=TRUE, ncol=5)
opt = getopt(spec)

old_man = read.table('sample_selection/samples.manifest')

r1 = paste0(opt$dir, '/', basename(old_man[,1]))
r2 = paste0(opt$dir, '/', basename(old_man[,3]))

new_man = paste(r1, 0, r2, 0, old_man[,5], sep='\t')

writeLines(new_man, paste0(opt$dir, '/samples.manifest'))
