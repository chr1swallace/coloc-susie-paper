#!/home/cew54/localc/bin/ruby
require 'zlib'
require 'multiple_files_gzip_reader'
require 'tempfile'
require_relative 'dirs.rb'

outfile="ctsh_data/haps.vcf.gz"

## get and check list of snps
vcffile="#{REFDIR}/ALL.chr19.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
puts "subsetting vcf file #{vcffile}"

samplefile="/home/cew54/share/Data/reference/1000GP_Phase3/sparse_basis/EUR.sample" # for now, can make an option later

# vcftemp= `mktemp /tmp/tempvcf.XXXXXX`.chomp
command = "zcat #{vcffile} | " +
          "sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | " +
          "#{ENV['HOME']}/localc/bin/vcftools " +
          " --gzvcf - --IMPUTE  --out #{outfile} " +
          " --chr 19 --from-bp 78900000 --to-bp 79000000 " +
          " --remove-indels --remove-filtered-all --keep #{samplefile} " +
          " --maf 0.01 --max-alleles 2 "
system(command)

# ## next: make plink
# command = "#{ENV['HOME']}/localc/bin/plink --vcf #{vcftemp} --out #{outfile}"
# system(command)

# ## cleanup
# File.unlink(vcftemp)

exit 0;
