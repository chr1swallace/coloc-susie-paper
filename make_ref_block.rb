#!/home/cew54/localc/bin/ruby
require 'zlib'
require 'multiple_files_gzip_reader'
require 'tempfile'
require_relative 'dirs.rb'

block=ARGV[0]
outfile=ARGV[1] # plink root

## get and check list of snps
vcffile="#{DIR}/reference/byblock/#{block}.vcf.gz"
puts "subsetting vcf file #{vcffile}"

samplefile="/home/cew54/share/Data/reference/1000GP_Phase3/sparse_basis/EUR.sample" # for now, can make an option later

# vcftemp= `mktemp /tmp/tempvcf.XXXXXX`.chomp
command = "zcat #{vcffile} | " +
          "sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | " +
          "#{ENV['HOME']}/localc/bin/vcftools " +
          " --gzvcf - --IMPUTE  --out #{outfile} " +
          " --remove-indels --remove-filtered-all --keep #{samplefile} " +
          " --maf 0.01 --max-alleles 2 "
system(command)

# ## next: make plink
# command = "#{ENV['HOME']}/localc/bin/plink --vcf #{vcftemp} --out #{outfile}"
# system(command)

# ## cleanup
# File.unlink(vcftemp)

exit 0;
