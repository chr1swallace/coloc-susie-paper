require 'rake'
require '/home/cew54/slurmer/Qsub.rb'
require 'rake/clean'
require 'fileutils'
require 'date'
require_relative 'dirs.rb' # defines DIR, ROOTDIR, REFDIR

## CLEAN

# "clean temporary files"
# task :clean do
# end

CLEAN.include('machine.*')
CLEAN.include('slurm-*.out')
CLEAN.include('slurm-*.sh*')
CLEAN.include('rake-2021-*.sh*')

QFILE="rake-#{DateTime.now}.sh"
qfile_=open(QFILE,'w')
qrun=false
qstr=" -r " # q options

at_exit do
    qfile_.close
    if qrun
        n=`cat #{QFILE} | wc -l`.chomp
        puts "\nTo run #{n} jobs:\nqlines_asarray.rb #{qstr} #{QFILE}"
    else
        File.unlink(QFILE)
    end
end

args = {#:job=>"impute",
    :time=>'4:00:00',
    :tasks=>1,
    :cpus=>1,
    :autorun=>true,
    :interactiverun=> ! ENV['INTERACT'].nil?,
    :excl=>" "}

require 'bundler/setup'
$:.unshift File.expand_path('../lib', __FILE__)

# rake spec
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) { |t| t.verbose = false   }

## read in blocks
require_relative '/home/cew54/E/blocks.rb'
BLOCKS=Blocks.new()

# rake console
task :console do
    require 'pry'
    ARGV.clear
    Pry.start
end

namespace "coloc" do

    task :run do
        qrun=true
        qstr += " -j simcoloc -g 50 "
    end
    # NSIMDATA=100 # how many data files
    NSIMCOLOC=2000 # how many coloc runs
    NPERSIM=10 # how many simulations per output file
    # desc "sim coloc data"
    # task :data do
    #     qrun=true
    #     qstr += " -j simdata -g 2 "
    #     # how many data files?
    #     simfiles=Dir.glob("*", base: SIMDATA)
    #     i=0
    #     while i < (NSIMDATA - simfiles.length) do
    #         simblock=BLOCKS.label.sample
    #         # locs=BLOCKS.label2be(simblock)
    #         qfile_.puts("./sim_data.R --args block=#{simblock}")
    #         i=i+1
    #     end
    # end

    # A is causal for trait 1, B is causal for trait 2.
    # +1 adds the other snp, -1 switches to the other snp
    # SCENARIOS indicates whhether A is shared, B is shared
    SCENARIOS=[
        [1,-1], # A,A shared
        [0,0], # A,B, no sharing
        [1,0], # A,AB one shared
        [1,1], # AB,AB both shared
    ]

    SCENARIOS.each do |s|
        # how many data files?
        [1,2,3].each do |ncopies|
            d="#{SIMCOLOC}/sim_ncopies_#{ncopies}_A_#{s[0]}_B_#{s[1]}"

            desc "clean earlier output files"
            task :clean do
                FileUtils.rm_r d if File.exist?(d);
            end

            desc "run coloc on sim data"
            task :run do
                FileUtils.mkdir_p(d) # does nothing if d exists
                simfiles=Dir.glob("*", base: d)
                i=0
                while i < ((NSIMCOLOC/NPERSIM - simfiles.length)).ceil do
                    qfile_.puts("Rscript ./sim_coloc.R --args Ashared=#{s[0]} Bshared=#{s[1]} nsim=#{NPERSIM} ncopies=#{ncopies}")
                    i=i+1
                end
            end
        end
    end

    desc "summarise sims to date"
    task :summary do
        qrun=true
        qstr += " -j simsumm "
        qfile_.puts("./summary_coloc.R")
    end
end

namespace :approx do
    NAPPROX=200

    desc "run approximation simulations"
    task :run do
        qrun=true
        qstr += " -j simsusie "
    end

    [1,2,3].each do |n|
        d="#{SIMCOLOC}/approx_#{n}"

        desc "clean earlier output files"
        task :clean do
            FileUtils.rm_r d if File.exist?(d);
        end

        desc "run approximation simulations"
        task :run do
            FileUtils.mkdir_p(d) # does nothing if d exists
            simfiles=Dir.glob("*", base: d)
            i=0
            while i < (NAPPROX - simfiles.length) do
            # while i < (NAPPROX) do
                qfile_.puts("Rscript ./sim_susie.R --args Bshared=1 nsim=10 ncopies=#{n}")
                qfile_.puts("Rscript ./sim_susie.R --args Bshared=0 nsim=10 ncopies=#{n}")
                i=i+1
            end
        end
    end

# BaseQTL is a Bayesian method to map molecular QTL affecting allele-specific expression even when no genotypes are available. It is well suited to discover eQTLs hidden in a wealth of RNA-seq data to unravel molecular mechanisms underpinning disease.

    desc "summarise susie sims to date"
    task :summary do
        qrun=true
        qstr += " -j sussum "
        qfile_.puts("./summary_susie.R")
    end
end
