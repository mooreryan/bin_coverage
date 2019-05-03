#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "abort_if"
require "fileutils"
require "aai"
require "optimist"
require "set"

module Utils
  extend Aai::CoreExtensions::Time
  extend Aai::CoreExtensions::Process
end

include AbortIf

VERSION = "v0.3.1"

opts = Optimist.options do
  version VERSION

  banner <<-EOS

  Version: #{VERSION}

  Names file is two tab delimited columns.  First column is bin name,
  second column is contig name.

  --fast will probably break with large set of input contigs.  If you
    get some weird errors or things aren't working with --fast
    enabled, and you have a lot of contigs, try without --fast.

  Options:
  EOS

  opt(:samtools, "Path to samtools", default: `which samtools`.chomp)
  opt(:bam, "Bam file", type: :string)
  opt(:names, "Names file", type: :string)
  opt(:outdir, "Output directory", default: ".")
  opt(:outbase, "Basename for output", default: "snazzy")
  opt(:fast, "Use this option to filter contigs with samtools view")
end

abort_unless opts[:bam] && File.exist?(opts[:bam]),
             "#{opts[:bam]} does not exist. Try --help"

abort_unless opts[:names] && File.exist?(opts[:names]),
             "#{opts[:names]} does not exist. Try --help"

FileUtils.mkdir_p opts[:outdir]

bin2contigs = {}
contig_cov = {}
contig2start_posn = {}
all_contigs = Set.new

Utils.time_it "Reading names file", AbortIf::logger do
  File.open(opts[:names], "rt").each_line do |line|
    bin, contig = line.chomp.split "\t"

    all_contigs << contig

    unless bin2contigs.has_key? bin
      bin2contigs[bin] = Set.new
    end

    bin2contigs[bin] << contig
  end
end

depth_fname = File.join opts[:outdir], opts[:outbase] + ".depth.txt"

if opts[:fast]
  # If bam index can't be found, make one
  unless File.exist? "#{opts[:bam]}.bai"
    cmd = "#{opts[:samtools]} index #{opts[:bam]}"
    Utils.run_and_time_it! "Generating bam index", cmd
  end

  contig_str = all_contigs.to_a.join " "

  cmd = "#{opts[:samtools]} view -h #{opts[:bam]} #{contig_str} " \
        "| #{opts[:samtools]} depth -aa - > #{depth_fname}"
else
  cmd = "#{opts[:samtools]} depth -aa #{opts[:bam]} > #{depth_fname}"
end
Utils.run_and_time_it! "Calculating coverage", cmd


Utils.time_it "Reading depth file", AbortIf::logger do
  File.open(depth_fname, "rt").each_line do |line|
    contig, posn, cov = line.chomp.split "\t"

    if all_contigs.include? contig
      unless contig_cov.has_key? contig
        contig_cov[contig] = []
      end

      # The assumption is that samtools depth -aa gives a coverage for
      # every position.
      contig_cov[contig] << cov
    end
  end
end

Utils.time_it "Writing graph lines", AbortIf::logger do
  bin2contigs.keys.each_with_index do |bin, idx|
    contigs = bin2contigs[bin]
    graph_lines_fname = File.join opts[:outdir], opts[:outbase] + ".#{bin}_lines"
    start_posns_fname = File.join opts[:outdir], opts[:outbase] + ".#{bin}_start_posns"
    png_fname = File.join opts[:outdir], opts[:outbase] + ".#{bin}_coverage_plot.png"
    rscript_fname = File.join opts[:outdir], opts[:outbase] + ".#{bin}_plotter.r"

    overall_posn = 1

    File.open(graph_lines_fname, "w") do |lines_f|
      File.open(start_posns_fname, "w") do |posns_f|
        lines_f.puts "position\tcoverage"
        posns_f.puts "contig\tstart"

        contigs.each do |contig|
          posns_f.puts "#{contig}\t#{overall_posn}"

          abort_unless contig_cov.has_key?(contig),
                       "Missing #{contig} from contig_cov hash table"

          contig_cov[contig].each do |cov|
            lines_f.puts "#{overall_posn}\t#{cov}"
            overall_posn += 1
          end
        end
      end
    end

    rscript = %Q{
dat <- read.table("#{graph_lines_fname}", sep="\\t", header=T);
start.posns <- read.table("#{start_posns_fname}", sep="\\t", header=T);
width.mult <- nrow(dat) / 50000
png("#{png_fname}", units = "in", res=120, width = max(1 * width.mult, 11), height=5);
plot(dat, type = "l", xlab="Position", ylab = "Coverage", main="#{bin}");
invisible(lapply(start.posns$start, FUN=function(pos){abline(v=pos, col="blue")}));
names <- start.posns$contig
for(idx in 1:length(names)){text(cex = 0.75, srt=90, col="red", pos=4, x=start.posns$start[idx], y = max(dat$coverage)/2, names[idx])};
invisible(dev.off());
}

    File.open(rscript_fname, "w") do |f|
      f.puts rscript
    end

    Utils.run_it! "Rscript #{rscript_fname}"
  end
end
