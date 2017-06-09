# Bin Coverage

Make chart with coverage of contigs in a bin.

## Install

### From source

To compile from source, you need the [Chicken Scheme compiler](http://code.call-cc.org/).

Clone the repository, `cd` into the directory and type `make` into the terminal.  Optionally, you can type `make install` to install.

### Binaries

Check the [GitHub release tab](https://github.com/mooreryan/bin_coverage/releases) to get the latest pre-compiled binaries.

### Run the tests

Run `make test` to see if everything went right.

## Usage

### Requirements

- A version of Samtools that has the `-aa` option.
- `Rscript` program (to run the R scripts)

### Input files

- SAM or BAM file with recruitment to contigs
- Text file with one contig name per line.  These contigs will be treated as a single "bin" and included in the same coverage plot.

### Synopsis

```
USAGE: bin-cov /path/to/samtools recruitment.bam contig-names.txt
```