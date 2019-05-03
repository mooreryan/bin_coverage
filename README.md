# Bin Coverage

Make chart with coverage of contigs in bunches of bins.

## Install

### Dependencies

The following RubyGems are required.

- aai
- abort_if
- optimist

Install them like so

```bash
$ gem install aai abort_if optimist
```

Also, you need a version of Samtools that supports the `-aa` option, and the `Rscript` program to run R from the command line.

### Get the code

Using git

```bash
$ git clone https://github.com/mooreryan/bin_coverage.git
```

In that folder you will have a script called `bin_coverage.rb`.  That's all you need.  If you want, you can symlink it somewhere on your path.  Perhaps like this

```bash
$ ln -s $PWD/bin_coverage.rb /usr/local/bin/bin_coverage.rb
```

## Usage

To see the help message, type

```bash
$ bin_coverage.rb -h
```

### Input files

- BAM file with recruitment to contigs
- Two column, tab delimited text file.  The first column is bin name, the second column is contig name.  Contigs in the same bin will show up in the same coverage plot.  Don't include the `>` in the name of the sequence. Note, if the names don't seem to match, try to include only up to the first space in the names file.
  - The bam file can have more contigs than are listed in the names file.  Any extra contigs in the bam file that are not in the names file will be ignored.
