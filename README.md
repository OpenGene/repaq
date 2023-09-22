[![install with conda](
https://anaconda.org/bioconda/repaq/badges/version.svg)](https://anaconda.org/bioconda/repaq)
# repaq
A tool to compress FASTQ files with ultra-high compression ratio and high speed. `repaq` supports compressing the FASTQ to `.rfq` or `.rfq.xz` formats. Compressing to `.rfq` is ultra fast, while compressing to `.rfq.xz` provides very high compression ratio. 

For NovaSeq data, as an example:  
* the `.rfq` file can be much smaller than `.fq.gz`, and the compressing time is usually less than 1/5 of gzip compression.
* The `.rfq.xz` file can be as small as 5% of the original FASTQ file, or smaller than 30% of the `.fq.gz` file.

For paired-end FASTQ files, `repaq` compresses them into one single file to provide higher compression ratio.

This tool also supports non-Illumina format FASTQ (i.e. the BGI-SEQ format), but the compression ratio is not as good Illumina format FASTQ.

*Citation: Chen S, Chen Y, Wang Z, Qin W, Zhang J, Nand H, Zhang J, Li J, Zhang X, Liang X and Xu M (2023) Efficient sequencing data compression and FPGA acceleration based on a two-step framework. Front. Genet. 14:1260531. doi: 10.3389/fgene.2023.1260531*   

# take a look at the compression ratio
Here we demonstrate the compression ratio of two paired-end NovaSeq data. You can download these files and test locally.
* `nova.R1.fq`: 1704 MB, the original read1 file, http://opengene.org/repaq/testdata/nova.R1.fq
* `nova.R2.fq`: 1704 MB, the original read2 file, http://opengene.org/repaq/testdata/nova.R2.fq
* `nova.R1.fq.gz`: 308 MB (CR 18.08%), the gzipped read1, http://opengene.org/repaq/testdata/nova.R1.fq.gz
* `nova.R2.fq.gz`: 325 MB (CR 19.07%), the gzipped read2, http://opengene.org/repaq/testdata/nova.R2.fq.gz
* `nova.rfq`: 333 MB (CR 9.77%), the repacked file of read1+read2, http://opengene.org/repaq/testdata/nova.rfq
* `nova.rfq.xz`: 134 MB (CR 3.93%), the xz compressed `nova.rfq`, http://opengene.org/repaq/testdata/nova.rfq.xz

See? The size of final `nova.rfq.xz` is only 3.39% of the original FASTQ files! You can decompress it and check the md5 to see whether they are identical! 

Typically with one single CPU core, it takes less than 1 minute to convert `nova.R1.fq + nova.R2.fq` to `nova.rfq`, and takes less than 5 minutes to compress the `nova.rfq` to `nova.rfq.xz` by xz.

# get repaq
## install with Bioconda
[![install with conda](
https://anaconda.org/bioconda/repaq/badges/version.svg)](https://anaconda.org/bioconda/repaq)
```shell
conda install -c bioconda repaq
```
## download binary 
This binary is only for Linux systems: http://opengene.org/repaq/repaq
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/repaq/repaq
chmod a+x ./repaq
```
## or compile from source
```shell
# get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/repaq.git

# build
cd repaq
make

# Install
sudo make install
```

# usage
For single-end mode:
```shell
# compress to .rfq.xz
repaq -c -i in.fq -o out.rfq.xz

# decompress from .rfq.xz
repaq -d -i in.rfq.xz -o out.fq
```

For paired-end mode:
```shell
# compress to .rfq.xz
repaq -c -i in.R1.fq -I in.R2.fq -o out.rfq.xz

# decompress from .rfq.xz
repaq -d -i in.rfq.xz -o out.R1.fq -O out.R2.fq
```

Tips:
* `-i` and `-I` always denote the first and second input files, while `-o` and `-O` always denote the first and second output files.
* the FASTQ input/output files can be gzipped if their names are ended with `.gz`.
* for paired-end data. the .rfq file created in paired-end mode is usually much smaller than the sum of the .rfq files created in single-end mode for R1 and R2 respectively. To obtain high compression rate, please always use PE mode for PE data.
* if you want higher speed and are not concern with compression ratio, replace `xxx.rfq.xz` with `xxx.rfq`, then repaq will compress or decompress `.rfq` format.

# system requirements
* Memory: 16G RAM
* CPU: 4 cores

# verify the compressed file
repaq offers a `compare` mode to check the consistency of the original FASTQ file(s) and the compressed .rfq or .rfq.xz file. 
* set `--compare` to enable the `compare` mode
* specify the .rfq or .rfq.xz file by `-r` option
* specify the FASTQ files by `-i` and `-I` options.

Examples:
```shell
# for single-end data
repaq --compare -i original.R1.fq  -r compressed.rfq.xz

# for paired-end data
repaq --compare -i original.R1.fq.gz -I original.R2.fq.gz  -r compressed.rfq.xz
```
Without any expection, you will get an output of a JSON like:
```json
{
	"result":"passed",
	"msg":"",
	"fastq_reads":50000,
	"rfq_reads":50000,
	"fastq_bases":7419082,
	"rfq_bases":7419082
}
```
The `result` will be "failed" if the compressed file is not consistent with the original FASTQ files.

# STDIN and STDOUT
repaq can read the input from STDIN, and write the output to STDOUT.
* specify `--stdin` if you want to read the STDIN for compression or decompression.
* specify `--stdout` if you want to output to the STDOUT for compression or decompression
* in decompression mode, if `--stdout` is specified, the output will be interleaved PE stream.
* if the STDIN is an interleaved paired-end stream, specify `--interleaved_in` to indicate that.
* be noted that STDIN cannot be read when the input is a .xz file, and STDOUT cannot be written when the output is a .xz file

Here gives you an example of compressing the interleaved PE output from fastp by directly using pipes:
```shell
fastp -i R1.fq -I R2.fq --stdout | repaq -c --interleaved_in --stdin -o out.rfq.xz
```

# FASTQ Format compatibility  
repaq was initially designed for compressing Illumina data, but it also works with data from other platforms, like BGI-Seq. To work with repaq, the FASTQ format should meet following condidtions:
* only has bases A/T/C/G/N.
* each FASTQ record has, and only has four lines (name, sequence, strand, quality).
* the name and strand line cannot be longer than 255 bytes.
* the number of different quality characters cannot be more than 127.

`repaq` works best for Illumina data directly output by `bcl2fastq`.

# all options
```shell
options:
  -i, --in1                    input file name (string [=])
  -o, --out1                   output file name (string [=])
  -I, --in2                    read2 input file name when encoding paired-end FASTQ files (string [=])
  -O, --out2                   read2 output file name when decoding to paired-end FASTQ files (string [=])
  -c, --compress               compress input to output
  -d, --decompress             decompress input to output
  -k, --chunk                  the chunk size (kilo bases) for encoding, default 1000=1000kb. (int [=1000])
      --stdin                  input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
      --stdout                 write to STDOUT. When decompressing PE data, this option will result in interleaved FASTQ output for paired-end input. Disabled by defaut.
      --interleaved_in         indicate that <in1> is an interleaved paired-end FASTQ which contains both read1 and read2. Disabled by defaut.
  
# following options are used to check the consistency of the compressed data
  -p, --compare                compare the files read by read to check the compression consistency. <rfq_to_compare> should be specified in this mode.
  -r, --rfq_to_compare         the RFQ file to be compared with the input. This option is only used in compare mode. (string [=])
  -j, --json_compare_result    the file to store the comparison result. This is optional since the result is also printed on STDOUT. (string [=])

# options for .xz output
  -t, --thread                 thread number for xz compression. Higher thread num means higher speed and lower compression ratio (1~16), default 1. (int [=1])
  -z, --compression            compression level. Higher level means higher compression ratio, and more RAM usage (1~9), default 4. (int [=4])

  -?, --help                   print this message
```

# external dependency
`repaq` makes a system call in order to run the xz compression tool available on GNU/Linux systems. If xz isn't installed, `repaq` will fail with the message: 

```
failed to call xz, please confirm that xz is installed in your system
```
