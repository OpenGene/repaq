# repaq
A tool to repack Illumina format FASTQ to a smaller binary file (.rfq), which can be further compressed by xz or pxz (.rfq.xz).   

For NovaSeq data, the .rfq file can be much smaller than .fq.gz, and the compressing time is usually less than 1/5 of gzip compression. 

The biggest advantage is that the .rfq file can be further compressed with xz, which is based on LZMA algorithm. The .rfq.xz file can be as small as 5% of the original FASTQ file, or smaller than 30% of the .fq.gz file. Note that usually the gz files are not compressible by xz.

This tool also supports non-Illumina format FASTQ (i.e. the BGI-SEQ format), but the compression ratio is not as good Illumina format FASTQ.

***WARNING: be careful about using repaq for production before v1.0 is released, since its spec v1.0 has not been frozen.***

# take a look of the compression ratio
Here we demonstrate the compression ratio of two paired-end NovaSeq data. You can download these files and test locally.
* `nova.R1.fq`: 1704 MB, the original read1 file, http://opengene.org/repaq/testdata/nova.R1.fq
* `nova.R2.fq`: 1704 MB, the original read2 file, http://opengene.org/repaq/testdata/nova.R2.fq
* `nova.R1.fq.gz`: 308 MB (CR 18.08%), the gzipped read1, http://opengene.org/repaq/testdata/nova.R1.fq.gz
* `nova.R2.fq.gz`: 325 MB (CR 19.07%), the gzipped read2, http://opengene.org/repaq/testdata/nova.R2.fq.gz
* `nova.rfq`: 341 MB (CR 10.01%), the repacked file of read1+read2, http://opengene.org/repaq/testdata/nova.rfq
* `nova.rfq.xz`: 134 MB (CR 3.93%), the xz compressed `nova.rfq`, http://opengene.org/repaq/testdata/nova.rfq.xz

See? The size of final `nova.rfq.xz` is only 3.39% of the original FASTQ files! You can decompress it and check the md5 to see whether they are identical! 

Typically with one single CPU core, it takes less than 1 minute to convert `nova.R1.fq + nova.R2.fq` to `nova.rfq`, and takes less than 5 minutes to compress the `nova.rfq` to `nova.rfq.xz` by xz.

# get repaq
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
# compress
repaq -c -i in.fq -o out.rfq

# decompress
repaq -c -i in.rfq -o out.fq
```

For paired-end mode:
```shell
# compress
repaq -d -i in.R1.fq -I in.R2.fq -o out.rfq

# decompress
repaq -d -i in.rfq -o out.R1.fq -O out.R2.fq
```

Tips:
* `-i` and `-I` always denote the first and second input files, while `-o` and `-O` always denote the first and second output files.
* the FASTQ input/output files can be gzipped if their names are ended with `.gz`.
* for paired-end data. the .rfq file created in paired-end mode is usually much smaller than the sum of the .rfq files created in single-end mode for R1 and R2 respectively. To obtain high compression rate, please always use PE mode for PE data.

# compress .rfq to .rfq.xz with xz
To get highest compression ratio (need at least 16G RAM):
```
xz --lzma2="dict=1000000000" in.rfq
```

To get normal ratio (need at least 1G RAM):
```
xz -9 in.rfq
```

The latest version of xz supports multithreading, so you can specify the thread number with `-T` option:
```
xz -T4 -9 in.rfq
```

You can also use pxz for parallel xz compressing:
```
pxz -9 in.rfq
```

Tips:
* lower compression ratio than `-9` is not recommended, since it will not be faster. The difference is the RAM requirement.
```
```

# FASTQ Format compatibility  
repaq was initially designed for compressing Illumina data, but it also works with data from other platforms, like BGI-Seq. To work with repaq, the FASTQ format should meet following condidtions:
* only has bases A/T/C/G/N.
* each FASTQ record has, and only has four lines (name, sequence, strand, quality).
* the name and strand line cannot be longer than 255 bytes.
* the number of different quality characters cannot be more than 127.

`repaq` works best for Illumina data directly output by `bcl2fastq`.
