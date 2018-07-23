# repaq
A tool to repack FASTQ to a smaller binary file (.rfq), which can be further compressed by xz  (.rfq.xz).   

For NovaSeq data, the .rfq file can be much smaller than .fq.gz, and the compressing time is usually less than 1/5 of gzip compression. 

The biggest advantage is that the .rfq file can be further compressed with xz, which is based on LZMA algorithm. The .rfq.xz file can be as small as 5% of the original FASTQ file, or smaller than 30% of the .fq.gz file. Note that usually the gz files are not compressible by xz.

* WARNING The spec v1.0 has not been locked. So please don't use it for production until v1.0 is released.*

# Get repaq
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
repaq -i in.fq -o out.rfq

# decompress
repaq -i in.rfq -o out.fq
```

For paired-end mode:
```shell
# compress
repaq -i in.R1.fq -I in.R2.fq -o out.rfq

# decompress
repaq -i in.rfq -o out.R1.fq -O out.R2.fq
```

Tips:
* `-i` and `-I` always denote the first and second input files, while `-o` and `-O` always denote the first and second output files.
* the FASTQ input/output files can be gzipped if their names are ended with `.gz`.
* for paired-end data. the .rfq file created in paired-end mode is usually smaller than the sum of the .rfq files created in single-end mode for R1 and R2 respectively.

# compress .rfq to .rfq.xz with xz
To get highest compression ratio (need at least 16G RAM)
```
xz --lzma2="dict=1000000000" in.rfq
```

To get normal ratio (need at least 1G RAM)
```
xz -9 in.rfq
```
Tips:
* you can add `-T4` to use 4 threads for parallel compression if your system has enough RAM.
* lower compression ratio than `-9` is not recommended, since it will not be faster. The difference is the RAM requirement.
```
```

# FASTQ Format compatibility  
repaq was initially designed for compressing Illumina data, but it also works with data from other platforms, like BGI-Seq. To work with repaq, the FASTQ format should meet following condidtions:
* only has bases A/T/C/G/N.
* each FASTQ record has, and only has four lines (name, sequence, strand, quality).
* the name and strand line cannot be longer than 255 bytes.
* if one quality is with N base, it cannot be with other bases. For example, if `#` is with N base, it cannot be with `A/T/C/G`.
* the number of different quality characters cannot be more than 64.

`repaq` works best for Illumina data directly output by `bcl2fastq`.
