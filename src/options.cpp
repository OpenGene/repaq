#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string.h>

Options::Options(){
    in1 = "";
    out1 = "";
    in2 = "";
    out2 = "";
    chunkSize = 1000;
    compressMode = true;
}

bool Options::isFastqFile(string filename) {
    if(ends_with(filename, ".fq") || ends_with(filename, ".fastq") || ends_with(filename, ".fq.gz")|| ends_with(filename, ".fastq.gz"))
        return true;

    return false;
}

bool Options::isRfqFile(string filename) {
    if(ends_with(filename, ".rfq") || ends_with(filename, ".rfq.xz"))
        return true;

    return false;
}

bool Options::validate() {
    if(in1.empty()) {
        error_exit("read1 input should be specified by --in1");
    } else {
        check_file_valid(in1);
    }

    if(!in2.empty()) 
        check_file_valid(in2);

    if(out1.empty()) {
        error_exit("read1 output should be specified by --out1");
    }

    if(compressMode == true) {
        if(!out2.empty())
            error_exit("In compress mode, only one RFQ output file is allowed, but you specified <out2>");
        if(!isRfqFile(out1))
            error_exit("In compress mode, the output should be a RFQ file. Expect a .rfq file, but got " + out1);
        if(isFastqFile(out1))
            error_exit("In compress mode, the output should not be a FASTQ file. Expect a .rfq or .rfq.xz file, but got " + out1);
        if(isRfqFile(in1))
            error_exit("In compress mode, the input should not be a RFQ file. Expect a .fq or .fq.gz file, but got " + in1);
        if(!in2.empty() && isRfqFile(in2))
            error_exit("In compress mode, the read2 input should not be a RFQ file. Expect a .fq or .fq.gz file, but got " + in2);
    }

    if(compressMode == false) {
        if(!in2.empty())
            error_exit("In decompress mode, only one RFQ input file is allowed, but you specified <in2>");
        if(!isRfqFile(in1))
            error_exit("In decompress mode, the input should be a RFQ file. Expect a .rfq file, but got " + out1);
        if(isFastqFile(in1))
            error_exit("In decompress mode, the input should not be a FASTQ file. Expect a .rfq or .rfq.xz file, but got " + in1);
        if(isRfqFile(out1))
            error_exit("In decompress mode, the output should not be a RFQ file. Expect a .fq or .fq.gz file, but got " + out1);
        if(!out2.empty() && isRfqFile(out2))
            error_exit("In decompress mode, the read2 output should not be a RFQ file. Expect a .fq or .fq.gz file, but got " + out2);
    }

    if(!out2.empty()) {
        if(compressMode)
            error_exit("In compress mode, only one rfq output file is allowed");
    }

    if(ends_with(in1, ""))

    if(chunkSize < 10000) {
        error_exit("chunk size cannot be less than 10 kb");
    }
    if(chunkSize >500000000) {
        error_exit("chunk size cannot be greater than 500,000 kb");
    }

    return true;
}