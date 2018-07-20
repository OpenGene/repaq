#include <stdio.h>
#include <time.h>
#include "cmdline.h"
#include "common.h"
#include "options.h"
#include "unittest.h"
#include <sstream>
#include "repaq.h"
#include "util.h"

string command;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cout << "repaq: repack FASTQ to a smaller binary file (.rfq)" << endl << "version " << VERSION_NUM << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in1", 'i', "input file name", true, "");
    cmd.add<string>("out1", 'o', "output file name", true, "");
    cmd.add<string>("in2", 'I', "read2 input file name when encoding paired-end FASTQ files", false, "");
    cmd.add<string>("out2", 'O', "read2 output file name when decoding to paired-end FASTQ files", false, "");
    cmd.add("compress", 'c', "compress input to output");
    cmd.add("decompress", 'd', "decompress input to output");
    cmd.add<int>("chunk", 0 , "the chunk size (kilo bases) for encoding, default 100=100kb.", false, 100);

    cmd.parse_check(argc, argv);
    
    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    Options opt;

    opt.in1 = cmd.get<string>("in1");
    opt.out1 = cmd.get<string>("out1");
    opt.in2 = cmd.get<string>("in2");
    opt.out2 = cmd.get<string>("out2");
    opt.chunkSize = cmd.get<int>("chunk") * 1000;

    if(cmd.exist("compress") && cmd.exist("decompress"))
        error_exit("You cannot specify both compress and decompress, please choose either one");
    else if(cmd.exist("decompress")) {
        opt.compressMode = false;
    }
    else {
        opt.compressMode = true;
    }

    opt.validate();

    Repaq repaq(&opt);
    repaq.run();

    return 0;
}