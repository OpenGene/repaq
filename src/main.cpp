#include <stdio.h>
#include <time.h>
#include "cmdline.h"
#include "common.h"
#include "options.h"
#include "unittest.h"
#include <sstream>
#include <iostream>
#include <unistd.h>
#include "repaq.h"
#include "util.h"

string command;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cerr << "repaq: repack FASTQ to a smaller binary file (.rfq)" << endl << "version " << VERSION_NUM << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in1", 'i', "input file name", false, "");
    cmd.add<string>("out1", 'o', "output file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name when encoding paired-end FASTQ files", false, "");
    cmd.add<string>("out2", 'O', "read2 output file name when decoding to paired-end FASTQ files", false, "");
    cmd.add("compress", 'c', "compress input to output");
    cmd.add("decompress", 'd', "decompress input to output");
    cmd.add<int>("chunk", 'k' , "the chunk size (kilo bases) for encoding, default 1000=1000kb.", false, 1000);
    cmd.add("stdin", 0, "input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.");
    cmd.add("stdout", 0, "write to STDOUT. When decompressing PE data, this option will result in interleaved FASTQ output for paired-end input. Disabled by defaut.");
    cmd.add("interleaved_in", 0, "indicate that <in1> is an interleaved paired-end FASTQ which contains both read1 and read2. Disabled by defaut.");
    cmd.add("compare", 'p', "compare the files read by read to check the compression consistency. <rfq_to_compare> should be specified in this mode.");
    cmd.add<string>("rfq_to_compare", 'r', "the RFQ file to be compared with the input. This option is only used in compare mode.", false, "");
    cmd.add<string>("json_compare_result", 'j', "the file to store the comparison result. This is optional since the result is also printed on STDOUT.", false, "");

    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }
    

    time_t t1 = time(NULL);

    Options opt;

    opt.in1 = cmd.get<string>("in1");
    opt.out1 = cmd.get<string>("out1");
    opt.in2 = cmd.get<string>("in2");
    opt.out2 = cmd.get<string>("out2");
    opt.rfqCompare = cmd.get<string>("rfq_to_compare");
    opt.jsonFileForCompare = cmd.get<string>("json_compare_result");
    opt.chunkSize = max(100, cmd.get<int>("chunk")) * 1000; // use at least 100Kb chunk
    opt.inputFromSTDIN = cmd.exist("stdin");
    opt.outputToSTDOUT = cmd.exist("stdout");
    opt.interleavedInput = cmd.exist("interleaved_in");

    int modeNum = 0;
    if(cmd.exist("compress"))
        modeNum++;
    if(cmd.exist("decompress"))
        modeNum++;
    if(cmd.exist("compare"))
        modeNum++;
    if(modeNum > 1)
        error_exit("repaq can run in compress/decompress/compare mode, you can only choose any one mode.");
    
    if(cmd.exist("decompress"))  {
        opt.mode = REPAQ_DECOMPRESS;
    }
    else if(cmd.exist("compare"))  {
        opt.mode = REPAQ_COMPARE;
    } else {
        // compress is the default mode
        opt.mode = REPAQ_COMPRESS;
    }

    if(opt.mode == REPAQ_COMPRESS && opt.outputToSTDOUT && !opt.out1.empty()) {
        cerr << "Output to STDOUT, ignore --out1 = " << opt.out1 << endl;
        opt.out1 = "";
    }

    if(opt.mode == REPAQ_DECOMPRESS && opt.inputFromSTDIN && !opt.in1.empty()) {
        cerr << "Input from STDIN, ignore --in1 = " << opt.in1 << endl;
        opt.in1 = "";
    }

    if(opt.mode == REPAQ_COMPARE && opt.inputFromSTDIN && !opt.rfqCompare.empty()) {
        cerr << "Input from STDIN, ignore --rfq_to_compare = " << opt.rfqCompare << endl;
        opt.rfqCompare = "";
    }

    opt.validate();

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    if(ends_with(opt.out1, ".xz") || ends_with(opt.in1, ".xz") || ends_with(opt.rfqCompare, ".xz")) {
        if(opt.inputFromSTDIN)
            error_exit("STDIN cannot be read in xz mode");
        if(opt.outputToSTDOUT)
            error_exit("STDOUT cannot be written in xz mode");
    }

    // deal with xz
    if(opt.mode == REPAQ_COMPRESS && ends_with(opt.out1, ".xz")) {
        replaceAll(command, opt.out1, "");
        replaceAll(command, "-o ", "");
        replaceAll(command, "--out1=", "");
        command = command + " --stdout | xz -z -c -9 > " + opt.out1;
        int ret = system(command.c_str());
        if(ret != 0) {
            error_exit("failed to call xz, please confirm that xz is installed in your system");
        }
    } else if(opt.mode == REPAQ_DECOMPRESS && ends_with(opt.in1, ".xz")) {
        replaceAll(command, opt.in1, "");
        replaceAll(command, "-i", "");
        replaceAll(command, "--in1=", "");
        command = "xz -d -c " + opt.in1 + "| " + command + " --stdin " ;
        int ret = system(command.c_str());
        if(ret != 0) {
            error_exit("failed to call xz, please confirm that xz is installed in your system");
        }
    }  else if(opt.mode == REPAQ_COMPARE && ends_with(opt.rfqCompare, ".xz")) {
        replaceAll(command, opt.rfqCompare, "");
        replaceAll(command, "-r", "");
        replaceAll(command, "--rfq_to_compare=", "");
        command = "xz -d -c " + opt.rfqCompare + "| " + command + " --stdin " ;
        int ret = system(command.c_str());
        if(ret != 0) {
            error_exit("failed to call xz, please confirm that xz is installed in your system");
        }
    } else {

        Repaq repaq(&opt);
        repaq.run();

        return 0;
    }
}