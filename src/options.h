#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

#define REPAQ_COMPRESS 0
#define REPAQ_DECOMPRESS 1
#define REPAQ_COMPARE 2

class Options{
public:
    Options();
    bool validate();

    static bool isFastqFile(string filename);
    static bool isRfqFile(string filename);

public:
    // IO
    string in1;
    string out1;
    string in2;
    string out2;
    string rfqCompare;
    string jsonFileForCompare;
    bool inputFromSTDIN;
    bool outputToSTDOUT;
    // the input R1 file is interleaved
    bool interleavedInput;

    // chunk
    int chunkSize;
    int mode;

};

#endif