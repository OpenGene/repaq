#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

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

    // chunk
    int chunkSize;
    bool compressMode;

};

#endif