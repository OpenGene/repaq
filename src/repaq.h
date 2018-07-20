#ifndef REPAQ_H
#define REPAQ_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "rfqcodec.h"
#include "options.h"

using namespace std;

class Repaq{
public:
    Repaq(Options* opt);
    void run();
    void compress();
    void compressPE();
    void decompress();
    void decompressPE();

private:
    bool hasLineBreakAtEnd(string& filename);

private:
    Options* mOptions;
};

#endif