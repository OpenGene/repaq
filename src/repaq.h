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
    void compare();
    void comparePE();

private:
    bool hasLineBreakAtEnd(string& filename);
    void reportCompareResult(bool passed, string message, long fqReads, long fqBases, long rfqReads, long rfqBases);
    bool doubleCheckAndOutput(RfqChunk* chunk, RfqCodec& codec4check, RfqHeader* header4check, vector<Read*>& reads, ostream& out);
    bool doubleCheckAndOutput(RfqChunk* chunk, RfqCodec& codec4check, RfqHeader* header4check, vector<ReadPair*>& pairs, ostream& out);

private:
    Options* mOptions;
};

#endif