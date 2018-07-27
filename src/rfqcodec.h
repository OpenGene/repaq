#ifndef RFQCODEC_H
#define RFQCODEC_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include <iostream>
#include <fstream>
#include "rfqheader.h"
#include "rfqchunk.h"
#include "read.h"

using namespace std;

class RfqCodec{
public:
    RfqCodec();
    ~RfqCodec();
    void setHeader(RfqHeader* header);
    RfqHeader* makeHeader(vector<Read*>& reads);
    RfqHeader* makeHeader(vector<ReadPair*>& pairs);
    RfqChunk* encodeChunk(vector<Read*>& reads, bool isPE = false);
    RfqChunk* encodeChunk(vector<ReadPair*>& pairs);
    vector<Read*> decodeChunk(RfqChunk* chunk);

private:
    uint32 encodeSeqQual(char* seq, uint8* qual, char* seqEncoded, char* qualEncoded, uint32 totalLen);
    void decodeSeqQual(RfqChunk* chunk, string& seq, string& qual, uint32 len);
    int overlap(string& r1, string& r2);

private:
    RfqHeader* mHeader;
};

#endif