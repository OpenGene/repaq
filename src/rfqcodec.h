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
    uint32 encodeSeqQual(char* seq, uint8* qual, char* seqEncoded, char* qualEncoded, uint32 seqLen, uint32 quaLen);
    uint32 encodeQualRunLenCoding(char* seq, uint8* qual, char* seqEncoded, char* qualEncoded, uint32 seqLen, uint32 quaLen);
    uint32 encodeQualByCol(char* seq, uint8* qual, char* seqEncoded, char* qualEncoded, uint32 seqLen, uint32 quaLen);
    uint32 encodeSingleQualByCol(uint8* qual, uint8 q, uint8* encoded, uint32 seqLen, uint32 quaLen);
    void decodeSeqQual(RfqChunk* chunk, string& seq, string& qual, uint32 len, uint32* readLenBuf);
    void decodeQualByRunLenCoding(RfqChunk* chunk, string& seq, string& qual, uint32 len);
    void decodeQualByCol(RfqChunk* chunk, string& seq, string& qual, uint32 len);
    void decodeSingleQualByCol(uint8* buf, uint32 bufLen, uint8 q, string& seq, string& qual);
    int overlap(string& r1, string& r2);

private:
    RfqHeader* mHeader;
};

#endif