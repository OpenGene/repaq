#ifndef RFQCHUNK_H
#define RFQCHUNK_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include <iostream>
#include <fstream>
#include "rfqheader.h"

using namespace std;

/*
* the query name of FASTQ is usually like
* @A00250:26:H3YTWDSXX:1:1101:16776:17989 1:N:0:ACTGTTCC
* can be divided into 
* <QNAME> = <NAME1>:<LANE>:<TILE>:<X>:<Y><NAME2>
* <NAME1> = @A00250:26:H3YTWDSXX:
* <NAME2> =  1:N:0:ACTGTTCC
*/

// flags
// if set, all reads share same length, so length array only has one element
#define BIT_READ_LEN_SAME (1<<0)
// if set, all name1 share same length, so name1 length array only has one element
#define BIT_NAME1_LEN_SAME (1<<1)
// if set, all name2 share same length, so name2 length array only has one element
#define BIT_NAME2_LEN_SAME (1<<2)
// if set, all strand share same length, so strand length array only has one element
#define BIT_STRAND_LEN_SAME (1<<3)
// if set, all lane number are same, so lane array only has one element
#define BIT_LANE_SAME (1<<4)
// if set, all tile number are same, so tile array only has one element
#define BIT_TILE_SAME (1<<5)
// if set, all name1 are same, so name1 array only has one element
#define BIT_NAME1_SAME (1<<6)
// if set, all name2 are same, so name2 array only has one element
#define BIT_NAME2_SAME (1<<7)
// if set, all strand (3rd line in FASTQ) are same, so strand array only has one element
#define BIT_STRAND_SAME (1<<8)
// in this BIT_PE_INTERLEAVED mode
// the reads are paired-end, interleaved, read2 sequence is reversed and complemented
// tile, X and Y information are only stored once in odd read1
// name2 are only stored once in read1, and name2 can be got by transfering read1 to read2 by replace 1 with 2 
#define BIT_PE_INTERLEAVED (1<<9)

class RfqChunk{
public:
    RfqChunk(RfqHeader* header);
    ~RfqChunk();
    void read(ifstream& ifs);
    void write(ofstream& ofs);
    void calcTotalBufSize();

private:
    void readReadLenBuf(ifstream& ifs);
    void readName1LenBuf(ifstream& ifs);
    void readName2LenBuf(ifstream& ifs);
    void readStrandLenBuf(ifstream& ifs);
    void readLaneBuf(ifstream& ifs);
    void readTileBuf(ifstream& ifs);

public:
    // the entire buffer size of this chunk
    uint32 mSize;
    // how many reads in this chunk
    uint32 mReads;
    // read length bit...
    uint16 mFlags;
    // size of encoded sequence buffer
    uint32 mSeqBufSize;
    // size of encoded quality buffer
    uint32 mQualBufSize;
    // buffers
    uint8* mReadLenBuf;
    uint8* mName1LenBuf;
    uint8* mName2LenBuf;
    uint8* mStrandLenBuf;
    uint8* mLaneBuf;
    uint16* mTileBuf;
    uint16* mXBuf;
    uint16* mYBuf;
    char* mName1Buf;
    char* mName2Buf;
    char* mSeqBuf;
    uint8* mQualBuf;
    char* mStrandBuf;

    // buffers
    uint32 mReadLenBufSize;
    uint32 mName1LenBufSize;
    uint32 mName2LenBufSize;
    uint32 mStrandLenBufSize;
    uint32 mLaneBufSize;
    uint32 mTileBufSize;
    uint32 mName1BufSize;
    uint32 mName2BufSize;
    uint32 mStrandBufSize;

    RfqHeader* mHeader;
};

#endif