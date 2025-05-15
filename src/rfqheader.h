#ifndef RFQHEADER_H
#define RFQHEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include <iostream>
#include "read.h"

using namespace std;

/*
* the query name of FASTQ is usually like
* @A00250:26:H3YTWDSXX:1:1101:16776:17989 1:N:0:ACTGTTCC
* can be divided into 
* <QNAME> = <NAME1>:<LANE>:<TILE>:<X>:<Y><NAME2>
* <NAME1> = @A00250:26:H3YTWDSXX:
* <NAME2> =  1:N:0:ACTGTTCC
*/

// if set, the query name has lane
#define BIT_HAS_LANE (1<<0)
// if set, the query name has tile
#define BIT_HAS_TILE (1<<1)
// if set, the query name has X
#define BIT_HAS_X (1<<2)
// if set, the query name has Y
#define BIT_HAS_Y (1<<3)
// if set, the query name has NAME2
#define BIT_HAS_NAME2 (1<<4)
// if set, the data is paired end
#define BIT_PAIRED_END (1<<5)
// if set, the quality is encoded by col
#define BIT_ENCODE_PE_BY_OVERLAP (1<<6)
// if set, the quality is encoded by col
#define BIT_ENCODE_QUAL_BY_COL (1<<7)
// if set, the quality string will not be encoded
#define BIT_DONT_ENCODE_QUAL (1<<8)
// if set, the positions of N bases in the sequence will be encoded, which means the quality of N is not unique
#define BIT_ENCODE_N_POS (1<<9)

class RfqHeader{
public:
    RfqHeader();
    void read(istream& ifs);
    void write(ostream& ofs);
    bool hasLane();
    bool hasTile();
    bool hasX();
    bool hasY();
    bool hasName2();

    char qual2bit(char qual);
    char bit2qual(char qual);
    char majorQual();
    char nBaseQual();
    char nBaseBit();
    int majorQualNumBits();
    int normalQualNumBits();

    bool supportInterleaved();

    void makeQualityTable(vector<Read*>& reads, bool hasLaneTileXY);
    void setNBaseQual(char qual);

    uint8 qualBins();
    uint8* qualBuf();
    uint8 normalQualBins();
    uint8* normalQualBuf();

    bool identicalWith(RfqHeader* other);

private:
    void makeQualBitTable();
    void computeNormalQualBits();

public:
    // the flag should be "repaq"
    char mRepaqFlag[3];
    char mRepaqVersion[5];
    // to support backward compatibility
    char mAlgorithmVersion;
    uint8 mReadLengthBytes;
    uint16 mFlags;
    char mOverlapShift;

    /*
    * PE processing
    * if mSupportInterleaved == true, the chunk can be in BIT_PE_INTERLEAVED mode
    * read2name2 = read1name2;
    * read2name2[mName2DiffPos] = mName2DiffChar;
    */
    bool mSupportInterleaved;
    uint8 mName2DiffPos;
    char mName2DiffChar;

    // quality table
    uint8 mQualBins;
    uint8* mQualBuf;
    char mQual2BitTable[256];
    char mBit2QualTable[256];
    int mNormalQualNumBits;
    char mNBaseQual;

};

#endif