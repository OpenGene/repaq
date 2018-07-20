#include "rfqheader.h"
#include <time.h>
#include "util.h"

RfqHeader::RfqHeader(){
    memset(this, 0, sizeof(RfqHeader));
    mRepaqFlag[0] = 'R';
    mRepaqFlag[1] = 'F';
    mRepaqFlag[2] = 'Q';
    mAlgorithmVersion = ALGORITHM_VER;
    mReadLengthBytes = 1;
}

void RfqHeader::read(ifstream& ifs) {
    ifs.read(mRepaqFlag, 3);
    ifs.read(&mAlgorithmVersion, 1);
    if(ALGORITHM_VER < mAlgorithmVersion) {
        error_exit("The software is too old to read this file, please update repaq. \nSee: https://github.com/OpenGene/repaq");
    }
    ifs.read((char*)&mReadLengthBytes, 1);
    ifs.read((char*)&mFlags, 2);
    ifs.read((char*)&mName2DiffPos, 1);
    ifs.read(&mName2DiffChar, 1);
    ifs.read((char*)&mQualBins, 1);

    mQualBuf = new uint8[mQualBins];
    ifs.read((char*)mQualBuf, mQualBins);

    makeQualBitTable();

    if(mRepaqFlag[0] != 'R' || mRepaqFlag[1] != 'F' || mRepaqFlag[2] != 'Q') {
        error_exit("Not a valid repaq file!");
    }
}

void RfqHeader::write(ofstream& ofs) {
    ofs.write(mRepaqFlag, 3);
    ofs.write(&mAlgorithmVersion, 1);
    ofs.write((const char*)&mReadLengthBytes, 1);
    ofs.write((const char*)&mFlags, 2);
    ofs.write((const char*)&mName2DiffPos, 1);
    ofs.write(&mName2DiffChar, 1);
    ofs.write((const char*)&mQualBins, 1);
    ofs.write((const char*)mQualBuf, mQualBins);
}

void RfqHeader::makeQualBitTable() {
    // qual bit = 0, 1, 3, 5, 7, 9, ...
    for(int i=0; i<mQualBins; i++){
        uint8 q = mQualBuf[i];
        int bit = i;
        if(i > 0)
            bit = 2*i - 1;
        mQual2BitTable[q] = bit;
        mBit2QualTable[bit] = q;
    }

    computeNormalQualBits();
}

void RfqHeader::computeNormalQualBits() {
    // except the major qual
    int normalQualBins = mQualBins - 1;
    if(normalQualBins>=64) normalQualBins = 7;
    else if(normalQualBins>=32) normalQualBins = 6;
    else if(normalQualBins>=16) normalQualBins = 5;
    else if(normalQualBins>=8) normalQualBins = 4;
    else if(normalQualBins>=4) normalQualBins = 3;
    else if(normalQualBins>=2) normalQualBins = 2;
    else normalQualBins = 1;

    mNormalQualNumBits = 8 - normalQualBins - 2;
}

void RfqHeader::makeQualityTable(string& qualStr) {
    int table[128];
    memset(table, 0, sizeof(int)*128);

    for(int i=0; i<qualStr.length(); i++) {
        char q = qualStr[i];
        if(q < 0) {
            error_exit("bad quality value: " + to_string(q));
        }
        table[q]++;
    }

    mQualBins = 0;
    int maxNum = 0;
    char majorQual = '\0';

    bool hasN = false;

    for(int i=0; i<128; i++) {
        if(table[i] > 0) {
            mQualBins++;
            if(i == '#')
                hasN = true;
        }

        if(table[i] > maxNum) {
            maxNum = table[i];
            majorQual = i;
        }

    }

    if(mQualBins == 0)
        error_exit("bad quality string, is this a valid FASTQ file?");
    else if(mQualBins >= 64) {
        cerr << "WARNING: this FASTQ file's quality bins are too complicated, which may affect the compression ratio." << endl;
        cerr << "Please confirm this is a valid FASTQ file." << endl;

        mFlags |= BIT_DONT_ENCODE_QUAL;
    }

    if(!hasN)
        mQualBins += 1;

    mQualBuf = new uint8[mQualBins];
    mQualBuf[0] = majorQual;
    int cur = 1;
    for(int i=0; i<128; i++) {
        if(i == majorQual)
            continue;
        if(table[i] > 0) {
            mQualBuf[cur] = i;
            cur++;
        }
    }

    if(!hasN)
        mQualBuf[mQualBins - 1] = '#';

    makeQualBitTable();
}

char RfqHeader::qual2bit(char qual) {
    return mQual2BitTable[qual];
}

char RfqHeader::bit2qual(char bit) {
    return mBit2QualTable[bit];
}

int RfqHeader::majorQualNumBits() {
    return 7;
}

int RfqHeader::normalQualNumBits() {
    return mNormalQualNumBits;
}


char RfqHeader::majorQual() {
    return mBit2QualTable[0];
}

char RfqHeader::nBaseQual() {
    return '#';
}

char RfqHeader::nBaseBit() {
    return mQual2BitTable['#'];
}

bool RfqHeader::hasLane() {
    return mFlags & BIT_HAS_LANE;
}

bool RfqHeader::hasTile() {
    return mFlags & BIT_HAS_TILE;
}

bool RfqHeader::hasX() {
    return mFlags & BIT_HAS_X;
}

bool RfqHeader::hasY() {
    return mFlags & BIT_HAS_Y;
}

bool RfqHeader::hasName2() {
    return mFlags & BIT_HAS_NAME2;
}

bool RfqHeader::hasLineBreakAtEnd() {
    return mFlags & BIT_HAS_LINE_BREAK_AT_END;
}

bool RfqHeader::hasLineBreakAtEndR2() {
    return mFlags & BIT_HAS_LINE_BREAK_AT_END_R2;
}

bool RfqHeader::supportInterleaved() {
    return mSupportInterleaved;
}