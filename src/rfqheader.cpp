#include "rfqheader.h"
#include <time.h>
#include "util.h"
#include <memory.h>
#include "endian.h"

RfqHeader::RfqHeader(){
    memset(this, 0, sizeof(RfqHeader));
    mRepaqFlag[0] = 'R';
    mRepaqFlag[1] = 'F';
    mRepaqFlag[2] = 'Q';
    mAlgorithmVersion = ALGORITHM_VER;
    memcpy(mRepaqVersion, VERSION_NUM, min(5, (int)sizeof(VERSION_NUM)));
    mReadLengthBytes = 1;
    mNBaseQual = '#';
    mOverlapShift = -24;
}

void RfqHeader::read(ifstream& ifs) {
    ifs.read(mRepaqFlag, 3);
    ifs.read(mRepaqVersion, 5);
    ifs.read(&mAlgorithmVersion, 1);
    if(ALGORITHM_VER != mAlgorithmVersion) {
        error_exit("The data is encoded by different version of repaq, please try repaq v" + string(mRepaqVersion, 5) + ". \nSee: https://github.com/OpenGene/repaq/releases");
    }
    ifs.read((char*)&mReadLengthBytes, 1);
    //ifs.read((char*)&mFlags, 2);
    mFlags = readLittleEndian16(ifs);
    ifs.read((char*)&mName2DiffPos, 1);
    ifs.read(&mName2DiffChar, 1);
    ifs.read(&mNBaseQual, 1);
    ifs.read(&mOverlapShift, 1);
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
    ofs.write(mRepaqVersion, 5);
    ofs.write(&mAlgorithmVersion, 1);
    ofs.write((const char*)&mReadLengthBytes, 1);
    //ofs.write((const char*)&mFlags, 2);
    writeLittleEndian(ofs, mFlags);
    ofs.write((const char*)&mName2DiffPos, 1);
    ofs.write(&mName2DiffChar, 1);
    ofs.write(&mNBaseQual, 1);
    ofs.write(&mOverlapShift, 1);
    ofs.write((const char*)&mQualBins, 1);
    ofs.write((const char*)mQualBuf, mQualBins);
}

void RfqHeader::setNBaseQual(char qual) {
    mNBaseQual = qual;
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
    // 0, 1, 3, 5, 7, 9...
    int maxQualVal = max(1, mQualBins*2 - 3);
    if(maxQualVal>=64) mNormalQualNumBits = 1;
    else if(maxQualVal>=32) mNormalQualNumBits = 2;
    else if(maxQualVal>=16) mNormalQualNumBits = 3;
    else if(maxQualVal>=8) mNormalQualNumBits = 4;
    else if(maxQualVal>=4) mNormalQualNumBits = 5;
    else if(maxQualVal>=2) mNormalQualNumBits = 6;
    else mNormalQualNumBits = 7;
}

void RfqHeader::makeQualityTable(vector<Read*>& reads, bool hasLaneTileXY) {
    int table[128];
    memset(table, 0, sizeof(int)*128);

    bool nBaseQualMade = false;
    bool sharpIsNQual = true;
    for(int r=0; r<reads.size(); r++) {
        Read* read = reads[r];
        const char* seq = read->mSeq.mStr.c_str();
        for(int i=0; i<read->length(); i++) {
            char q = read->mQuality[i];
            if(q < 0) {
                error_exit("bad quality value: " + to_string(q));
            }
            table[q]++;
            char base = seq[i];
            if(base == 'N') {
                if(!nBaseQualMade) {
                    mNBaseQual = q;
                    nBaseQualMade = true;
                } else if(mNBaseQual != q) {
                    //error_exit("The quality score of N bases are different.");
                    mFlags |= BIT_ENCODE_N_POS;
                    mNBaseQual = -1;
                }
            }
            if(base != 'A' && base != 'T' && base != 'C' && base != 'G' && base != 'N' ) {
                if(base == 'a' || base =='t' || base == 'c' || base == 't') {
                    string errmsg("repaq doesn't support FASTQ with lowercase bases (a/t/c/g)");
                    errmsg += "\nbut we get:\n" + read->mSeq.mStr;
                    error_exit(errmsg);
                }
                else {
                    string errmsg("repaq only supports FASTQ with uppercase bases (A/T/C/G/N)");
                    errmsg += "\nbut we get:\n" + read->mSeq.mStr;
                    error_exit(errmsg);
                }
            }
            if(q=='#' && base!='N')
                sharpIsNQual = false;
            if(nBaseQualMade && q==mNBaseQual && base!='N'){
                    //string errmsg("the quality " + string(1, q) + " should be with N base, but we get " + string(1, base));
                    //errmsg += "\nThe read is:\n" + read->mSeq.mStr + "\n" + read->mQuality;
                    //error_exit(errmsg);
                    // in this case, we have to encode N bases
                    mFlags |= BIT_ENCODE_N_POS;
                    mNBaseQual = -1;
                }

        }
    }

    // no N base found, and is not Illumina data
    if(!nBaseQualMade) {
        if(!sharpIsNQual || !hasLaneTileXY)
            mNBaseQual = 0;
    }

    mQualBins = 0;
    int maxNum = 0;
    char majorQual = '\0';

    bool hasN = false;

    for(int i=0; i<128; i++) {
        if(table[i] > 0) {
            mQualBins++;
            if(i == mNBaseQual)
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
        mQualBuf[mQualBins - 1] = mNBaseQual;

    // make the quality encoded by col for better compression ratio
    if(mQualBins <= 64)
        mFlags |= BIT_ENCODE_QUAL_BY_COL;

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
    return mNBaseQual;
}

char RfqHeader::nBaseBit() {
    return mQual2BitTable[mNBaseQual];
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

bool RfqHeader::supportInterleaved() {
    return mSupportInterleaved;
}

uint8 RfqHeader::qualBins() {
    return mQualBins;
}

uint8* RfqHeader::qualBuf() {
    return mQualBuf;
}

uint8 RfqHeader::normalQualBins() {
    return mQualBins - 1;
}

uint8* RfqHeader::normalQualBuf() {
    uint8 bins = normalQualBins();
    uint8* buf = new uint8[bins];
    int count = 0;
    for(int i=0; i<mQualBins; i++) {
        if(mQualBuf[i] != majorQual()) {
            buf[count] = mQualBuf[i];
            count++;
            if(count > bins)
                break;
        }
    }
    return buf;
}