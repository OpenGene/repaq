#include "rfqcodec.h"
#include <time.h>
#include "util.h"
#include "fastqmeta.h"
#include <memory.h>
#include <sstream>

RfqCodec::RfqCodec(){
    mHeader = NULL;
}


void RfqCodec::setHeader(RfqHeader* header) {
    mHeader = header;
}

RfqHeader* RfqCodec::makeHeader(vector<Read*>& reads) {
    if(reads.size() == 0)
        return NULL;

    RfqHeader* header = new RfqHeader();
    bool hasLaneTileXY = true;
    int maxReadLen = 0;

    string qualStr;

    for(int i=0; i<reads.size(); i++) {
        Read* r = reads[i];
        FastqMeta meta = FastqMeta::parse(r->mName);
        hasLaneTileXY &= meta.hasLaneTileXY;
        maxReadLen = max(maxReadLen, r->length());
        qualStr += r->mQuality;
    }

    if(hasLaneTileXY) {
        header->mFlags |= BIT_HAS_LANE;
        header->mFlags |= BIT_HAS_TILE;
        header->mFlags |= BIT_HAS_X;
        header->mFlags |= BIT_HAS_Y;
        header->mFlags |= BIT_HAS_NAME2;
    }

    header->makeQualityTable(qualStr);

    if(maxReadLen>65535)
        header->mReadLengthBytes = 4;
    if(maxReadLen>255)
        header->mReadLengthBytes = 2;
    else
        header->mReadLengthBytes = 1;

    mHeader = header;
    return header;
}

RfqHeader* RfqCodec::makeHeader(vector<ReadPair*>& pairs) {
    if(pairs.size() == 0)
        return NULL;

    RfqHeader* header = new RfqHeader();
    bool hasLaneTileXY = true;
    int maxReadLen = 0;

    bool supportInterleaved = true;
    int name2DiffPos = 0;
    char name2DiffChar = '\0';

    string qualStr;

    for(int i=0; i<pairs.size(); i++) {
        Read* r1 = pairs[i]->mLeft;
        Read* r2 = pairs[i]->mRight;
        qualStr += r1->mQuality;
        qualStr += r2->mQuality;
        FastqMeta meta1 = FastqMeta::parse(r1->mName);
        FastqMeta meta2 = FastqMeta::parse(r2->mName);
        hasLaneTileXY &= meta1.hasLaneTileXY;
        hasLaneTileXY &= meta2.hasLaneTileXY;
        maxReadLen = max(maxReadLen, r1->length());
        maxReadLen = max(maxReadLen, r2->length());

        if(!hasLaneTileXY)
            supportInterleaved = false;
        else if(supportInterleaved){
            if(i == 0) {
                if(meta1.namePart2.length() != meta2.namePart2.length())
                    supportInterleaved = false;

                for(int p=0; p<meta1.namePart2.length(); p++) {
                    if(meta1.namePart2[p] != meta2.namePart2[p]) {
                        name2DiffPos = p;
                        name2DiffChar = meta2.namePart2[p];
                        break;
                    }
                }
            }

            if(meta1.namePart2.length() < name2DiffPos)
                supportInterleaved = false;
            else {
                // there is one, and just one different char
                if(name2DiffChar != '\0')
                    meta1.namePart2[name2DiffPos] = name2DiffChar;
                if(meta1.namePart2 !=  meta2.namePart2)
                    supportInterleaved = false;
            }
        }
    }

    if(supportInterleaved) {
        header->mSupportInterleaved = supportInterleaved;
        header->mName2DiffPos = name2DiffPos;
        header->mName2DiffChar = name2DiffChar;
    }

    header->makeQualityTable(qualStr);

    if(hasLaneTileXY) {
        header->mFlags |= BIT_HAS_LANE;
        header->mFlags |= BIT_HAS_TILE;
        header->mFlags |= BIT_HAS_X;
        header->mFlags |= BIT_HAS_Y;
        header->mFlags |= BIT_HAS_NAME2;
    }

    header->mFlags |= BIT_PAIRED_END;

    if(maxReadLen>65535)
        header->mReadLengthBytes = 4;
    if(maxReadLen>255)
        header->mReadLengthBytes = 2;
    else
        header->mReadLengthBytes = 1;

    mHeader = header;
    return header;
}

RfqChunk* RfqCodec::encodeChunk(vector<ReadPair*>& pairs) {
    if(mHeader == NULL)
        makeHeader(pairs);

    if(mHeader == NULL)
        return NULL;
    
    vector<Read*> reads;
    for(int i=0; i<pairs.size(); i++) {
        ReadPair* pair = pairs[i];
        reads.push_back(pair->mLeft);
        reads.push_back(pair->mRight);
    }
    return encodeChunk(reads, true);
}

RfqChunk* RfqCodec::encodeChunk(vector<Read*>& reads, bool isPE) {
    int s = reads.size();
    if(s == 0)
        return NULL;

    if(mHeader == NULL)
        makeHeader(reads);

    bool readLenSame = true;
    bool name1LenSame = true;
    bool name2LenSame = true;
    bool strandLenSame = true;
    bool strandSame = true;
    bool laneSame = true;
    bool tileSame = true;
    bool name1Same = true;
    bool name2Same = true;

    Read* r0 = reads[0];
    FastqMeta meta0 = FastqMeta::parse(r0->mName);

    int readLen0 = r0->length();
    int name1Len0 = meta0.namePart1.length();
    int name2Len0 = meta0.namePart2.length();
    int strandLen0 = r0->mStrand.length();
    string strand0 = r0->mStrand;
    uint8 lane0 = meta0.lane;
    uint16 tile0 = meta0.tile;
    string name10 = meta0.namePart1;
    string name20 = meta0.namePart2;

    uint32 totalReadLen = 0;
    uint32 totalName1Len = 0;
    uint32 totalName2Len = 0;
    uint32 totalStrandLen = 0;

    uint8* laneBuf = new uint8[s];
    memset(laneBuf, 0, sizeof(uint8)*s);
    uint16* tileBuf = new uint16[s];
    memset(tileBuf, 0, sizeof(uint16)*s);
    uint16* xBuf = new uint16[s];
    memset(xBuf, 0, sizeof(uint16)*s);
    uint16* yBuf = new uint16[s];
    memset(yBuf, 0, sizeof(uint16)*s);

    // in this BIT_PE_INTERLEAVED mode
    // the reads are paired-end, interleaved, read2 sequence is reversed and complemented
    // tile, X and Y information are only stored once in odd read1
    // name2 are only stored once in read1, and name2 can be got by transfering read1 to read2 by replace 1 with 2 
    bool canBePeInterleaved = isPE && mHeader->supportInterleaved();
    string lastName2;
    uint16 lastX;
    uint16 lastY;
    uint16 lastTile;
    uint8 lastLane;
    for(int i=0; i<reads.size(); i++) {
        Read* r = reads[i];
        int rlen = r->length();
        FastqMeta meta = FastqMeta::parse(r->mName);

        readLenSame &= readLen0 == rlen;
        name1LenSame &= name1Len0 == meta.namePart1.length();
        name2LenSame &= name2Len0 == meta.namePart2.length();
        strandLenSame &= strandLen0 == r->mStrand.length();
        strandSame  &= strand0 == r->mStrand;
        laneSame  &= lane0 == meta.lane;
        tileSame  &= tile0 == meta.tile;
        name1Same  &= name10 == meta.namePart1;
        if(!canBePeInterleaved)
            name2Same  &= name20 == meta.namePart2;
        else {
            // read2, check its consistent with read1
            if(i%2 == 1) {
                string replacedName2 = lastName2;
                // !='\0' means it has only one difference
                if(mHeader->mName2DiffChar != '\0')
                    replacedName2[mHeader->mName2DiffPos] = mHeader->mName2DiffChar;
                if(replacedName2 != meta.namePart2) {
                    canBePeInterleaved = false;
                    name2Same  &= name20 == meta.namePart2;
                }
            } else {
                lastName2 = meta.namePart2;
            }
        }

        laneBuf[i] = meta.lane;
        tileBuf[i] = meta.tile;
        xBuf[i] = meta.x;
        yBuf[i] = meta.y;

        if(canBePeInterleaved) {
            // read2, check its consistent with read1
            if(i%2 == 1) {
                canBePeInterleaved &= lastLane == meta.lane;
                canBePeInterleaved &= lastTile == meta.tile;
                canBePeInterleaved &= lastX == meta.x;
                canBePeInterleaved &= lastY == meta.y;
            } else {
                lastLane = meta.lane;
                lastTile = meta.tile;
                lastX = meta.x;
                lastY = meta.y;
            }
        }

        totalReadLen += rlen;
        totalName1Len += meta.namePart1.length();
        totalName2Len += meta.namePart2.length();
        totalStrandLen += r->mStrand.length();
    }


    if(canBePeInterleaved) {
        for(int p=0; p<s/2; p++) {
            // only keep read1
            laneBuf[p] = laneBuf[p*2];
            tileBuf[p] = tileBuf[p*2];
            xBuf[p] = xBuf[p*2];
            yBuf[p] = yBuf[p*2];
        }
    }

    uint8* readLenBuf = NULL;
    if(!readLenSame) 
        readLenBuf = new uint8[mHeader->mReadLengthBytes * s];
    uint16* readLenBuf16 = (uint16*)readLenBuf;
    uint32* readLenBuf32 = (uint32*)readLenBuf;

    uint8* name1LenBuf = NULL;
    if(!name1LenSame)
        name1LenBuf = new uint8[s];
    uint8* name2LenBuf = NULL;
    if(!name2LenSame)
        name2LenBuf = new uint8[s];
    uint8* strandLenBuf = NULL;
    if(!strandLenSame)
        strandLenBuf = new uint8[s];

    char* name1Buf = NULL;
    if(!name1Same)
        name1Buf = new char[totalName1Len];
    char* name2Buf = NULL;
    if(!name2Same)
        name2Buf = new char[totalName2Len];
    char* strandBuf = NULL;
    if(!strandSame)
        strandBuf = new char[totalStrandLen];

    char* seqBufOriginal = new char[totalReadLen];
    memset(seqBufOriginal, 0, totalReadLen);
    uint8* qualBufOriginal = new uint8[totalReadLen];
    memset(qualBufOriginal, 0, totalReadLen);

    int name1Copied = 0;
    int name2Copied = 0;
    int strandCopied = 0;
    int seqCopied = 0;

    for(int i=0; i<reads.size(); i++) {
        Read* r = reads[i];
        int rlen = r->length();

        if(!readLenSame) {
            if(mHeader->mReadLengthBytes == 1)
                readLenBuf[i] = rlen;
            else if(mHeader->mReadLengthBytes == 2)
                readLenBuf16[i] = rlen;
            else
                readLenBuf32[i] = rlen;
        }

        if(!name1Same || !name2Same) {
            FastqMeta meta = FastqMeta::parse(r->mName);
            if(!name1Same) {
                int name1len = meta.namePart1.length();
                memcpy(name1Buf + name1Copied, meta.namePart1.c_str(), name1len);
                name1Copied += name1len;
                if(!name1LenSame)
                    name1LenBuf[i] = name1len;
            }
            if(!name2Same) {
                int name2len = meta.namePart2.length();
                memcpy(name2Buf + name2Copied, meta.namePart2.c_str(), name2len);
                name2Copied += name2len;
                if(!name2LenSame)
                    name2LenBuf[i] = name2len;
            }
        }

        if(!strandSame) {
            int strandlen = r->mStrand.length();
            memcpy(strandBuf + strandCopied, r->mStrand.c_str(), strandlen);
            strandCopied += strandlen;
            if(!strandLenSame)
                strandLenBuf[i] = strandlen;
        }

        if(canBePeInterleaved)
            r->changeToReverseComplement();

        memcpy(seqBufOriginal + seqCopied, r->mSeq.mStr.c_str(), rlen);
        memcpy(qualBufOriginal + seqCopied, r->mQuality.c_str(), rlen);
        seqCopied += rlen;
    }

    int encodedSeqBufLen = (totalReadLen + 3) / 4;
    char* seqBufEncoded = new char[encodedSeqBufLen];
    memset(seqBufEncoded, 0, encodedSeqBufLen);
    char* qualBufEncoded = new char[totalReadLen];
    memset(qualBufEncoded, 0, totalReadLen);

    uint32 encodedQualBufLen = encodeSeqQual(seqBufOriginal, qualBufOriginal, seqBufEncoded, qualBufEncoded, totalReadLen);

    delete seqBufOriginal;
    seqBufOriginal = NULL;
    delete qualBufOriginal;
    qualBufOriginal = NULL;

    RfqChunk* chunk = new RfqChunk(mHeader);

    chunk->mReads = reads.size();

    if(canBePeInterleaved)
        chunk->mFlags |= BIT_PE_INTERLEAVED;

    if(readLenSame) chunk->mFlags |= BIT_READ_LEN_SAME;
    if(name1LenSame) chunk->mFlags |= BIT_NAME1_LEN_SAME;
    if(name2LenSame) chunk->mFlags |= BIT_NAME2_LEN_SAME;
    if(strandLenSame) chunk->mFlags |= BIT_STRAND_LEN_SAME;
    if(strandSame) chunk->mFlags |= BIT_STRAND_SAME;
    if(laneSame) chunk->mFlags |= BIT_LANE_SAME;
    if(tileSame) chunk->mFlags |= BIT_TILE_SAME;
    if(name1Same) chunk->mFlags |= BIT_NAME1_SAME;
    if(name2Same) chunk->mFlags |= BIT_NAME2_SAME;

    chunk->mSeqBufSize = encodedSeqBufLen;
    chunk->mQualBufSize = encodedQualBufLen;

    if(readLenSame) {
        chunk->mReadLenBuf = new uint8[mHeader->mReadLengthBytes];
        memcpy(chunk->mReadLenBuf, &readLen0, mHeader->mReadLengthBytes);
        chunk->mReadLenBufSize = mHeader->mReadLengthBytes;
    } else {
        chunk->mReadLenBuf = readLenBuf;
        chunk->mReadLenBufSize = mHeader->mReadLengthBytes * s;
    }

    if(name1LenSame) {
        chunk->mName1LenBuf = new uint8[1];
        chunk->mName1LenBuf[0] = name1Len0;
        chunk->mName1LenBufSize = 1;
    } else {
        chunk->mName1LenBuf = name1LenBuf;
        chunk->mName1LenBufSize = s;
    }

    if(name2LenSame) {
        chunk->mName2LenBuf = new uint8[1];
        chunk->mName2LenBuf[0] = name2Len0;
        chunk->mName2LenBufSize = 1;
    } else {
        chunk->mName2LenBuf = name2LenBuf;
        chunk->mName2LenBufSize = s;
    }

    if(strandLenSame) {
        chunk->mStrandLenBuf = new uint8[1];
        chunk->mStrandLenBuf[0] = strandLen0;
        chunk->mStrandLenBufSize = 1;
    } else {
        chunk->mStrandLenBuf = strandLenBuf;
        chunk->mStrandLenBufSize = s;
    }

    if(laneSame) {
        chunk->mLaneBuf = new uint8[1];
        chunk->mLaneBuf[0] = lane0;
        delete[] laneBuf;
        laneBuf = NULL;
        chunk->mLaneBufSize = 1;
    } else {
        chunk->mLaneBuf = laneBuf;
        if(canBePeInterleaved)
            chunk->mLaneBufSize = s/2;
        else
            chunk->mLaneBufSize = s;
    }

    if(tileSame) {
        chunk->mTileBuf = new uint16[1];
        chunk->mTileBuf[0] = tile0;
        delete[] tileBuf;
        tileBuf = NULL;
        chunk->mLaneBufSize = sizeof(uint16);
    } else {
        chunk->mTileBuf = tileBuf;
        if(canBePeInterleaved)
            chunk->mLaneBufSize = sizeof(uint16) * s/2;
        else
            chunk->mLaneBufSize = sizeof(uint16) * s;
    }

    // TODO
    chunk->mXBuf = xBuf;
    chunk->mYBuf = yBuf;

    if(name1Same) {
        chunk->mName1Buf = new char[name1Len0];
        memcpy(chunk->mName1Buf, name10.c_str(), name1Len0);
        chunk->mName1BufSize = name1Len0;
    } else {
        chunk->mName1Buf = name1Buf;
        chunk->mName1BufSize = totalName1Len;
    }

    if(name2Same) {
        chunk->mName2Buf = new char[name2Len0];
        memcpy(chunk->mName2Buf, name20.c_str(), name2Len0);
        chunk->mName2BufSize = name2Len0;
    } else {
        chunk->mName2Buf = name2Buf;
        chunk->mName2BufSize = totalName2Len;
    }

    if(strandSame) {
        chunk->mStrandBuf = new char[strandLen0];
        memcpy(chunk->mStrandBuf, strand0.c_str(), strandLen0);
        chunk->mStrandBufSize = strandLen0;
    } else {
        chunk->mStrandBuf = strandBuf;
        chunk->mStrandBufSize = totalStrandLen;
    }

    chunk->mSeqBuf = seqBufEncoded;

    chunk->mQualBuf = new uint8[encodedQualBufLen];
    memcpy(chunk->mQualBuf, qualBufEncoded, encodedQualBufLen);
    delete qualBufEncoded;
    qualBufEncoded = NULL;

    chunk->calcTotalBufSize();

    return chunk;
}

uint32 RfqCodec::encodeSeqQual(char* seq, uint8* qual, char* seqEncoded, char* qualEncoded, uint32 totalLen) {
    // encode seq first
    for(int i=0; i<totalLen; i++) {
        char c = seq[i];
        uint8 val = 0;
        switch(c) {
            case 'G': val=0; break;
            case 'A': val=1; break;
            case 'T': val=2; break;
            case 'C': val=3; break;
            default: break;
        }
        int pos = i>>2;
        int offset = i & 0x03;
        val = val << (offset*2);
        seqEncoded[pos] |= val;

        // mask N base in quality string
        /*if(c == 'N') {
            qual[i] = -1;
        }*/
    }

    // dont encode qual
    if(mHeader->mFlags & BIT_DONT_ENCODE_QUAL) {
        memcpy(qualEncoded, qual, totalLen);
        return totalLen;
    }

    // encode quality
    uint32 qualBufLen = 0;
    char mq = mHeader->majorQual();
    int mqNumBits = mHeader->majorQualNumBits();
    int nqNumBits = mHeader->normalQualNumBits();
    int mqNumMax = (0x01 << mqNumBits);
    int nqNumMax = (0x01 << nqNumBits);
    char curQual = qual[0];
    int first = 0;
    int i=1;
    while(i < totalLen) {
        char q = qual[i];

        bool needRestart = false;
        if(q != curQual)
            needRestart = true;
        else {
            if(curQual == mq && i-first >= mqNumMax)
                needRestart = true;
            if(curQual != mq && i-first >= nqNumMax)
                needRestart = true;
        }

        if(needRestart) {
            // num = 0 means 1 continue
            uint8 num = i - first - 1;
            char qualBit = mHeader->qual2bit(curQual);
            uint8 enQual;
            if(curQual == mq) {
                enQual = qualBit | (num << (8 - mqNumBits));
            } else {
                enQual = qualBit | (num << (8 - nqNumBits));
            }
            qualEncoded[qualBufLen] = enQual;
            qualBufLen++;

            first = i;
            curQual = q;
        }

        i++;
    }
    // encode last one
    // num = 0 means 1 continue
    uint8 num = totalLen - first - 1;
    char qualBit = mHeader->qual2bit(curQual);
    uint8 enQual;
    if(curQual == mq) {
        enQual = qualBit | (num << (8 - mqNumBits));
    } else {
        enQual = qualBit | (num << (8 - nqNumBits));
    }
    qualEncoded[qualBufLen] = enQual;
    qualBufLen++;

    return qualBufLen;
}

void RfqCodec::decodeSeqQual(RfqChunk* chunk, string& seq, string& qual, uint32 len) {
    int mqNumBits = mHeader->majorQualNumBits();
    int nqNumBits = mHeader->normalQualNumBits();
    char nBaseQual = mHeader->nBaseQual();
    char nBaseBit = mHeader->nBaseBit();
    char nqMask = 0;
    for(int b=0; b<(8-nqNumBits); b++)
        nqMask |= (0x01 << b);
    
    uint32 decoded = 0;
    // decode sequence
    for(int i=0; i<chunk->mSeqBufSize; i++) {
        char seqEncoded = chunk->mSeqBuf[i];
        for(int b=0; b<4; b++) {
            char basebits = (seqEncoded & (0x03<<(b*2))) >> (b*2);
            char base = 'N';
            switch(basebits) {
                case 0: base = 'G'; break;
                case 1: base = 'A'; break;
                case 2: base = 'T'; break;
                case 3: base = 'C'; break;
                default: break;
            }
            if(decoded < len) {
                seq[decoded] = base;
                decoded++;
            } else 
                break;
        }
        if(decoded >= len)
            break;
    }

    // qual is not encoded
    if(mHeader->mFlags & BIT_DONT_ENCODE_QUAL) {
        for(int i=0; i<chunk->mQualBufSize; i++) {
            qual[i] = chunk->mQualBuf[i];
            if(qual[i] == nBaseQual)
                seq[i] = 'N';
        }
        return;
    }

    decoded = 0;
    while(decoded < len) {
        for(int i=0; i<chunk->mQualBufSize; i++) {
            uint8 qualEncoded = chunk->mQualBuf[i];
            char q = 0;
            uint8 num = 0;
            // major qual
            if( (qualEncoded & 0x01) == 0) {
                q = 0;
                num = qualEncoded >> (8 - mqNumBits);
            } else {
                q = qualEncoded & nqMask;
                num = qualEncoded >> (8 - nqNumBits);
            }
            // num = 0 means 1 continue, so we increase it by 1 here
            num += 1;
            char qualVal = mHeader->bit2qual(q);
            // this is a N base, restore its quality
            if(q == -1)
                qualVal = nBaseQual;
            for(int fill = decoded; fill < decoded + num && fill < len; fill++) {
                qual[fill] = qualVal;
                // write N in sequence
                if(q == nBaseBit)
                    seq[fill] = 'N';
            }
            decoded += num;
            if(decoded >= len)
                break;
        }
    }
}

vector<Read*> RfqCodec::decodeChunk(RfqChunk* chunk) {
    vector<Read*> ret;

    if(mHeader == NULL)
        return ret;

    bool peInterleaved = chunk->mFlags & BIT_PE_INTERLEAVED;

    uint32 readLen0 = 0;
    uint16* rlenBuf16 = (uint16*)chunk->mReadLenBuf;
    uint32* rlenBuf32 = (uint32*)chunk->mReadLenBuf;
    switch(mHeader->mReadLengthBytes) {
        case 1: readLen0 = chunk->mReadLenBuf[0]; break;
        case 2: readLen0 = rlenBuf16[0]; break;
        case 4: readLen0 = rlenBuf32[0]; break;
        default: error_exit("header incorrect: read length bytes should be 1/2/4");
    }

    uint32 seqLen = 0;
    if(chunk->mFlags & BIT_READ_LEN_SAME)
        seqLen = readLen0 * chunk->mReads;
    else {
        for(int i=0; i<chunk->mReads; i++) {
            switch(mHeader->mReadLengthBytes) {
                case 1: seqLen += chunk->mReadLenBuf[i]; break;
                case 2: seqLen += rlenBuf16[i]; break;
                case 4: seqLen += rlenBuf32[i]; break;
                default: error_exit("header incorrect: read length bytes should be 1/2/4");
            }
        }
    }

    string allSeq(seqLen, 'N');
    string allQual(seqLen, mHeader->nBaseQual());

    decodeSeqQual(chunk, allSeq, allQual, seqLen);
    int name1Len0 = chunk->mName1LenBuf[0];
    string name10(chunk->mName1Buf, name1Len0);
    int strandLen0 = chunk->mStrandLenBuf[0];
    string strand0(chunk->mStrandBuf, strandLen0);

    int name2Len0 = 0;
    string name20;
    uint8 lane0 = 0;
    uint16 tile0 = 0;

    char* curName1 = chunk->mName1Buf;
    char* curName2 = chunk->mName2Buf;
    char* curStrand = chunk->mStrandBuf;
    uint32 curSeq = 0;

    if(mHeader->hasName2()) {
        name2Len0 = chunk->mName2LenBuf[0];
        name20 = string(chunk->mName2Buf, name2Len0);
    }

    if(mHeader->hasLane())
        lane0 = chunk->mLaneBuf[0];
    if(mHeader->hasTile())
        tile0 = chunk->mTileBuf[0];

    for(int r=0; r<chunk->mReads; r++) {
        uint32 rlen = readLen0;
        if( (chunk->mFlags & BIT_READ_LEN_SAME) == false) {
            switch(mHeader->mReadLengthBytes) {
                case 1: rlen = chunk->mReadLenBuf[r]; break;
                case 2: rlen = rlenBuf16[r]; break;
                case 4: rlen = rlenBuf32[r]; break;
                default: error_exit("header incorrect: read length bytes should be 1/2/4");
            }
        }

        string sequence = allSeq.substr(curSeq, rlen);
        string quality = allQual.substr(curSeq, rlen);
        curSeq += rlen;

        // get name;
        stringstream ss;
        // get name1
        string name1;
        if(chunk->mFlags & BIT_NAME1_SAME)
            name1 = name10;
        else if(chunk->mFlags & BIT_NAME1_LEN_SAME) {
            name1 = string(curName1, name1Len0);
            curName1 += name1Len0;
        }
        else {
            int name1Len = chunk->mName1LenBuf[r];
            name1 = string(curName1, name1Len);
            curName1 += name1Len;
        }

        ss << name1;

        int xyPos = r;
        if(peInterleaved)
            xyPos = r/2;

        if(mHeader->hasLane()) {
            uint8 lane = lane0;
            if((chunk->mFlags & BIT_LANE_SAME) == false)
                lane = chunk->mLaneBuf[xyPos];

            ss << ":" << (int)lane;
        }

        if(mHeader->hasTile()) {
            uint16 tile = tile0;
            if((chunk->mFlags & BIT_TILE_SAME) == false)
                tile = chunk->mTileBuf[xyPos];

            ss << ":" << tile;
        }

        if(mHeader->hasX()) {
            uint16 x = chunk->mXBuf[xyPos];
            ss << ":" << x;
        }

        if(mHeader->hasY()) {
            uint16 y = chunk->mYBuf[xyPos];
            ss << ":" << y;
        }

        if(mHeader->hasName2()) {
            // get name2
            string name2;
            if(chunk->mFlags & BIT_NAME2_SAME) {
                name2 = name20;
                if(peInterleaved) {
                    // is read2
                    if(r % 2 == 1) {
                        // !='\0' means it has only one difference
                        if(mHeader->mName2DiffChar != '\0')
                            name2[mHeader->mName2DiffPos] = mHeader->mName2DiffChar;
                    }
                }
            }
            else if(chunk->mFlags & BIT_NAME2_LEN_SAME) {
                name2 = string(curName2, name2Len0);
                curName2 += name2Len0;
            }
            else {
                int name2Len = chunk->mName2LenBuf[r];
                name2 = string(curName2, name2Len);
                curName2 += name2Len;
            }

            ss << name2;
        }

        string name = ss.str();

        // get strand
        string strand;
        if(chunk->mFlags & BIT_STRAND_SAME)
            strand = strand0;
        else if(chunk->mFlags & BIT_STRAND_LEN_SAME) {
            strand = string(curStrand, strandLen0);
            curStrand += strandLen0;
        }
        else {
            int strandLen = chunk->mStrandLenBuf[r];
            strand = string(curStrand, strandLen);
            curStrand += strandLen;
        }

        Read* read = new Read(name, sequence, strand, quality);
        if(peInterleaved)
            read->changeToReverseComplement();
        ret.push_back(read);
    }

    return ret;
}