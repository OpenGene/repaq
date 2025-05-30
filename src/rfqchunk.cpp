#include "rfqchunk.h"
#include <memory.h>
#include "endian.h"

RfqChunk::RfqChunk(RfqHeader* header){
    memset(this, 0, sizeof(RfqChunk));
    mHeader = header;
}

RfqChunk::~RfqChunk(){
    if(mReadLenBuf)
        delete[] mReadLenBuf;
    if(mName1LenBuf)
        delete[] mName1LenBuf;
    if(mName2LenBuf)
        delete[] mName2LenBuf;
    if(mLaneBuf)
        delete[] mLaneBuf;
    if(mTileBuf)
        delete[] mTileBuf;
    if(mXBuf)
        delete[] mXBuf;
    if(mYBuf)
        delete[] mYBuf;
    if(mName1Buf)
        delete[] mName1Buf;
    if(mName2Buf)
        delete[] mName2Buf;
    if(mSeqBuf)
        delete[] mSeqBuf;
    if(mQualBuf)
        delete[] mQualBuf;
    if(mStrandBuf)
        delete[] mStrandBuf;
    if(mOverlapBuf)
        delete[] mOverlapBuf;
    if(mNPosBuf)
        delete[] mNPosBuf;
}

void RfqChunk::readReadLenBuf(istream& ifs) {
    // read length
    int readLenCount = 1;
    if((mFlags & BIT_READ_LEN_SAME) == false) {
        readLenCount = mReads;
    }
    int bytes = mHeader->mReadLengthBytes;
    mReadLenBuf = new uint8[readLenCount * bytes];
    ifs.read((char*)mReadLenBuf, readLenCount * bytes);
    mReadLenBufSize = 0;
    for(int i=0; i<readLenCount; i++) {
        uint8* data = mReadLenBuf + i*bytes;
        if(bytes == 1) {
            mReadLenBufSize += *data;
        } else if(bytes == 2) {
            mReadLenBufSize += *((uint16*)data);
        } else {
            mReadLenBufSize += *((uint32*)data);
        }
    }
}

void RfqChunk::readName1LenBuf(istream& ifs) {
    // name1 length
    mName1LenBufSize = 1;
    if((mFlags & BIT_NAME1_LEN_SAME) == false) {
        mName1LenBufSize = mReads;
    }
    mName1LenBuf = new uint8[mName1LenBufSize];
    ifs.read((char*)mName1LenBuf, mName1LenBufSize);
    mName1BufSize = 0;
    for(int i=0; i<mName1LenBufSize; i++) {
        mName1BufSize += mName1LenBuf[i];
    }
    if( (mFlags & BIT_NAME1_LEN_SAME) &&  (mFlags & BIT_NAME1_SAME)==false)
        mName1BufSize *= mReads;
}

void RfqChunk::readName2LenBuf(istream& ifs) {
    // name2 length
    mName2LenBufSize = 1;
    if((mFlags & BIT_NAME2_LEN_SAME) == false) {
        mName2LenBufSize = mReads;
    }
    mName2LenBuf = new uint8[mName2LenBufSize];
    ifs.read((char*)mName2LenBuf, mName2LenBufSize);
    mName2BufSize = 0;
    for(int i=0; i<mName2LenBufSize; i++) {
        mName2BufSize += mName2LenBuf[i];
    }
    if( (mFlags & BIT_NAME2_LEN_SAME) &&  (mFlags & BIT_NAME2_SAME)==false)
        mName2BufSize *= mReads;
}

void RfqChunk::readStrandLenBuf(istream& ifs) {
    // strand length
    mStrandLenBufSize = 1;
    if((mFlags & BIT_STRAND_LEN_SAME) == false) {
        mStrandLenBufSize = mReads;
    }
    mStrandLenBuf = new uint8[mStrandLenBufSize];
    ifs.read((char*)mStrandLenBuf, mStrandLenBufSize);
    mStrandBufSize = 0;
    for(int i=0; i<mStrandLenBufSize; i++) {
        mStrandBufSize += mStrandLenBuf[i];
    }
    if( (mFlags & BIT_STRAND_LEN_SAME) &&  (mFlags & BIT_STRAND_SAME)==false)
        mStrandBufSize *= mReads;
}

void RfqChunk::readLaneBuf(istream& ifs) {
    int laneCount = 1;
    if((mFlags & BIT_LANE_SAME) == false) {
        if(mFlags & BIT_PE_INTERLEAVED)
            laneCount = mReads / 2;
        else
            laneCount = mReads;
    }
    mLaneBuf = new uint8[laneCount];
    ifs.read((char*)mLaneBuf, sizeof(uint8)*laneCount);
}

void RfqChunk::readTileBuf(istream& ifs) {
    int tileCount = 1;
    if((mFlags & BIT_TILE_SAME) == false) {
        if(mFlags & BIT_PE_INTERLEAVED)
            tileCount = mReads / 2;
        else
            tileCount = mReads;
    }
    mTileBuf = new uint16[tileCount];
    ifs.read((char*)mTileBuf, sizeof(uint16)*tileCount);
    // convert from little endian if the system is big endian
    if(!isLittleEndian()) {
        for(int i=0; i<tileCount; i++) {
            mTileBuf[i] = adaptToLittleEndian(mTileBuf[i]);
        }
    }
}

void RfqChunk::calcTotalBufSize() {
    mSize = sizeof(mSize) + sizeof(mReads) + sizeof(mFlags) + sizeof(mSeqBufSize) + sizeof(mQualBufSize);
    mSize +=  mReadLenBufSize + mName1LenBufSize + mName2LenBufSize + mStrandLenBufSize;
    mSize += mLaneBufSize + mTileBufSize + mName1BufSize + mName2BufSize + mStrandBufSize;
    mSize += mSeqBufSize + mQualBufSize;
    // overlap buf size;
    if( (mFlags & BIT_PE_INTERLEAVED) && (mHeader->mFlags & BIT_ENCODE_PE_BY_OVERLAP))
        mSize += mReads/2;
    if(mHeader->encodeNPos()) {
        mSize += sizeof(mNPosBufSize);
        mSize += mNPosBufSize;
    }
    if(mHeader->hasX()) {
        mSize += sizeof(mXBufSize) + mXBufSize;
    }
    if(mHeader->hasY()) {
        mSize += sizeof(mYBufSize) + mYBufSize;
    }
}

void RfqChunk::read(istream& ifs) {
    //ifs.read((char*)&mSize, sizeof(uint32));
    mSize = readLittleEndian32(ifs);
    //ifs.read((char*)&mReads, sizeof(uint32));
    mReads = readLittleEndian32(ifs);
    //ifs.read((char*)&mFlags, sizeof(uint16));
    mFlags = readLittleEndian16(ifs);
    //ifs.read((char*)&mSeqBufSize, sizeof(uint32));
    mSeqBufSize = readLittleEndian32(ifs);
    //ifs.read((char*)&mQualBufSize, sizeof(uint32));
    mQualBufSize = readLittleEndian32(ifs);
    // mNPosBufSize
    if(mHeader->encodeNPos())
        mNPosBufSize = readLittleEndian32(ifs);

    readReadLenBuf(ifs);
    readName1LenBuf(ifs);
    if(mHeader->hasName2())
        readName2LenBuf(ifs);

    readStrandLenBuf(ifs);

    if(mHeader->hasLane())
        readLaneBuf(ifs);
    if(mHeader->hasTile())
        readTileBuf(ifs);

    int xyCount = mReads;
    if(mFlags & BIT_PE_INTERLEAVED)
        xyCount = mReads / 2;
    if(mHeader->hasX()) {
        mXBufSize = readLittleEndian32(ifs);
        mXBuf = new uint8[mXBufSize];
        ifs.read((char*)mXBuf, mXBufSize);
    }
    if(mHeader->hasY()) {
        mYBufSize = readLittleEndian32(ifs);
        mYBuf = new uint8[mYBufSize];
        ifs.read((char*)mYBuf, mYBufSize);
    }

    mName1Buf = new char[mName1BufSize];
    ifs.read(mName1Buf, mName1BufSize);

    if(mHeader->hasName2()) {
        mName2Buf = new char[mName2BufSize];
        ifs.read(mName2Buf, mName2BufSize);
    }

    mStrandBuf = new char[mStrandBufSize];
    ifs.read(mStrandBuf, mStrandBufSize);

    mSeqBuf = new char[mSeqBufSize];
    ifs.read(mSeqBuf, mSeqBufSize);

    mQualBuf = new uint8[mQualBufSize];
    ifs.read((char*)mQualBuf, mQualBufSize);

    if( (mFlags & BIT_PE_INTERLEAVED) && (mHeader->mFlags & BIT_ENCODE_PE_BY_OVERLAP)) {
        mOverlapBuf = new char[mReads/2];
        ifs.read(mOverlapBuf, mReads/2);
    }

    if(mHeader->encodeNPos()) {
        mNPosBuf = new uint8[mNPosBufSize];
        ifs.read((char*)mNPosBuf, mNPosBufSize);
    }
}

void RfqChunk::write(ostream& ofs) {
    //ofs.write((const char*)&mSize, sizeof(uint32));
    writeLittleEndian(ofs, mSize);
    //ofs.write((const char*)&mReads, sizeof(uint32));
    writeLittleEndian(ofs, mReads);
    //ofs.write((char*)&mFlags, sizeof(uint16));
    writeLittleEndian(ofs, mFlags);
    //ofs.write((const char*)&mSeqBufSize, sizeof(uint32));
    writeLittleEndian(ofs, mSeqBufSize);
    //ofs.write((const char*)&mQualBufSize, sizeof(uint32));
    writeLittleEndian(ofs, mQualBufSize);
    if(mHeader->encodeNPos())
        writeLittleEndian(ofs, mNPosBufSize);

    ofs.write((const char*)mReadLenBuf, mReadLenBufSize);
    ofs.write((const char*)mName1LenBuf, mName1LenBufSize);

    if(mHeader->hasName2())
        ofs.write((const char*)mName2LenBuf, mName2LenBufSize);

    ofs.write((const char*)mStrandLenBuf, mStrandLenBufSize);

    if(mHeader->hasLane()) {
        int laneCount = 1;
        if((mFlags & BIT_LANE_SAME) == false) {
            if(mFlags & BIT_PE_INTERLEAVED)
                laneCount = mReads/2;
            else
                laneCount = mReads;
        }
        ofs.write((const char*)mLaneBuf, sizeof(uint8)*laneCount);
    }

    if(mHeader->hasTile()) {
        int tileCount = 1;
        if((mFlags & BIT_TILE_SAME) == false) {
            if(mFlags & BIT_PE_INTERLEAVED)
                tileCount = mReads/2;
            else
                tileCount = mReads;
        }
        // convert to little endian if the system is big endian
        if(!isLittleEndian()) {
            for(int i=0; i<tileCount; i++) {
                mTileBuf[i] = adaptToLittleEndian(mTileBuf[i]);
            }
        }
        ofs.write((const char*)mTileBuf, sizeof(uint16)*tileCount);
        // then restore it back
        if(!isLittleEndian()) {
            for(int i=0; i<tileCount; i++) {
                mTileBuf[i] = adaptToLittleEndian(mTileBuf[i]);
            }
        }
    }

    int xyCount = mReads;
    if(mFlags & BIT_PE_INTERLEAVED)
        xyCount = mReads/2;
    if(mHeader->hasX()) {
        writeLittleEndian(ofs, mXBufSize);
        ofs.write((const char*)mXBuf, mXBufSize);
    }
    if(mHeader->hasY()) {
        writeLittleEndian(ofs, mYBufSize);
        ofs.write((const char*)mYBuf, mYBufSize);
    }

    ofs.write(mName1Buf, mName1BufSize);
    if(mHeader->hasName2())
        ofs.write(mName2Buf, mName2BufSize);
    ofs.write(mStrandBuf, mStrandBufSize);
    ofs.write(mSeqBuf, mSeqBufSize);
    ofs.write((char*)mQualBuf, mQualBufSize);

    if( (mFlags & BIT_PE_INTERLEAVED) && (mHeader->mFlags & BIT_ENCODE_PE_BY_OVERLAP)) {
        ofs.write(mOverlapBuf, mReads/2);
    }

    if(mHeader->encodeNPos()) {
        ofs.write((char*)mNPosBuf, mNPosBufSize);
    }
}