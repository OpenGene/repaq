#include <stdint.h>
#include "endian.h"

bool isLittleEndian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 4; 
}

uint16 swapEndianess(uint16 num) {
    return ((num>>8) | (num<<8));
}

uint32 swapEndianess(uint32 num) {
    return ((num>>24)&0xff) | ((num<<8)&0xff0000) |  ((num>>8)&0xff00) | ((num<<24)&0xff000000); 
}

uint16 adaptToLittleEndian(uint16 num) {
    if(isLittleEndian())
        return num;
    else
        return swapEndianess(num);
}

uint32 adaptToLittleEndian(uint32 num) {
    if(isLittleEndian())
        return num;
    else
        return swapEndianess(num);
}

uint16 adaptToBigEndian(uint16 num) {
    if(!isLittleEndian())
        return num;
    else
        return swapEndianess(num);
}

uint32 adaptToBigEndian(uint32 num) {
    if(!isLittleEndian())
        return num;
    else
        return swapEndianess(num);
}

void writeLittleEndian(ostream& ofs, uint16 num) {
    uint16 data = adaptToLittleEndian(num);
    ofs.write((char*)&data, sizeof(uint16));
}

void writeLittleEndian(ostream& ofs, uint32 num) {
    uint32 data = adaptToLittleEndian(num);
    ofs.write((char*)&data, sizeof(uint32));
}

void writeBigEndian(ostream& ofs, uint16 num) {
    uint16 data = adaptToBigEndian(num);
    ofs.write((char*)&data, sizeof(uint16));
}

void writeBigEndian(ostream& ofs, uint32 num) {
    uint32 data = adaptToBigEndian(num);
    ofs.write((char*)&data, sizeof(uint32));
}

uint16 readLittleEndian16(istream& ifs) {
    uint16 data=0;
    ifs.read((char*)&data, sizeof(uint16));
    return adaptToLittleEndian(data);
}

uint32 readLittleEndian32(istream& ifs) {
    uint32 data=0;
    ifs.read((char*)&data, sizeof(uint32));
    return adaptToLittleEndian(data);
}

uint16 readBigEndian16(istream& ifs) {
    uint16 data=0;
    ifs.read((char*)&data, sizeof(uint16));
    return adaptToBigEndian(data);
}

uint32 readBigEndian32(istream& ifs) {
    uint32 data=0;
    ifs.read((char*)&data, sizeof(uint32));
    return adaptToBigEndian(data);
}
