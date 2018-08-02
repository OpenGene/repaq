#ifndef REPAQ_ENDIAN_H
#define REPAQ_ENDIAN_H

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include <iostream>
#include <fstream>

using namespace std;

bool isLittleEndian();
uint16 swapEndianess(uint16 num);
uint32 swapEndianess(uint32 num);
uint16 adaptToLittleEndian(uint16 num);
uint32 adaptToLittleEndian(uint32 num);
uint16 adaptToBigEndian(uint16 num);
uint32 adaptToBigEndian(uint32 num);
void writeLittleEndian(ofstream& ofs, uint16 num);
void writeLittleEndian(ofstream& ofs, uint32 num);
void writeBigEndian(ofstream& ofs, uint16 num);
void writeBigEndian(ofstream& ofs, uint32 num);
uint16 readLittleEndian16(ifstream& ifs);
uint32 readLittleEndian32(ifstream& ifs);
uint16 readBigEndian16(ifstream& ifs);
uint32 readBigEndian32(ifstream& ifs);
#endif