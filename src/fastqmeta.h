#ifndef FASTQ_META_H
#define FASTQ_META_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "common.h"

using namespace std;

/*
* the query name of FASTQ is usually like
* @A00250:26:H3YTWDSXX:1:1101:16776:17989 1:N:0:ACTGTTCC
* can be divided into 
* <QNAME> = <NAME1>:<LANE>:<TILE>:<X>:<Y><NAME2>
* <NAME1> = @A00250:26:H3YTWDSXX:
* <NAME2> =  1:N:0:ACTGTTCC
*/

class FastqMeta{
public:
    FastqMeta();
    void print();
    static FastqMeta parse(string str);
    static bool test();

public:
    string namePart1;
    string namePart2;
    uint8 lane;
    uint16 tile;
    uint32 x;
    uint32 y;
    bool hasLaneTileXY;
};

#endif