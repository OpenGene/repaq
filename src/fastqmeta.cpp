#include "fastqmeta.h"
#include <iostream>

FastqMeta::FastqMeta(){
    lane = 0;
    tile = 0;
    x = 0;
    y = 0;
    hasLaneTileXY = false;
}

void FastqMeta::print() {
    if(hasLaneTileXY)
        cerr << namePart1 << ":" << (int)lane << ":" << tile << ":" << x << ":" << y << namePart2 << endl;
    else
        cerr << namePart1;
}

//illumina sequence name line format
//@<instrument>:<run number>:<flowcell ID>:<lane>:<tile_no>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>

FastqMeta FastqMeta::parse(string str) {
    int colon = 0;
    int lastColonPos = 0;
    int coordsStartAt = 0;
    int coordsEndAt = 0;
    uint8 lane = 0;
    uint16 tile = 0;
    uint32 x = 0;
    uint32 y = 0;

    for(int i=0; i<str.length(); i++) {
        char c = str[i];
        if(c == ':') {
            colon++;
        }
        if(c == ':' || c==' ') {
            if(colon >= 4 && colon<=7) {
                string item = str.substr(lastColonPos+1, i-lastColonPos-1);
                int val = atoi(item.c_str());
                switch(colon) {
                    case 4: 
                        lane = val;
                        coordsStartAt = lastColonPos + 1;
                        break;
                    case 5: tile = val; break;
                    case 6: 
                        if(c == ':')
                            x = val;
                        break;
                    case 7: y = val; break;
                    default: break;
                }
                if(c==' ' && colon == 6)
                    y = val;
            }
        }
        if(c == ':') {
            lastColonPos = i;
        }
        if(c == ' ' || (c == ':' && colon==7)) {
            coordsEndAt = i;
            break;
        }
    }

    FastqMeta ret;
    if(coordsStartAt>0 && coordsEndAt>0) {
        ret.lane = lane;
        ret.tile = tile;
        ret.x = x;
        ret.y = y;
        ret.hasLaneTileXY = true;
        ret.namePart1 = str.substr(0, coordsStartAt - 1);
        ret.namePart2 = str.substr(coordsEndAt, str.length() - coordsEndAt);
    } else {
        ret.namePart1 = str;
    }
    return ret;
}

bool FastqMeta::test() {
    FastqMeta meta = FastqMeta::parse("@A00251:28:H3YV7DSXX:40:1101:2356:1000 1:N:0:TAAGTGGC");
    meta.print();
    if(meta.namePart1 != "@A00251:28:H3YV7DSXX") {
        cerr << "namePart1 is wrong:" << meta.namePart1 << endl;
        return false;
    }
    if(meta.lane != 40) {
        cerr << "lane is wrong:" << meta.lane << endl;
        return false;
    }
    if(meta.tile != 1101) {
        cerr << "tile is wrong:" << meta.tile << endl;
        return false;
    }
    if(meta.x != 2356) {
        cerr << "x is wrong:" << meta.x << endl;
        return false;
    }
    if(meta.y != 1000) {
        cerr << "y  is wrong:" << meta.y << endl;
        return false;
    }
    if(meta.namePart2 != " 1:N:0:TAAGTGGC") {
        cerr << "namePart2 is wrong:" << meta.namePart2 << endl;
        return false;
    }
    return true;
}