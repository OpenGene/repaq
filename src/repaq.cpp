#include "repaq.h"
#include "fastqreader.h"
#include "util.h"
#include "writer.h"
#include <stdio.h>

Repaq::Repaq(Options* opt){
    mOptions = opt;
}

void Repaq::run() {
    if(mOptions->compressMode) {
        if(mOptions->in2.empty())
            compress();
        else
            compressPE();
    }
    else {
        if(mOptions->out2.empty())
            decompress();
        else
            decompressPE();
    }
}

void Repaq::decompress(){
    RfqCodec codec;

    ifstream input;
    input.open(mOptions->in1, ios::in | ios::binary);

    Writer writer(mOptions->out1);

    RfqHeader* header = new RfqHeader();
    header->read(input);

    if(header->mFlags & BIT_PAIRED_END) {
        error_exit("The input RFQ file was encoded by paired-end FASTQ, you should specify <out1> and <out2>");
    }

    codec.setHeader(header);

    string outstr;
    while(!input.eof()) {
        RfqChunk* chunk = new RfqChunk(header);
        chunk->read(input);
        vector<Read*> reads = codec.decodeChunk(chunk);
        delete chunk;
        if(reads.size() == 0)
            break;

        if(!outstr.empty()) {
            // if it is the last chunk, we don't write the last line break
            if(input.eof() && !header->hasLineBreakAtEnd()) {
                string nobreak = outstr.substr(0, outstr.length()-1);
                writer.writeString(nobreak);
            }
            else
                writer.writeString(outstr);
            outstr.clear();
        }

        for(int r=0; r<reads.size(); r++) {
            outstr += reads[r]->toString();
            delete reads[r];
        }
    }

    if(!outstr.empty()) {
        if(!header->hasLineBreakAtEnd()) {
            string nobreak = outstr.substr(0, outstr.length()-1);
            writer.writeString(nobreak);
        }
        else
            writer.writeString(outstr);
    }
}

void Repaq::decompressPE(){
    RfqCodec codec;

    ifstream input;
    input.open(mOptions->in1, ios::in | ios::binary);

    Writer writer1(mOptions->out1);
    Writer writer2(mOptions->out2);

    RfqHeader* header = new RfqHeader();
    header->read(input);

    codec.setHeader(header);

    string outstr1;
    string outstr2;
    while(!input.eof()) {
        RfqChunk* chunk = new RfqChunk(header);
        chunk->read(input);
        vector<Read*> reads = codec.decodeChunk(chunk);
        delete chunk;
        if(reads.size() == 0)
            break;

        if(!outstr1.empty()) {
            // if it is the last chunk, we don't write the last line break
            if(input.eof() && !header->hasLineBreakAtEnd()) {
                string nobreak = outstr1.substr(0, outstr1.length()-1);
                writer1.writeString(nobreak);
            }
            else
                writer1.writeString(outstr1);
            outstr1.clear();

            // if it is the last chunk, we don't write the last line break
            if(input.eof() && !header->hasLineBreakAtEndR2()) {
                string nobreak = outstr2.substr(0, outstr2.length()-1);
                writer2.writeString(nobreak);
            }
            else
                writer2.writeString(outstr2);
            outstr2.clear();
        }

        for(int r=0; r<reads.size(); r++) {
            if(r%2==0)
                outstr1 += reads[r]->toString();
            else
                outstr2 += reads[r]->toString();
            delete reads[r];
        }
    }

    if(!outstr1.empty()) {
        if(!header->hasLineBreakAtEnd()) {
            string nobreak = outstr1.substr(0, outstr1.length()-1);
            writer1.writeString(nobreak);
        }
        else
            writer1.writeString(outstr1);
    }

    if(!outstr2.empty()) {
        if(!header->hasLineBreakAtEndR2()) {
            string nobreak = outstr2.substr(0, outstr2.length()-1);
            writer2.writeString(nobreak);
        }
        else
            writer2.writeString(outstr2);
    }
}

bool Repaq::hasLineBreakAtEnd(string& filename) {
    FILE* fp = fopen(filename.c_str(), "r");
    if(!fp)
        return false;
    fseek(fp, -1, SEEK_END);
    char c;
    fread(&c, 1, 1, fp);
    fclose(fp);
    return c=='\n';
}

void Repaq::compress(){
    RfqCodec codec;
    FastqReader reader(mOptions->in1);

    bool hasN = hasLineBreakAtEnd(mOptions->in1);

    ofstream out;
    out.open(mOptions->out1, ios::out | ios::binary);

    vector<Read*> reads;
    RfqHeader* header = NULL;
    uint32 totalBses = 0;
    while(true){
        Read* read = reader.read();
        if(!read){
            break;
        }
        reads.push_back(read);
        totalBses += read->length();
        if(totalBses >= mOptions->chunkSize) {
            if(header == NULL) {
                header = codec.makeHeader(reads);
                if(hasN)
                    header->mFlags |= BIT_HAS_LINE_BREAK_AT_END;
                header->write(out);
            }
            if(header == NULL)
                error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
            RfqChunk* chunk = codec.encodeChunk(reads);
            if(chunk) {
                chunk->write(out);
                delete chunk;
            }
            for(int r=0; r<reads.size(); r++)
                delete reads[r];
            reads.clear();
            totalBses = 0;
        }
    }
    if(reads.size() > 0) {
        if(header == NULL) {
            header = codec.makeHeader(reads);
            if(hasN)
                header->mFlags |= BIT_HAS_LINE_BREAK_AT_END;
            header->write(out);
        }
        if(header == NULL)
            error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
        RfqChunk* chunk = codec.encodeChunk(reads);
        if(chunk) {
            chunk->write(out);
            delete chunk;
        }
        for(int r=0; r<reads.size(); r++)
            delete reads[r];
        reads.clear();
    }
    out.flush();
    out.close();

    if(header) {
        delete header;
        header = NULL;
    }
}

void Repaq::compressPE(){
    RfqCodec codec;
    FastqReaderPair reader(mOptions->in1, mOptions->in2);

    bool hasN = hasLineBreakAtEnd(mOptions->in1);
    bool hasN2 = hasLineBreakAtEnd(mOptions->in2);

    ofstream out;
    out.open(mOptions->out1, ios::out | ios::binary);

    vector<ReadPair*> reads;
    RfqHeader* header = NULL;
    uint32 totalBses = 0;
    while(true){
        ReadPair* read = reader.read();
        if(!read){
            break;
        }
        reads.push_back(read);
        totalBses += read->mLeft->length() + read->mRight->length();
        if(totalBses >= mOptions->chunkSize) {
            if(header == NULL) {
                header = codec.makeHeader(reads);
                if(hasN)
                    header->mFlags |= BIT_HAS_LINE_BREAK_AT_END;
                if(hasN2)
                    header->mFlags |= BIT_HAS_LINE_BREAK_AT_END_R2;
                header->write(out);
            }
            if(header == NULL)
                error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
            RfqChunk* chunk = codec.encodeChunk(reads);
            if(chunk) {
                chunk->write(out);
                delete chunk;
            }
            for(int r=0; r<reads.size(); r++)
                delete reads[r];
            reads.clear();
            totalBses = 0;
        }
    }
    if(reads.size() > 0) {
        if(header == NULL) {
            header = codec.makeHeader(reads);
            if(hasN)
                header->mFlags |= BIT_HAS_LINE_BREAK_AT_END;
            if(hasN2)
                header->mFlags |= BIT_HAS_LINE_BREAK_AT_END_R2;
            header->write(out);
        }
        if(header == NULL)
            error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
        RfqChunk* chunk = codec.encodeChunk(reads);
        if(chunk) {
            chunk->write(out);
            delete chunk;
        }
        for(int r=0; r<reads.size(); r++)
            delete reads[r];
        reads.clear();
    }
    out.flush();
    out.close();

    if(header) {
        delete header;
        header = NULL;
    }
}