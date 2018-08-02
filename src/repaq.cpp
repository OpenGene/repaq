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
        if(mOptions->in2.empty() && !mOptions->interleavedInput)
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

    /*if(header->mFlags & BIT_PAIRED_END) {
        error_exit("The input RFQ file was encoded by paired-end FASTQ, you should specify <out1> and <out2>");
    }*/

    codec.setHeader(header);

    bool hasNoLineBreakAtEnd = true;

    RfqChunk* chunk = NULL;
    while(!input.eof() || chunk != NULL) {
        if(!chunk) {
            chunk = new RfqChunk(header);
            chunk->read(input);
        }
        vector<Read*> reads = codec.decodeChunk(chunk);
        if(reads.size() == 0) {
            delete chunk;
            chunk = NULL;
            break;
        }

        string outstr;

        for(int r=0; r<reads.size(); r++) {
            outstr += reads[r]->toString();
            delete reads[r];
        }

        bool isLastOne = false;
        bool hasNoLineBreakAtEnd = chunk->mFlags & BIT_HAS_NO_LINE_BREAK_AT_END;

        if(hasNoLineBreakAtEnd) {
            delete chunk;
            chunk = NULL;
            isLastOne = input.eof();
            // eof sometime doesn't work, we need to check the last chunk to see if it's valid
            if(!isLastOne) {
                chunk = new RfqChunk(header);
                chunk->read(input);
                isLastOne = chunk->mReads == 0;
            }
        }

        // if it's last chunk, we should check the last line break
        if(hasNoLineBreakAtEnd ) {
            if(isLastOne) {
                string nobreak = outstr.substr(0, outstr.length()-1);
                writer.writeString(nobreak);
                break;
            } else {
                writer.writeString(outstr);
                continue;
            }
        }
        else
            writer.writeString(outstr);

        delete chunk;
        chunk = NULL;
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

    if( (header->mFlags & BIT_PAIRED_END) == false) {
        error_exit("The input RFQ file was encoded by single-end FASTQ, you should not specify <out2>");
    }

    codec.setHeader(header);

    string outstr1;
    string outstr2;
    while(!input.eof()) {
        RfqChunk* chunk = new RfqChunk(header);
        chunk->read(input);
        vector<Read*> reads = codec.decodeChunk(chunk);
        if(reads.size() == 0) {
            delete chunk;
            break;
        }

        string outstr1;
        string outstr2;

        for(int r=0; r<reads.size(); r++) {
            if(r%2==0)
                outstr1 += reads[r]->toString();
            else
                outstr2 += reads[r]->toString();
            delete reads[r];
        }

        bool isLastOne = false;
        bool hasNoLineBreakAtEndR1 = chunk->mFlags & BIT_HAS_NO_LINE_BREAK_AT_END;
        bool hasNoLineBreakAtEndR2 = chunk->mFlags & BIT_HAS_NO_LINE_BREAK_AT_END_R2;

        if(hasNoLineBreakAtEndR1 || hasNoLineBreakAtEndR2) {
            delete chunk;
            chunk = NULL;
            isLastOne = input.eof();
            // eof sometime doesn't work, we need to check the last chunk to see if it's valid
            if(!isLastOne) {
                chunk = new RfqChunk(header);
                chunk->read(input);
                isLastOne = chunk->mReads == 0;
            }
        }

        if(hasNoLineBreakAtEndR1) {
            if(isLastOne) {
                string nobreak = outstr1.substr(0, outstr1.length()-1);
                writer1.writeString(nobreak);
            } else {
                writer1.writeString(outstr1);
                continue;
            }
        }
        else
            writer1.writeString(outstr1);

        if(hasNoLineBreakAtEndR2) {
            if(isLastOne) {
                string nobreak = outstr2.substr(0, outstr2.length()-1);
                writer2.writeString(nobreak);
            } else {
                writer2.writeString(outstr2);
                continue;
            }
        }
        else
            writer2.writeString(outstr2);

        delete chunk;
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
                header->write(out);
            }
            if(header == NULL)
                error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
            RfqChunk* chunk = codec.encodeChunk(reads);
            if(chunk) {
                if(reader.hasNoLineBreakAtEnd())
                    chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;
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
            header->write(out);
        }
        if(header == NULL)
            error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
        RfqChunk* chunk = codec.encodeChunk(reads);
        if(chunk) {
            if(reader.hasNoLineBreakAtEnd())
                chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;
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
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, false, mOptions->interleavedInput);

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
                header->write(out);
            }
            if(header == NULL)
                error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
            RfqChunk* chunk = codec.encodeChunk(reads);
            if(chunk) {
                if(reader.mLeft->hasNoLineBreakAtEnd())
                    chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;
                if(reader.mRight->hasNoLineBreakAtEnd())
                    chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END_R2;
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
            header->write(out);
        }
        if(header == NULL)
            error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
        RfqChunk* chunk = codec.encodeChunk(reads);
        if(chunk) {
            if(reader.mLeft->hasNoLineBreakAtEnd())
                chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;
            if(reader.mRight->hasNoLineBreakAtEnd())
                chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END_R2;
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