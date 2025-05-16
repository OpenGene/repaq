#include "repaq.h"
#include "fastqreader.h"
#include "util.h"
#include "writer.h"
#include <stdio.h>
#include <sstream>

Repaq::Repaq(Options* opt){
    mOptions = opt;
}

void Repaq::run() {
    if(mOptions->mode == REPAQ_COMPRESS) {
        if(mOptions->in2.empty() && !mOptions->interleavedInput)
            compress();
        else
            compressPE();
    }
    else if(mOptions->mode == REPAQ_DECOMPRESS) {
        if(mOptions->out2.empty())
            decompress();
        else
            decompressPE();
    }
    else if(mOptions->mode == REPAQ_COMPARE) {
        if(mOptions->in2.empty())
            compare();
        else
            comparePE();
    }
    else {
        error_exit("no mode specified, you should specify one of compress/decompress/compare mode");
    }
}

void Repaq::compare(){
    RfqCodec codec;

    ifstream input;
    input.open(mOptions->rfqCompare, ios::in | ios::binary);

    FastqReader reader(mOptions->in1);

    RfqHeader* header = new RfqHeader();
    header->read(input);

    codec.setHeader(header);

    long fqReads = 0;
    long fqBases = 0;
    long rfqReads = 0;
    long rfqBases = 0;

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

        for(int r=0; r<reads.size(); r++) {
            Read* rfq = reads[r];
            rfqBases += rfq->length();
            rfqReads++;

            Read* fq = reader.read();
            if(!fq) {
                string msg = "The RFQ file has more reads than the FASTQ file.";
                msg += " The RFQ file has >= " + to_string(rfqReads);
                msg += " reads, while the FASTQ file only has " + to_string(fqReads) + " reads";
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }

            fqReads++;
            fqBases+=fq->length();

            //check read
            if(rfq->mName != fq->mName) {
                string msg = "The RFQ file and FASTQ file have different name in the " + to_string(rfqReads) + " read. ";
                msg += rfq->mName + " | " + fq->mName;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }
            else if(rfq->mSeq.mStr != fq->mSeq.mStr) {
                string msg = "The RFQ file and FASTQ file have different sequence in the " + to_string(rfqReads) + " read. ";
                msg += rfq->mSeq.mStr + " | " + fq->mSeq.mStr;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }
            else if(rfq->mStrand != fq->mStrand) {
                string msg = "The RFQ file and FASTQ file have different strand in the " + to_string(rfqReads) + " read. ";
                msg += rfq->mStrand + " | " + fq->mStrand;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }
            else if(rfq->mQuality != fq->mQuality) {
                string msg = "The RFQ file and FASTQ file have different quality in the " + to_string(rfqReads) + " read. ";
                msg += rfq->mQuality + " | " + fq->mQuality;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }

            delete rfq;
            delete fq;
        }

        delete chunk;
        chunk = NULL;
    }

    if(reader.read()) {
        fqReads++;
        string msg = "The FASTQ file has more reads than the RFQ file.";
        msg += " The FASTQ file has >= " + to_string(fqReads);
        msg += " reads, while the RFQ file only has " + to_string(rfqReads) + " reads";
        reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
        return;
    }

    reportCompareResult(true, "", fqReads, fqBases, rfqReads, rfqBases);
}

void Repaq::comparePE(){
    RfqCodec codec;

    ifstream input;
    input.open(mOptions->rfqCompare, ios::in | ios::binary);

    FastqReaderPair reader(mOptions->in1, mOptions->in2);

    RfqHeader* header = new RfqHeader();
    header->read(input);

    codec.setHeader(header);

    long fqReads = 0;
    long fqBases = 0;
    long rfqReads = 0;
    long rfqBases = 0;

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

        ReadPair* pair=NULL;
        for(int r=0; r<reads.size(); r++) {
            Read* rfq = reads[r];
            rfqBases += rfq->length();
            rfqReads++;

            if(pair == NULL)
                pair = reader.read();
            if(!pair) {
                string msg = "The RFQ file has more reads than the FASTQ file.";
                msg += " The RFQ file has >= " + to_string(rfqReads/2);
                msg += " pairs, while the FASTQ file only has " + to_string(fqReads/2) + " pairs";
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }

            Read* fq;
            if(rfqReads%2==1)
                fq = pair->mLeft;
            else
                fq = pair->mRight;

            fqReads++;
            fqBases+=fq->length();

            //check read
            if(rfq->mName != fq->mName) {
                string msg = "The RFQ file and FASTQ file have different name in the " + to_string(rfqReads/2) + " pair. ";
                msg += rfq->mName + " | " + fq->mName;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }
            else if(rfq->mSeq.mStr != fq->mSeq.mStr) {
                string msg = "The RFQ file and FASTQ file have different sequence in the " + to_string(rfqReads/2) + " pair. ";
                msg += rfq->mSeq.mStr + " | " + fq->mSeq.mStr;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }
            else if(rfq->mStrand != fq->mStrand) {
                string msg = "The RFQ file and FASTQ file have different strand in the " + to_string(rfqReads/2) + " pair. ";
                msg += rfq->mStrand + " | " + fq->mStrand;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }
            else if(rfq->mQuality != fq->mQuality) {
                string msg = "The RFQ file and FASTQ file have different quality in the " + to_string(rfqReads/2) + " pair. ";
                msg += rfq->mQuality + " | " + fq->mQuality;
                reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
                return;
            }

            delete rfq;
            if(rfqReads%2==0) {
                delete pair;
                pair = NULL;
            }
        }

        delete chunk;
        chunk = NULL;
    }

    if(reader.read()) {
        fqReads++;
        string msg = "The FASTQ file has more reads than the RFQ file.";
        msg += " The FASTQ file has >= " + to_string(fqReads/2);
        msg += " pairs, while the RFQ file only has " + to_string(rfqReads/2) + " pairs";
        reportCompareResult(false, msg, fqReads, fqBases, rfqReads, rfqBases);
        return;
    }

    reportCompareResult(true, "", fqReads, fqBases, rfqReads, rfqBases);
}

void Repaq::reportCompareResult(bool passed, string msg, long fqReads, long fqBases, long rfqReads, long rfqBases) {
    string json = "{\n";
    if(passed)
        json += "\t\"result\":\"passed\",\n";
    else
        json += "\t\"result\":\"failed\",\n";

    json += "\t\"msg\":\"" + msg + "\",\n";

    json += "\t\"fastq_reads\":" + to_string(fqReads) + ",\n";
    json += "\t\"rfq_reads\":" + to_string(rfqReads) + ",\n";
    json += "\t\"fastq_bases\":" + to_string(fqBases) + ",\n";
    json += "\t\"rfq_bases\":" + to_string(rfqBases) + "\n";

    json += "}\n";

    if(!mOptions->jsonFileForCompare.empty()) {
        ofstream out(mOptions->jsonFileForCompare);
        out.write(json.c_str(), json.length());
        out.flush();
        out.close();
    }

    cout << json;
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

bool Repaq::doubleCheckAndOutput(RfqChunk* chunk, RfqCodec& codec4check, RfqHeader* header4check, vector<Read*>& reads, ostream& out) {
    ostringstream oss;
    istringstream iss;
    chunk->write(oss);
    out<<oss.str();
    iss.str(oss.str());
    RfqChunk* chunk4check = new RfqChunk(header4check);
    chunk4check->read(iss);
    vector<Read*> reads4check = codec4check.decodeChunk(chunk4check);
    //todo: check reads identity
    if(reads4check.size() != reads.size()) {
        error_exit("encoding error in chunk, the output will be wrong, quit now!");
        return false;
    }
    bool identical = true;
    for(int i=0; i<reads4check.size(); i++) {
        Read* rfq = reads4check[i];
        Read* fq = reads[i];
        if(rfq->mName != fq->mName) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mName << endl << "got:\n" << rfq->mName << endl;
            identical = false;
        }
        else if(rfq->mSeq.mStr != fq->mSeq.mStr) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mSeq.mStr << endl << "got:\n" << rfq->mSeq.mStr << endl;
            identical = false;
        }
        else if(rfq->mStrand != fq->mStrand) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mStrand << endl << "got:\n" << rfq->mStrand << endl;
            identical = false;
        }
        else if(rfq->mQuality != fq->mQuality) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mQuality << endl << "got:\n" << rfq->mQuality << endl;
            identical = false;
        }
    }
    delete chunk4check;
    for(int r=0; r<reads4check.size(); r++)
        delete reads4check[r];
    reads4check.clear();
    return identical;
}

bool Repaq::doubleCheckAndOutput(RfqChunk* chunk, RfqCodec& codec4check, RfqHeader* header4check, vector<ReadPair*>& pairs, ostream& out) {
    ostringstream oss;
    istringstream iss;
    chunk->write(oss);
    out<<oss.str();
    iss.str(oss.str());
    RfqChunk* chunk4check = new RfqChunk(header4check);
    chunk4check->read(iss);
    vector<Read*> reads4check = codec4check.decodeChunk(chunk4check);
    //todo: check reads identity
    if(reads4check.size() != pairs.size()*2) {
        error_exit("encoding error in chunk, the output will be wrong, quit now!");
        return false;
    }
    bool identical = true;
    for(int i=0; i<reads4check.size(); i++) {
        Read* rfq = reads4check[i];

        Read* fq;
        if(i%2==0) {
            fq = pairs[i/2]->mLeft;
        } else {
            fq = pairs[i/2]->mRight;
            if(header4check->supportInterleaved())
                rfq->changeToReverseComplement();
        }
        if(rfq->mName != fq->mName) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mName << endl << "got:\n" << rfq->mName << endl;
            identical = false;
        }
        if(rfq->mSeq.mStr != fq->mSeq.mStr) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mSeq.mStr << endl << "got:\n" << rfq->mSeq.mStr << endl;
            identical = false;
        }
        if(rfq->mStrand != fq->mStrand) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mStrand << endl << "got:\n" << rfq->mStrand << endl;
            identical = false;
        }
        if(rfq->mQuality != fq->mQuality) {
            cerr << "integrity check failure \nexpected: " << endl << fq->mQuality << endl << "got:\n" << rfq->mQuality << endl;
            identical = false;
        }
    }
    delete chunk4check;
    for(int r=0; r<reads4check.size(); r++)
        delete reads4check[r];
    reads4check.clear();
    return identical;
}

void Repaq::compress(){
    RfqCodec codec;
    FastqReader reader(mOptions->in1);

    ofstream out;
    out.open(mOptions->out1, ios::out | ios::binary);

    // for double check
    RfqCodec codec4check;
    ostringstream ossHeader;
    RfqHeader* header4check = new RfqHeader();

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
                header->write(ossHeader);
                out<<ossHeader.str();
                // for double check
                istringstream issHeader;
                issHeader.str(ossHeader.str());
                header4check->read(issHeader);
                codec4check.setHeader(header4check);
                if(!header->identicalWith(header4check)) {
                    error_exit("encoding error in header, the output will be wrong, quit now!");
                }
            }
            if(header == NULL)
                error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
            RfqChunk* chunk = codec.encodeChunk(reads);
            if(chunk) {
                if(reader.hasNoLineBreakAtEnd())
                    chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;

                // for double check
                if(mOptions->doubleCheck) {
                    doubleCheckAndOutput(chunk, codec4check, header4check, reads, out);
                } else {
                    chunk->write(out);
                }

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
            header->write(ossHeader);
            out<<ossHeader.str();
            // for double check
            istringstream issHeader;
            issHeader.str(ossHeader.str());
            header4check->read(issHeader);
            codec4check.setHeader(header4check);
            if(!header->identicalWith(header4check)) {
                error_exit("encoding error in header, the output will be wrong, quit now!");
            }
        }
        if(header == NULL)
            error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
        RfqChunk* chunk = codec.encodeChunk(reads);
        if(chunk) {
            if(reader.hasNoLineBreakAtEnd())
                chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;

            // for double check
            if(mOptions->doubleCheck) {
                doubleCheckAndOutput(chunk, codec4check, header4check, reads, out);
            } else {
                chunk->write(out);
            }

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

    // for double check
    if(header4check) {
        delete header4check;
        header4check = NULL;
    }
}

void Repaq::compressPE(){
    RfqCodec codec;
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, false, mOptions->interleavedInput);

    ofstream out;
    out.open(mOptions->out1, ios::out | ios::binary);

    // for double check
    RfqCodec codec4check;
    ostringstream ossHeader;
    RfqHeader* header4check = new RfqHeader();

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
                header->write(ossHeader);
                out<<ossHeader.str();
                // for double check
                istringstream issHeader;
                issHeader.str(ossHeader.str());
                header4check->read(issHeader);
                // copy the mSupportInterleaved flag which is not stored in file
                header4check->mSupportInterleaved = header->mSupportInterleaved;
                codec4check.setHeader(header4check);
                if(!header->identicalWith(header4check)) {
                    error_exit("encoding error in header, the output will be wrong, quit now!");
                }
            }
            if(header == NULL)
                error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
            RfqChunk* chunk = codec.encodeChunk(reads);
            if(chunk) {
                bool noLineBreakAtEnd = reader.mLeft->hasNoLineBreakAtEnd();
                bool noLineBreakAtEndR2;
                if(!mOptions->interleavedInput)
                    noLineBreakAtEndR2 = reader.mRight->hasNoLineBreakAtEnd();
                else
                    noLineBreakAtEndR2 = noLineBreakAtEnd;
                if(noLineBreakAtEnd)
                    chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;
                if(noLineBreakAtEndR2)
                    chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END_R2;
                
                // for double check
                if(mOptions->doubleCheck) {
                    doubleCheckAndOutput(chunk, codec4check, header4check, reads, out);
                } else {
                    chunk->write(out);
                }

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
            header->write(ossHeader);
            out<<ossHeader.str();
            // for double check
            istringstream issHeader;
            issHeader.str(ossHeader.str());
            header4check->read(issHeader);
            // copy the mSupportInterleaved flag which is not stored in file
            header4check->mSupportInterleaved = header->mSupportInterleaved;
            codec4check.setHeader(header4check);
            if(!header->identicalWith(header4check)) {
                error_exit("encoding error in header, the output will be wrong, quit now!");
            }
        }
        if(header == NULL)
            error_exit("failed to encode, please confirm the input FASTQ file is valid and not empty");
        RfqChunk* chunk = codec.encodeChunk(reads);
        if(chunk) {
            bool noLineBreakAtEnd = reader.mLeft->hasNoLineBreakAtEnd();
            bool noLineBreakAtEndR2;
            if(!mOptions->interleavedInput)
                noLineBreakAtEndR2 = reader.mRight->hasNoLineBreakAtEnd();
            else
                noLineBreakAtEndR2 = noLineBreakAtEnd;
            if(noLineBreakAtEnd)
                chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END;
            if(noLineBreakAtEndR2)
                chunk->mFlags |= BIT_HAS_NO_LINE_BREAK_AT_END_R2;

            // for double check
            if(mOptions->doubleCheck) {
                doubleCheckAndOutput(chunk, codec4check, header4check, reads, out);
            } else {
                chunk->write(out);
            }
        
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