#include <iostream>
#include <cctype>
//#include <bits/stdc++.h>
#include <fstream>
#include <string>
#include <argparse/argparse.hpp>
#include <climits>
// DEBUG LVL 0-9
#define DEBUG 0 
// Nucliotide IUPAC defs
#define nuc_A 3
#define nuc_T 5
#define nuc_C 7
#define nuc_G 11
#define nuc_R (nuc_A * nuc_G)
#define nuc_Y (nuc_C * nuc_T)
#define nuc_K (nuc_G * nuc_T)
#define nuc_M (nuc_A * nuc_C)
#define nuc_S (nuc_C * nuc_G)
#define nuc_W (nuc_A * nuc_T)
#define nuc_B (nuc_C * nuc_G * nuc_T)
#define nuc_D (nuc_A * nuc_G * nuc_T)
#define nuc_H (nuc_A * nuc_C * nuc_T)
#define nuc_V (nuc_A * nuc_C * nuc_G)
#define nuc_N (nuc_A * nuc_C * nuc_T * nuc_G)
// Global vars


struct fasta_entry {
    std::string name;
    std::string description;
    std::string seq;
    void clear() {
        seq.clear();
        name.clear();
        description.clear();
    }
};

char clean_nucliotide(char c) {
    char nucliotide = toupper(c);
    if (nucliotide == 'U') {
        nucliotide = 'T';
    }
    return nucliotide;
}
int encode_nucliotide(char c) {
    if (c == 'A') {return nuc_A;}
    if (c == 'T') {return nuc_T;}
    if (c == 'C') {return nuc_C;}
    if (c == 'G') {return nuc_G;}
    if (c == 'R') {return nuc_R;}
    if (c == 'Y') {return nuc_Y;}
    if (c == 'K') {return nuc_K;}
    if (c == 'M') {return nuc_M;}
    if (c == 'S') {return nuc_S;}
    if (c == 'W') {return nuc_W;}
    if (c == 'B') {return nuc_B;}
    if (c == 'D') {return nuc_D;}
    if (c == 'H') {return nuc_H;}
    if (c == 'V') {return nuc_V;}
    if (c == 'N') {return nuc_N;}

}


char decode_nucliotide(const unsigned long long int &Code) {
    if (Code % nuc_N == 0 ){ return 'N';}
    if (Code % nuc_V == 0 ){ return 'V';}
    if (Code % nuc_H == 0 ){ return 'H';}
    if (Code % nuc_D == 0 ){ return 'D';}
    if (Code % nuc_B == 0 ){ return 'B';}
    if (Code % nuc_W == 0 ){ return 'W';}
    if (Code % nuc_S == 0 ){ return 'S';}
    if (Code % nuc_M == 0 ){ return 'M';}
    if (Code % nuc_K == 0 ){ return 'K';}
    if (Code % nuc_Y == 0 ){ return 'Y';}
    if (Code % nuc_R == 0 ){ return 'R';}
    if (Code % nuc_A == 0 ){ return 'A';}
    if (Code % nuc_T == 0 ){ return 'T';}
    if (Code % nuc_C == 0 ){ return 'C';}
    if (Code % nuc_G == 0 ){ return 'G';}
    std::cerr << "Failed to decode \"" << Code << "\"" << " as a multiple of "
    << nuc_A << " (A), "
    << nuc_T << " (T), "
    << nuc_C << " (C), and "
    << nuc_G << " (G)"
    << std::endl;
    exit(-6);
}

int parse_CLI(argparse::ArgumentParser *Parser, int &argument_number, const char *const* argument_string) {
    Parser->add_argument("--input","-i")
       .required()
       .help("This is your input Multiply sequence aliened fasta file");
    Parser->add_argument("--output","-o")
        .default_value("stdout")
        .help("fasta file to output our consensus to");
    Parser->add_argument("--disagree-mask","-dm")
        .default_value('N')
        .help("What character Should be used as the mask");
    Parser->add_argument("--use-ambiguity", "-a")
        .help("Use IUPAC ambiguity codes, takes a bit longer")
        .flag();
    Parser->add_argument("--force-fna","-fna")
        .help("This is a Nucliotide acid sequence file, don't use file ending")
        .flag();
    Parser->add_argument("--force-faa","-faa")
        .help("This is a Amino acid sequence file, don't use file ending")
        .flag();
    Parser->add_argument("--verbose","-v")
        .help("Be loud?")
        .flag();
    try {
        Parser->parse_args( argument_number, argument_string);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << *Parser ;
        std::exit(1);
    }
    if(Parser->get<bool>("-v")) {
        std::cerr << "Parsed arguments" << std::endl;
    }
     return 0;
}

int accumulate_seqs_from_fasta(std::string inPath,std::vector<fasta_entry> * Store,int is_nucliotides = -1){
    std::ifstream fileI;
    std::string ln;
    fasta_entry temp;
    int i=0 , j = 0;
    if(is_nucliotides == 0 or inPath.substr(inPath.find_last_of(".") + 1) == "faa"){
        std::cerr << "Opening: "
                  << inPath
                  << " as multiply aligned Protein Sequence input!"
                  << std::endl;
        is_nucliotides = true;
    }else if(is_nucliotides == 1 or inPath.substr(inPath.find_last_of(".") + 1) == "fna"){
            std::cerr << "Opening: "
                << inPath
                << " as multiply aligned Nucliotide Sequence input!"
                << std::endl;
        is_nucliotides = false;
    } else {
        std::cerr << "Found unrecognised file " << inPath << " please use the -fna or -faa flags or change your file endings" <<std::endl;
        return -1;
    }
    //open file
    fileI.open(inPath);
    //loop through file
    while(getline(fileI, ln)){
#if DEBUG > 0
            std::cerr <<  "LN-"<< i << " read line:\t" << ln << std::endl;
#endif

        if (ln[0] == '>') {
            //if if is a name line get the name, and put any thing after a space in decription
            if (j>0) {
                Store->push_back(temp);
                temp.clear();
            }
            auto cur_space = ln.find_first_of(' ');
            if (cur_space != std::string::npos) {
                // if there is a space anywhere we need to grab comments
                temp.name = ln.substr(1, cur_space);
                temp.description = ln.substr(cur_space + 1);
            }else{
                temp.name = ln.substr(1);
            }
            j++;
        } else {
            // we are in a sequence block just add it to the seq
            temp.seq.append(ln);
        }
        i++;
    }
    // add the last seq
    Store->push_back(temp);
    std::cerr << "Found: " << j << " sequense(s) in " << i << " lines" << std::endl;
    return 0;
}

int write_fasta(const std::string& pathO,std::vector<fasta_entry> &SEQ) {
    std::ofstream fileO;
    std::ostream *fo;
    if(pathO != "stdout"){
        fileO.open(pathO);
        fo = &fileO;
    }else{
        fo =  &std::cout ;
    }
    for(const fasta_entry  &temp : SEQ) {
        *fo << temp.name << " " << temp.description << std::endl;
        for (int j=0;j<temp.seq.size();j = j + 80) {
            *fo << temp.seq.substr(j,80) << std::endl;
        }
    }
    return 0;
}

int make_consensus(std::vector<fasta_entry>  Seqs,std::vector<fasta_entry> * Output,bool Strict = true) {
    fasta_entry consensus_entry;
    consensus_entry.name = ">Consensus_Sequence_of:";
    consensus_entry.description = "";
    unsigned long len_str = 0;
    if (Seqs.size()==1) {
        std::cerr << "Only found one seq, a consensus needs at least 2, is this the right file?" << std::endl;
        return -1;
    }
    for (int i=0;i<Seqs.size();i++) {
#if DEBUG > 2
        std::cerr << "Seq no "<< i << " is " << Seqs[i].name << "with size " << Seqs[i].seq.length()
#if DEBUG > 3
        << " and value: " << Seqs[i].seq
#endif
        << std::endl;
#endif
        if (i==0){
            len_str = Seqs[i].seq.length();
        }else if (len_str != Seqs[i].seq.length()) {
            std::cerr << "File has mismatched sequence lengths, is this a Aligned file?" << std::endl;
            return -1;
        }
    }
    for (int pos=0;pos<len_str;pos++) {
        char C = clean_nucliotide(Seqs[0].seq[pos]);
        // the loop below only checks the other sequences, so catch a gap in the
        // first one here, otherwise strict mode masks it to N instead
        if (C == '-') {
            consensus_entry.seq = consensus_entry.seq + '-';
            continue;
        }
        unsigned long long nuc_accumulative_encoding = encode_nucliotide(C);
        for (int seq_num = 1; seq_num < Seqs.size(); ++seq_num) {
#if DEBUG > 5
            std::cerr << "Sequence no. " << seq_num
            << " at position " << pos<< " is \""
            << Seqs[seq_num].seq[pos] << "\""
            << " Expected \"" << C << "\""
            << std::endl;
#endif
            if (clean_nucliotide(Seqs[seq_num].seq[pos]) == '-') {
                // if we hit a - we can just skip the rest of the sequences but we need to flag it
                C = '-';
                break;
            }
            /// no abiguity codes is much easier so just do it here
            if (Strict) {
                // if we are strict any non match goes to n mask
                if (clean_nucliotide(Seqs[seq_num].seq[pos]) != C) {
                    C = 'N';
                }
            // This is more difficult so we are multiplying primes here
            } else {
                if (nuc_accumulative_encoding < ULONG_LONG_MAX/nuc_N) {
                    nuc_accumulative_encoding = nuc_accumulative_encoding * encode_nucliotide(clean_nucliotide(Seqs[seq_num].seq[pos]));
                }else {
#if DEBUG > 1
                    std::cerr << "We had too many sequences, compressing" << std::endl;
#endif
                    int new_nuc_accumulator = 1;
                    if (nuc_accumulative_encoding % nuc_A == 0){new_nuc_accumulator = new_nuc_accumulator * nuc_A;}
                    if (nuc_accumulative_encoding % nuc_T == 0){new_nuc_accumulator = new_nuc_accumulator * nuc_T;}
                    if (nuc_accumulative_encoding % nuc_C == 0){new_nuc_accumulator = new_nuc_accumulator * nuc_C;}
                    if (nuc_accumulative_encoding % nuc_G == 0){new_nuc_accumulator = new_nuc_accumulator * nuc_G;}
                    nuc_accumulative_encoding = new_nuc_accumulator * encode_nucliotide(clean_nucliotide(Seqs[seq_num].seq[pos]));
                    }
            }
        }
        if (!Strict && C != '-') {
            C = decode_nucliotide(nuc_accumulative_encoding);
        }
        consensus_entry.seq = consensus_entry.seq + C;
    }
    Output->push_back(consensus_entry);
    return 0;
}

int main(int argc, char *argv[]) {
    // decs
    std::vector<fasta_entry> Seqs;
    auto INPUTS = argparse::ArgumentParser();
    std::vector<fasta_entry> Consensus_FASTA;
    //Parser
    parse_CLI( &INPUTS , argc, argv);

    // actual logic
    if (int acc_rt = accumulate_seqs_from_fasta(INPUTS.get("-i"),&Seqs); acc_rt  != 0) {return -1;}
    if (int cons_rt = make_consensus(Seqs,&Consensus_FASTA,!INPUTS.get<bool>("-a")); cons_rt  != 0) {return -2;}
    if(int write_rt = write_fasta(INPUTS.get("-o"), Consensus_FASTA); INPUTS.get<bool>("-v") && write_rt==0) {
        std::cerr <<" Consensus FASTA writen to :\t" << INPUTS.get("--output")<< std::endl;
    } else if (!write_rt) {return -3;}
    return 0;
}
