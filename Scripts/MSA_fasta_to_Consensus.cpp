#include <iostream>
//#include <bits/stdc++.h>
#include <fstream>
#include <string>
#include <argparse/argparse.hpp>
#define DEBUG 6
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
            std::cerr << "File has mismatched sequence lengths, is this a Alligned file?" << std::endl;
            return -1;
        }
    }
    for (int pos=0;pos<len_str;pos++) {
        char C = Seqs[0].seq[pos];
        for (int seq_num = 1; seq_num < Seqs.size(); ++seq_num) {
#if DEBUG > 5
            std::cerr << "Sequence no. " << seq_num
            << " at position " << pos<< " is \""
            << Seqs[seq_num].seq[pos] << "\""
            << " Expected \"" << C << "\""
            << std::endl;
#endif
            if (Strict) {
                // if we are strict any non match goes to n mask
                if (Seqs[seq_num].seq[pos] != C) {
                    C = 'n';
                }
            } else {
                //todo put IUPAC ambiguaty here
            }

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
    if (int cons_rt = make_consensus(Seqs,&Consensus_FASTA); cons_rt  != 0) {return -2;}
    if(int write_rt = write_fasta(INPUTS.get("-o"), Consensus_FASTA); INPUTS.get<bool>("-v") && write_rt==0) {
        std::cerr <<" Consensus FASTA writen to :\t" << INPUTS.get("--output")<< std::endl;
    } else if (!write_rt) {return -3;}
    return 0;
}
