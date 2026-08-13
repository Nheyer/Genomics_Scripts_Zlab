#include <iostream>
#include <fstream>
#include <string>
#include "argparse/argparse.hpp"

// Global vars


int parse_CLI(argparse::ArgumentParser  * Parser) {
    Parser->add_argument("--input","-i")
       .required()
       .help("This is your input Multiply sequence aliened fasta file");
    Parser->add_argument("--output","-o")
        .default_value(std::string("stdout"))
        .help("fasta file to output our consensus to");
    Parser->add_argument("--disagree-mask","-dm")
        .default_value('N')
        .help("What character Should be used as the mask");
    Parser->add_argument("--use-ambiguity", "-a")
        .default_value(false)
        .implicit_value(true)
        .help("Use IUPAC ambiguity codes");
    Parser->add_argument("--force-fna","-fna")
        .help("This is a fna file, don't use file ending")
        .default_value(false)
        .implicit_value(true);
     return 0;
}

int main(int argc, char *argv[]) {
    //Magic to get file redirection working well
    std::ofstream file;
    std::ostream *fo;
    //Parser
    auto INPUTS = argparse::ArgumentParser();
    parse_CLI( &INPUTS );
    INPUTS.parse_args(argc,argv);
    if(INPUTS.get("--output") != "stdout"){
        file.open(INPUTS.get("--output"));
        fo = &file;
    }else{
        fo =  &std::cout ;
    }



    if(INPUTS.get<bool>("-fna") or INPUTS.get("-i").substr(INPUTS.get("-i").find_last_of(".") + 1) == "fna"){
        std::cerr << "Opening: "
        << INPUTS.get("--input")
        << " as multiply aligned Nucleotide Sequence input!"
        << std::endl;
    }

     *fo << "Hello, World!" << std::endl;
    return 0;
}
