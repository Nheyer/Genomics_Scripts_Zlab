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
    Parser->add_argument("--verbose","-v")
        .help("Be loud?").flag()
     return 0;
}

int accumulate_seqs(argparse::ArgumentParser  * settings,std::vector<std::string> * Store){
    std::ifstream fileI;
    std::string ln;
    std::string inPath = settings->get("--inputs");
    int i = 0;
    if(settings->get<bool>("-v")){
        std::cerr << "we gor here wtf\n";
    }
    if(settings->get<bool>("-fna") or inPath.substr(inPath.find_last_of(".") + 1) == "fna"){
        std::cerr << "Opening: "
                  << settings->get("--i")
                  << " as multiply aligned Nucleotide Sequence input!"
                  << std::endl;
        fileI.open(settings->get("-i"));
    }
    while(getline(fileI, ln)){
        if(settings->get<bool>("-v")){
            std::cerr <<  "LN-"<< i << " read line :" << ln;
        }


        i++;
    }
    return 0;
}

int main(int argc, char *argv[]) {

    // decs
    std::ofstream fileO;
    std::ostream *fo;
    std::vector<std::string> Seqs;
    auto INPUTS = argparse::ArgumentParser();
    //Parser
    parse_CLI( &INPUTS );
    try {
        INPUTS.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << INPUTS;
        std::exit(1);
    }

    if(INPUTS.get<bool>("-v")) {
        std::cerr << "Parsed arguments" << std::endl;
    }
    //Magic to get file redirection working well
    if(INPUTS.get("--output") != "stdout"){
        fileO.open(INPUTS.get("--output"));
        fo = &fileO;
    }else{
        fo =  &std::cout ;
    }
    if(INPUTS.get<bool>("-v")) {
        std::cerr << INPUTS.get("--output")<<" Opened for output" << std::endl;
    }
    // actual logic
    accumulate_seqs(&INPUTS,&Seqs);


     *fo << "Hello, World!" << std::endl;
    return 0;
}
