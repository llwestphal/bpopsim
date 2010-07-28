#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "commandline.h"
#include <boost/program_options.hpp>

using namespace boost::program_options; 


void cline(int argc, char* argv[]) {
// setup and parse configuration options:
        options_description cmdline_options("Allowed options");
        cmdline_options.add_options()
        ("help,h", "produce this help message")
        ("bam,b", value<std::string>(), "bam file containing sequences to be aligned")
        ("fasta,f", value<std::string>(), "FASTA file of reference sequence")
        ("output,o", value<std::string>(), "output directory")
        ("readfile,r", value<std::vector<std::string> >(), "names of readfiles (no extension)");

        variables_map options;
        store(parse_command_line(argc, argv, cmdline_options), options);
        notify(options);

        // make sure that the config options are good:
        if(options.count("help")
                 || !options.count("bam")
                 || !options.count("fasta")
                 || !options.count("output")
                 || !options.count("readfile")) {
                std::cout << "Usage: error_count --bam <sequences.bam> --fasta <reference.fasta> --output <path> --readfile <filename>" << std::endl;
                std::cout << cmdline_options << std::endl;
        }
}

