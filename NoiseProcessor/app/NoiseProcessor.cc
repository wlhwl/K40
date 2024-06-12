#include "../include/Event_processor.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if(argc!=3) {
        std::cout<<"Usage: "<<argv[0]<<" <config_file> <root.file>"<<std::endl;
        return 1;
    }

    Event_processor ep(argv[1]);
    ep.set_up_rootfile(argv[2]);
    ep.Loop_events();
    return 0;
}
