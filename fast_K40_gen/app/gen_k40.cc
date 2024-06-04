#include "../include/K40_generation.hpp"

int main(){
    for(int i=0; i<24220; i++){
        K40_generation k40("../config/config.yaml");
        k40.generate_K40(i);
    }
    
    return 0;
}