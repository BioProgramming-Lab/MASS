#include "DataStructure.cpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

//int generateRLT(int )

int main (int argc, char* argv[]) { // argv[1] = directory, argv[2] = domain number of receptors, argv[3] = domain number of ligand 1, argv[4] = domain number of ligand 2, ...
    int n;
    int number_of_ligands = 0;
    std::string RLT_type = "simplified";
    std::string dir = "/usr/src/CodePile/RLT_files/";
    Set_of_Ligands LigandSuite;

    if (argc == 1) {    // if no arguments are given
        if (! std::filesystem::exists(dir)) {
            std::cout << "ERROR: /usr/src/CodePile/RLT_files/ directory does not exist.\n";
            exit(1);
        }

        std::cout << "input receptor-ligand topology. (e.g. \"4 2 1\")\n";
        
        std::cin >> n;
        dir += std::to_string(n);
        int a;    
        while(std::cin >> a) {
            Ligand* b = new Ligand(a);
            LigandSuite.suite.push_back(b);
            dir += "-"+std::to_string(a);
            number_of_ligands++;
            if (std::cin.get() == '\n')
                break;
        }
    }

    else {
        dir = argv[1];
        if (dir.back() != '/')
            dir += "/";
        if (! std::filesystem::exists(dir)) {   // check if directory exists
            std::cout << "ERROR: " << dir << " directory does not exist.\n";
            exit(1);    // terminate with error, add do you wonna create directory? --------------------------------------------------------------------------------------
        }

        RLT_type = argv[2];

        n = atoi(argv[3]);
        dir += std::to_string(n);
        for (int i = 4; i < argc; i++) {
            Ligand* b = new Ligand(atoi(argv[i]));
            LigandSuite.suite.push_back(b);
            dir += "-"+std::to_string(atoi(argv[i]));
            number_of_ligands++;
        }
    }

    std::ofstream basic(dir+"number_of_ligands.txt");    // number of ligands
    basic << number_of_ligands;

    //LigandSuite.suite.push_back(&c);
    relations crater(n, LigandSuite);
    crater.GenerateRelation(RLT_type);
    //crater.show_relations();
    crater.show_relations_in_dot(dir);    
    crater.show_relations_in_unframed_son_dot(dir);
    crater.show_states_in_file(dir);
    crater.show_relations_in_file(dir);
}