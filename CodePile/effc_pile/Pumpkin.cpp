// data structure of preprocessed data is called pumpkin
// basically, this is a handwriting 2D array of PDF*

// this file is for preprocess most of the convolution works
// memory recycle system didn't set up in this file ----------------------------------------------------------------------------

#ifndef Pumpkin
#define Pumpkin

#include <vector>
#include "PDFdef.cpp"
#include "GraphStructure.cpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

class pumpkin {
    public:
    std::vector<PDF*> seeds;    // handwritting  2D array
    std::vector<PDF*> stalks;   // build pdfs of linkers

    pumpkin (std::istream &InPut) {
        std::string line;
        // read in num_of_BDomains
        std::getline(InPut, line);
        this->num_of_BDomains = std::stoi (line, nullptr, 10);
        // read in stalks
        for (int i = 0; i < this->num_of_BDomains-1; i++) {
            std::getline(InPut, line);
            int bingze = std::stoi (line, nullptr, 10);
            if (bingze != i) {
                std::cout << "\nbingze error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            }
            this->stalks.push_back (new PDF(InPut));
        }
        // read in seeds
        for (int i = 0; i < this->num_of_BDomains; i++) {   // initiation
            for (int j = 0; j < this->num_of_BDomains; j++) {
                this->seeds.push_back(NULL);
            }
        }
        for (int i = 1; i <= this->num_of_BDomains; i++) {
            for (int j = 1; j <= this->num_of_BDomains; j++) {
                std::getline (InPut, line);
                this->gene_edit(i, j, new PDF(InPut));
            }
        }
    }

    pumpkin (int num_of_binding_domains) {  // initiate a pumpkin
        this->num_of_BDomains = num_of_binding_domains;
        for (int i = 0; i < num_of_binding_domains; i++) {
            for (int j = 0; j < num_of_binding_domains; j++) {
                this->seeds.push_back(NULL);
            }
        }
    }

    void fill_pumpkin (G_ligand ligand) {
        for (int i = 0; i < ligand.linker.size(); i++) {    // build linkers
            this->stalks.push_back (new PDF('l', ligand.linker[i].L_c, ligand.linker[i].L_p));
        }
        for (int i = 1; i <= ligand.domain.size(); i++) {   // build BDs
            this->gene_edit(i, i, new PDF('d', ligand.domain[i-1].diam));
        }
        for (int j = 1; j <= ligand.linker.size(); j++) {   // process the first time
            this->gene_edit(j, j+1, new PDF('c', this->get_seed(j, j), this->stalks[j-1]));
            this->gene_edit(j+1, j, new PDF('c', this->get_seed(j+1, j+1), this->stalks[j-1]));
        }
        for (int i = 2; i < ligand.domain.size(); i++) {    // i = number of right binding domain - number of left binding domain
            for (int j = 1; j <= ligand.domain.size()-i; j++) {      // j = number of right binding domain
                // ligand.domain[j] --> ligand.domain[j+i]
                PDF iron_man('c', this->get_seed(j+i-1, j), this->stalks[j+i-2]);
                this->gene_edit(j, j+i, new PDF('c', this->get_seed(j, j), &iron_man));
                this->gene_edit(j+i, j, new PDF('c', this->get_seed(j+i, j+i), &iron_man));
            }
        }
        return;
    }

    void output(std::ostream &OutPut) { // output PDF to stream (file)
        OutPut << this->num_of_BDomains << "\n";    // the first line is an integer, refers to the length of the ligand
        for (int i = 0; i < this->stalks.size(); i++) { // then follows the stalks data
            OutPut << i << "\n";                    // one line of the number of stalk
            this->stalks[i]->output(OutPut);        // then the PDF output
        }
        for (int i = 1; i <= this->num_of_BDomains; i++) {  // then the seeds
            for (int j = 1; j <= this->num_of_BDomains; j++) {
                OutPut << i << " \t" << j << "\n";  // each seed will start with a line of indexs
                this->get_seed(i, j)->output(OutPut);
            }
        }
    }

    PDF* get_seed (int x, int y) {            // the PDF between x num BD & the linker before y num BD        x,y: 1, 2, 3, 4, 5, ...
        return this->seeds[(x-1)*this->num_of_BDomains + y-1];
    }

    void gene_edit (int x, int y, PDF* value) { // set PDF between x num BD & the linker before y num BD      x,y: 1, 2, 3, 4, 5, ...
        this->seeds[(x-1)*this->num_of_BDomains + y-1] = value;
    }

    int size(void) {
        return this->num_of_BDomains;
    }

    private:
    int num_of_BDomains;    // length of one dimention
};

class squash_set {
    public:
    std::vector<pumpkin> set_of_squash;

    squash_set (void) {
        // declare an empty squash_set
    }

    void squash_plant (std::vector<G_ligand>& Ligand_set) {
        for (int i = 0; i < Ligand_set.size(); i++) {   // traverse all sets of ligands
            std::cout << "preprocessing ligand #" << i << "\n";
            this->set_of_squash.push_back (pumpkin(Ligand_set[i].domain.size()));
            this->set_of_squash.back().fill_pumpkin(Ligand_set[i]);
        }
    }

    void output (void) {    // test codes
        std::cout << "\nprinting preprocessed data:\n";
        for (int i = 0; i < this->set_of_squash.size(); i++) {
            std::cout << "pumpkin  " << i <<":\n";
            for (int j = 1; j <= this->set_of_squash[i].size(); j++) {
                for (int k = 1; k <= this->set_of_squash[i].size(); k++) {
                    std::cout << this->set_of_squash[i].get_seed(j, k) << "\t";
                }
                std::cout << "\n";
            }
        }
    }

    void output (std::string dir) {    // output squash_set into streams (files)
        std::cout << "\nsaving preprocessed data:\n";
        for (int i = 0; i < this->set_of_squash.size(); i++) {
            std::cout << "pumpkin  " << i <<":\n";
            std::ofstream file_out (dir + "/temp/Pumpkin_" + std::to_string(i) + ".ppd");
            this->set_of_squash[i].output (file_out);
            file_out.close();
        }
        std::cout << "finished\n";
    }

    void input (std::string dir, int num_of_ligands) {      // input squash_set from streams (files)
        std::cout << "\nreading preprocessed data from files:\n";
        for (int i = 0; i < num_of_ligands; i++) {
            std::cout << "reading in one ligand\n";
            std::ifstream file_in (dir + "/temp/Pumpkin_" + std::to_string(i) + ".ppd");
            this->set_of_squash.push_back (pumpkin(file_in));
            file_in.close();
        }
            
        std::cout << "finished\n";
    }
};

#endif