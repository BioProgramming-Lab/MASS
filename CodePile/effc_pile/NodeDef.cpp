#ifndef NodeDef
#define NodeDef

#include "PDFdef.cpp"
#include <iostream>
#include <vector>
#include <cmath>

class Node {
    public:
        // PDF pdf;
        bool is_passed = false;
        char type;
        double diam;
        double L_c;
        double L_p;
        int true_type = -1; // true type of the ligand
        int label = -1; // label starts from 0. binding domains and linkers are respectively labelled

        std::vector<Node*> friends;

        Node (char type, double diam, int label, int true_type) {
            if (type == 'd') {
                this->type = type;
                this->diam = diam;
                this->label = label;
                this->true_type = true_type;
                return;
            }
            std::cout << "\nwrong type in creating node\n";
        }
        Node (char type, double L_c, double L_p, int label, int true_type) {
            if (type == 'l') {
                this->type = type;
                this->L_c = L_c;
                this->L_p = L_p;
                this->label = label;
                this->true_type = true_type;
                return;
            }
            std::cout << "\nwrong type in creating node\n";
        }


    private:;


};

inline void connect (Node* a, Node* b) {
    a->friends.push_back(b);
    b->friends.push_back(a);
}

#endif