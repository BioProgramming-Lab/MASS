#ifndef GraphStructure
#define GraphStructure

#include <vector>

class G_Domain {
    public:
        double diam;
        G_Domain (double diam) {
            this->diam = diam;
        }
};

class G_Linker {
    public:
        double L_c;
        double L_p;
        G_Linker (double L_c, double L_p) {
            this->L_c = L_c;
            this->L_p = L_p;
        }
};

class G_ligand {
    public:
        std::vector<G_Domain> domain;   // 0 to n
        std::vector<G_Linker> linker;   // 0 to n-1
};

class G_ligandLocus {
    public:
        int LigandType;
        int DomainNumber;       // from 1 to length

        G_ligandLocus(int type, int num);
        G_ligandLocus();
        bool isVacant(void) {
            if (this->LigandType == -1)
                return true;
            else    return false;
        }
        bool isSame(G_ligandLocus target) {
            if (this->LigandType == target.LigandType && this->DomainNumber == target.DomainNumber)
                return true;
            else    return false;
        }
};

G_ligandLocus::G_ligandLocus(int type, int num) {
    this->DomainNumber = num;
    this->LigandType = type;
}

G_ligandLocus::G_ligandLocus() {
    this->DomainNumber = -1;    // -1 represents to vacancy
    this->LigandType = -1;
}

class G_microstate {
    public:
        std::vector<G_ligandLocus> holder;

};

#endif