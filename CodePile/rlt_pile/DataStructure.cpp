#include <vector>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <fstream>

int LigandTypeNumber = 1;

class Ligand {
    public:
        int length;
        int type;
        int merging_factor;     // the true TYPE of ligand. in this code, each type of ligands are duplicated for the sake of simulating several different ligands bind to a single receptor.
        
        Ligand (int length) {
            this->type = LigandTypeNumber;
            LigandTypeNumber++;
            this->length = length;
            this->merging_factor = -1;
        }

        void set_merging_factor(int factor) {
            if (this->merging_factor != -1) {
                std::cout << "error 3\n";
            }
            this->merging_factor = factor;
        }
};

class ligandLocus {
    public:
        int LigandType;
        int DomainNumber;       // from 1 to length

        ligandLocus(int type, int num);
        ligandLocus();
        bool isVacant(void) {
            if (this->LigandType == -1)
                return true;
            else    return false;
        }
        bool isSame(ligandLocus target) {
            if (this->LigandType == target.LigandType && this->DomainNumber == target.DomainNumber)
                return true;
            else    return false;
        }
        int get_merge_factor(std::vector<Ligand*>& suite) {
            return suite[this->LigandType-1]->merging_factor;
        }
};

ligandLocus::ligandLocus(int type, int num) {
    this->DomainNumber = num;
    this->LigandType = type;
}

ligandLocus::ligandLocus() {
    this->DomainNumber = -1;    // -1 represents to vacancy
    this->LigandType = -1;
}

class Set_of_Ligands {
    public:
        std::vector<Ligand*> suite;

        Set_of_Ligands (std::vector<Ligand*> suite) {   //------------------------------------------- to make a commandline software, here should be modified to be able to accept unlimited number of parameters
            this->suite = suite;
        }

        Set_of_Ligands () {
            // just do nothing
            // this makes it possible to claim a new variable of this class without any parameter.
        }
};

int StateNumber = 1;

class MicroState {
    public:
        std::vector<MicroState*> fathers;
        std::vector<MicroState*> sons;
        int stateNUM;
        int fatherNUM;
        int sonNUM;
        std::vector<ligandLocus> locus;
        MicroState(int length, MicroState* father);
        MicroState(int length, MicroState* father, std::vector<ligandLocus> locus);

        void show_state (void) {
            for (int i = 0; i < this->locus.size(); i++) {
                this->show_a_locus(i);
                std::cout << " ";
            }
        }
        void show_state_simplified (std::ofstream* file_out, int factor_num, Set_of_Ligands suite, int rec_length) {
            
            auto A_factorLocus = this->locus;   // copy the template

            std::vector<int> type_list_sub;     // subs of template
            for (int i = 1; i <= factor_num; i++) {     // traverse the factor
                type_list_sub.clear();
                for (int j = 0; j < this->length; j++) {   // traverse the locus
                    if (this->locus[j].isVacant()) {
                        continue;
                    }
                    if (this->locus[j].get_merge_factor(suite.suite) == i) {
                        int pos = -1;   // the position of the certain type
                        for (int i_sub = 0; i_sub < type_list_sub.size(); i_sub++) {
                            if (this->locus[ type_list_sub[i_sub] ].LigandType == this->locus[j].LigandType) {
                                if (i_sub == 0)
                                    A_factorLocus[j].LigandType = i;
                                else
                                    A_factorLocus[j].LigandType = factor_num + (i-1) * (rec_length-1) + i_sub;
                                pos = i_sub;
                                break;
                            }
                        }
                        if (pos == -1) {
                            type_list_sub.push_back(j);
                            if (type_list_sub.size() == 1)
                                A_factorLocus[j].LigandType = i;
                            else
                                A_factorLocus[j].LigandType = factor_num + (i-1) * (rec_length-1) + type_list_sub.size()-1;
                        }
                    }                    
                }
            }

            for (int i = 0; i < this->locus.size(); i++) {
                *file_out << A_factorLocus[i].LigandType;
                *file_out << "\t";
                *file_out << A_factorLocus[i].DomainNumber;
                *file_out << "\t";
            }/*
            for (int i = 0; i < this->locus.size(); i++) {
                *file_out << this->locus[i].LigandType;
                *file_out << "\t";
                *file_out << this->locus[i].DomainNumber;
                *file_out << "\t";
            }*/
        }
        void add_son (MicroState* son) {    // don't use this outer!  and need to add error check --------------------------------------
            if (this->is_son_contained(son)) {
                return;
            }
            this->sons.push_back(son);
            this->sonNUM++;
        }
        void add_fathers (std::vector<MicroState*>& father) {     // don't use this outer! one leaf node might have multiple fathers, so we need a vector here.
            for (int i = 0; i < father.size(); i++){
                if (this->is_father_contained(father[i]))
                    continue;
                this->fatherNUM ++;
                this->fathers.push_back(father[i]);
            }
        }
        void give_fathers_to (MicroState* best_friend) {    // it's logical to give fathers to your friend right? XD
            // std::cout << this->stateNUM << " -> " << best_friend->stateNUM << "\n";
            best_friend->add_fathers(this->fathers);
            for (auto i : this->fathers) {
                i->add_son(best_friend);
            }
            this->disappear();
        }
        bool isSymmetrical (MicroState* target) {
            if (this->locus.size() != this->length || this->length != target->length || target->locus.size() != target->length) {   //just check
                std::cout << "error 2";
            }
            for (int i = 0; i < this->length; i++) {
                if (! this->locus[i].isSame(target->locus[this->length-1-i]))
                    return false;
            }
            std::cout << "find" << this->stateNUM << " && " << target->stateNUM << "\n";
            return true;
        }
        bool isSame (MicroState* target) {
            if (this->locus.size() != this->length || this->length != target->length || target->locus.size() != target->length) {   // just check
                std::cout << "error 2";
            }
            for (int i = 0; i < this->length; i++) {
                if (! this->locus[i].isSame(target->locus[i]))
                    return false;
            }
            //std::cout << "find" << this->stateNUM << " && " << target->stateNUM << "\n";
            return true;
        }
        bool isFactorSame (MicroState* target, int factor_num, Set_of_Ligands suite) {
            if (this->locus.size() != this->length || this->length != target->length || target->locus.size() != target->length) {   // just check
                std::cout << "error 2";
            }
            for (int i = 0; i < this->length; i++) {    // proof filter
                if (this->locus[i].isVacant() != target->locus[i].isVacant())
                    return false;
            }

            auto A_factorLocus = this->locus;   // copy the template
            auto B_factorLocus = target->locus;   

            std::vector<int> type_list_sub;     // subs of template
            for (int i = 1; i <= factor_num; i++) {     // traverse the factor
                type_list_sub.clear();
                for (int j = 0; j < this->length; j++) {   // traverse the locus
                    if (this->locus[j].isVacant()) {
                        continue;
                    }
                    if (this->locus[j].get_merge_factor(suite.suite) == i) {
                        int pos = -1;   // the position of the certain type
                        for (int i_sub = 0; i_sub < type_list_sub.size(); i_sub++) {
                            if (this->locus[ type_list_sub[i_sub] ].LigandType == this->locus[j].LigandType) {
                                A_factorLocus[j].LigandType = i_sub + 1;
                                pos = i_sub;
                                break;
                            }
                        }
                        if (pos == -1) {
                            type_list_sub.push_back(j);
                            A_factorLocus[j].LigandType = type_list_sub.size();
                        }
                    }                    
                }
            }

            for (int i = 1; i <= factor_num; i++) {     // traverse the factor
                type_list_sub.clear();
                for (int j = 0; j < this->length; j++) {   // traverse the locus
                    if (target->locus[j].isVacant()) {
                        continue;
                    }
                    if (target->locus[j].get_merge_factor(suite.suite) == i) {
                        int pos = -1;   // the position of the certain type
                        for (int i_sub = 0; i_sub < type_list_sub.size(); i_sub++) {
                            if (target->locus[ type_list_sub[i_sub] ].LigandType == target->locus[j].LigandType) {
                                B_factorLocus[j].LigandType = i_sub + 1;
                                pos = i_sub;
                                break;
                            }
                        }
                        if (pos == -1) {
                            type_list_sub.push_back(j);
                            B_factorLocus[j].LigandType = type_list_sub.size();
                        }
                    }
                }
            }

            for (int i = 0; i < A_factorLocus.size(); i++) {
                if (this->locus[i].isVacant())
                    continue;
                if (! A_factorLocus[i].isSame(B_factorLocus[i])) {
                    return false;
                }
                if (this->locus[i].get_merge_factor(suite.suite) != target->locus[i].get_merge_factor(suite.suite))
                    return false;
            }
            return true;
        }
    private:
        int length;
        void show_a_locus (int num) {
            if (this->locus[num].LigandType == -1)
                std::cout << "_";
            else
                std::cout << this->locus[num].LigandType;
            std::cout << "&";
            if (this->locus[num].DomainNumber == -1)
                std::cout << "_";
            else
                std::cout << this->locus[num].DomainNumber;
        }

        void disappear (void) { // disappear from its fathers sight
            for (auto father : this->fathers) { // traverse all fathers
                for (int i = 0; i < father->sons.size(); i++) { // traverse the sons of the father
                    if (father->sons[i] == this) {
                        father->sons.erase(father->sons.begin()+i, father->sons.begin()+i+1);
                        father->sonNUM --;
                        break;
                    }
                }
            }
        }

        bool is_son_contained(MicroState* son) {
            for (auto i : this->sons) {
                if (son == i) {
                    return true;
                }
            }
            return false;
        }

        bool is_father_contained(MicroState* father) {
            for (auto i : this->fathers) {
                if (father == i) {
                    return true;
                }
            }
            return false;
        }
};

MicroState::MicroState(int length, MicroState* father) {
    this->length = length;
    this->stateNUM = StateNumber;
    StateNumber ++;
    this->fatherNUM = 1;
    this->sonNUM = 0;
    this->fathers.push_back(father);
    this->locus.resize(length, ligandLocus());
}

MicroState::MicroState(int length, MicroState* father, std::vector<ligandLocus> locus) {
    this->length = length;
    this->stateNUM = StateNumber;
    StateNumber ++;
    this->fatherNUM = 1;
    this->sonNUM = 0;
    this->fathers.push_back(father);
    this->locus = locus;
}

class relations {
    public:
        MicroState* head;   // this is a pointor to the first node
        int merge(int layer);
        relations(int length, Set_of_Ligands ligandsSuite);
        int GenerateRelation(std::string RLT_type);
        void show_relations(void) {
            std::cout << "\nshowing relations:\n";
            this->show_relations_recursion(this->head, 1);
            std::cout << "relations shown";
        }
        void show_relations_in_dot(std::string dir) {
            std::ofstream file_out(dir+"relations.dot");
            file_out << "digraph G {\n";            
            std::vector<MicroState*> nodes;
            for (int layer = 1; layer <= this->length+1; layer++) {
                nodes.clear();
                this->getLayerFathers(nodes, layer);
                for (auto j : nodes) {
                    for (int i = 0; i < j->sonNUM; i++) {    // output sons
                        file_out << "\t" << j->stateNUM << " -> " << j->sons[i]->stateNUM << " ;";
                        file_out << "\n";
                        file_out << "\t" << j->sons[i]->stateNUM << " -> " << j->stateNUM << " ;";
                        file_out << "\n";
                    }
                }
            }
            file_out << "}";    
            file_out.close();
        }
        void show_relations_in_unframed_son_dot(std::string dir) { // from father to son only, and no frame
            std::ofstream file_out (dir+"relations.udot");
            
            std::vector<MicroState*> nodes;
            for (int layer = 1; layer <= this->length+1; layer++) {
                nodes.clear();
                this->getLayerFathers(nodes, layer);
                for (auto j : nodes) {
                    for (int i = 0; i < j->sonNUM; i++) {    // output sons
                        file_out << j->stateNUM << "\t" << j->sons[i]->stateNUM;
                        file_out << "\n";
                    }
                }
            }
            file_out.close();
        }
        void show_states_in_file(std::string dir) {    // in each line, the first number is the stateNUM, followed by its state.
            std::ofstream file_out (dir+"states.states");
            std::vector<MicroState*> nodes;
            file_out << "# file explanation:" << "\n# the different microstates are stored in this file\n";
            file_out << "# lines with pound at first will be ignored\n";
            file_out << "# the rest lines form as: stateNumber [states].\n";
            file_out << "# each states are represented with two numbers, ligand type and domain number.\n";
            for (int layer = 1; layer <= this->length+1; layer++) {
                nodes.clear();
                this->getLayerFathers(nodes, layer);
                for (auto j : nodes) {
                    file_out << j->stateNUM << "\t\t";
                    j->show_state_simplified(&file_out, this->factor_num, this->ligandsSuite, this->length);
                    file_out << "\n";
                }
            }
            file_out.close();
        }
        void show_relations_in_file(std::string dir) {
            std::ofstream file_out (dir+"relations.rlt");

            std::vector<MicroState*> nodes;
            file_out << "# file explanation:" << "\n# the relationship between different microstates are stored in this file\n";
            file_out<< "# lines with pound at first will be ignored\n";
            file_out << "# the rest lines form as: stateNumber S [stateNumbers] F [stateNumbers].\n";
            for (int layer = 1; layer <= this->length+1; layer++) {
                nodes.clear();
                this->getLayerFathers(nodes, layer);
                for (auto j : nodes) {
                    file_out << j->stateNUM;
                    file_out << "\tS";
                    for (int i = 0; i < j->sonNUM; i++) {    // output sons
                        file_out << "\t" << j->sons[i]->stateNUM;
                    }
                    file_out << "\tF";
                    if (j->stateNUM == 1) {
                        file_out << "\n";
                        continue;
                    }
                    for (int i = 0; i < j->fatherNUM; i++) {    // output sons
                        file_out << "\t" << j->fathers[i]->stateNUM;
                    }
                    file_out << "\n";
                }
            }
            file_out.close();
        }

    private:
        int length;
        int factor_num;
        int GenerateLayer(int layer);
        Set_of_Ligands ligandsSuite;

        void show_relations_recursion(MicroState* nowstate, int step) {
            if (step > this->length) {
                return;
            }
            std::cout << "stateNUM:" << nowstate->stateNUM << "\t state:";
            nowstate->show_state();
            for (int i = 0; i < nowstate->sonNUM; i++) {    // output sons
                std::cout << "\n\t";
                std::cout << "SonStateNUM:" << nowstate->sons[i]->stateNUM << "\t SonState:";
                nowstate->sons[i]->show_state();
            }
            std::cout << "\n";
            if (step < this->length)
                for (int i = 0; i < nowstate->sonNUM; i++) {    // recursion
                    this->show_relations_recursion(nowstate->sons[i], step+1);                
                }
        }

        void duplicate_ligands(void) {  // 4-3-3 --> 1,2,1,1,1,2,2,2
            int original_num = this->ligandsSuite.suite.size();
            this->factor_num = original_num;
            for (int i = 0; i < original_num; i++) {
                this->ligandsSuite.suite[i]->set_merging_factor(this->ligandsSuite.suite[i]->type); // set the merging factor of original ligands first
                for (int j = 1; j < this->length; j++) {    // duplicate length-1 times
                    this->ligandsSuite.suite.push_back(new Ligand(this->ligandsSuite.suite[i]->length));
                    this->ligandsSuite.suite.back()->set_merging_factor(this->ligandsSuite.suite[i]->merging_factor);
                }
            }
        }
    
        void getLayerFathers(std::vector<MicroState*>& Fathers, int layer) {    // get the nodes of the layer. the parameter "Fathers" requires a container for the output.
            std::set<MicroState*> father_set;
            getLayerFathers_recursion(father_set, layer, 1, this->head);
            for (auto i = father_set.begin(); i != father_set.end(); i++) {
                Fathers.push_back(*i);
            }
        }

        void getLayerFathers_recursion(std::set<MicroState*>& father, int layer, int step, MicroState* nowstate) {
            if (step == layer) {
                father.insert(nowstate);
                return;
            }
            if (step > layer) {
                std::cout << "exception 1\n";   // error when step comes over certain layer
            }
            for (auto i : nowstate->sons) {     // a new gramma in C++11, hopefully it will work
                getLayerFathers_recursion(father, layer, step+1, i);
            }
        }

        int factor_merge(int layer) {
            std::vector<MicroState*> nodes;
            std::vector<int> tag;
            this->getLayerFathers(nodes, layer);
            tag.resize(nodes.size(), 0);
            int tag_num = 1;
            
            for (int i = 0; i < nodes.size(); i++) {    // set same tags for microstates to merge. tags start from 1, and 0 means no tag.
                if (tag[i] != 0)
                    continue;
                bool key = false;
                for (int j = i+1; j < nodes.size(); j++) {
                    if (nodes[i]->isFactorSame(nodes[j], this->factor_num, this->ligandsSuite)) {
                        if(key == false){
                            tag[i] = tag_num;
                            key = true;
                        }
                        tag[j] = tag_num;
                    }
                }
                if (key == true)
                    tag_num++;
            }

            for (int i = 1; i < tag_num; i++) { // traverse all tags to merge
                // correct the relations then
                int key = -1;
                if (i == 12){
                    i = 12;
                }
                for (int j = 0; j < nodes.size(); j++) {
                    if (tag[j] == i) {
                        if (key == -1){
                            key = j;
                            continue;
                        }
                        else {
                            nodes[j]->give_fathers_to(nodes[key]);
                            delete nodes[j];
                            nodes.erase(nodes.begin()+j, nodes.begin()+j+1);
                            tag.erase(tag.begin()+j, tag.begin()+j+1);
                            j--;
                        }
                    }
                }
            }
            return 0;
        }

        void stateNUM_reorder (void) {  // reorder the stateNUM after merging
            std::queue<MicroState*> processing;
            std::set<MicroState*> book; // use set to make sure no duplication in the queue
            int numberCount = 1;
            processing.push(this->head);
            book.insert(this->head);
            while (! processing.empty()) {
                for (auto i : processing.front()->sons) {
                    if ( book.find(i) == book.end() ){
                        processing.push(i);
                        book.insert(i);
                    }
                }
                processing.front()->stateNUM = numberCount;
                book.erase(processing.front());
                processing.pop();
                
                numberCount++;
            }
        }
};

int relations::merge(int layer) {
    // this->symmetrical_merge(layer);
    // the above func is not finished yet
    this->factor_merge(layer);
    // more merge choice could be done here ------------------------------------------------------------------------
    return 0;
}

relations::relations(int length, Set_of_Ligands ligandsSuite) {
    this->head = new MicroState(length, NULL);
    this->length = length;
    this->ligandsSuite.suite = ligandsSuite.suite;
}

int relations::GenerateRelation(std::string RLT_type) {
    this->factor_num = 0;
    if (RLT_type == "full")         // we do not duplicate ligands when connection should be simplified, which means some microstates are removed. Here, when RLT_type is "full", we duplicate ligands.
        this->duplicate_ligands();  // duplicate ligands first
    for(int i = 1; i <= this->length; i++) {
        this->GenerateLayer(i);
        this->merge(i+1);
    }
    this->stateNUM_reorder();
    //this->same_merge(this->length + 1);   not needed anymore, because factor-merge is a universal substitution.
    return 0;   // need proper return ---------------------------------------------------------
}

bool isContained(std::vector<ligandLocus>& usedLocus, ligandLocus locus) {
    for (auto i : usedLocus) {
        if (locus.isSame(i))
            return true;
    }
    return false;
}

int relations::GenerateLayer(int layer) {
    // get the array of fathers
    std::vector<MicroState*> fathers;
    this->getLayerFathers(fathers, layer);
    
    for (auto nowstate : fathers) {    // traverse all the father nodes
        // for each father nodes, we'll generate its son nodes here.
        std::vector<ligandLocus> locusSet = nowstate->locus;
        std::vector<ligandLocus> usedLocus;
        for (int i = 0; i < locusSet.size(); i++) {     // get the domains already used
            if (! locusSet[i].isVacant())
                usedLocus.push_back(locusSet[i]);
        }
        for (int i = 0; i < locusSet.size(); i++) {   // traverse all the locuses to find a vacancy
            if (locusSet[i].isVacant()) {
                
                for (int theLigandNum = 0; theLigandNum < this->ligandsSuite.suite.size(); theLigandNum++) {    //traverse the ligand suite
                    auto nowLigand = this->ligandsSuite.suite[theLigandNum];
                    for (int domain = 1; domain <= nowLigand->length; domain++) {   // traverse the domains to find a proper one
                        if (isContained(usedLocus, ligandLocus(nowLigand->type, domain)))
                            continue;   // we don't use the same domain already used
                        // create new node and add it to the tree
                        locusSet[i] = ligandLocus(nowLigand->type, domain);
                        MicroState* newMicrostate = new MicroState(this->length, nowstate, locusSet);
                        nowstate->add_son(newMicrostate);
                        locusSet[i] = ligandLocus();    // back to the origin
                    }
                }
            }
        }
    }
    return 0;   // need proper return -----------------------------------------------
}

