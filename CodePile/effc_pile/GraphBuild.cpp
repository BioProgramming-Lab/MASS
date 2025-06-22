#ifndef GraphBuild
#define GraphBuild

#include "GraphStructure.cpp"
#include "NodeDef.cpp"
#include "Pumpkin.cpp"
#include <vector>
#include <set>
#define AVOGADRO 0.0006022
// AVOGADRO = 6.022e23 / (1e9)^3,    (1e9)^3 is the constant when you transform Angstrom^3 to Liter

class NodesGraph {
    public:
        Node* origin;
        Node* Des_ligand;
        Node* Des_receptor;

        void test(void) {   // just for test, don't use it
            auto a = this->get_PDF_to_origin(this->Des_ligand, this->origin);
            freopen("testPDF.txt", "w", stdout);
            a->show_x_array();
            a->show_y_array();
        }

        void destory(void) {    // recycle the memory
            for (auto i : this->recycle_bin_node) {
                delete i;
            }
            for (auto i : this->recycle_bin_PDF) {
                delete i;
            }
        }

        double effc(void) {  // calculate the effective concentration
            auto a = this->get_PDF_to_origin(this->Des_receptor, this->origin);
            auto b = this->get_PDF_to_origin(this->Des_ligand, this->origin);

            if (a->not_exist || b->not_exist)
                std::cout << "graph build error\n";

            a->normalize();
            b->normalize();
            /*freopen("testPDF.txt", "w", stdout);
            a->show_x_array();
            a->show_y_array();*/
            int length = (a->length() < b->length() ? a->length(): b->length());
            double sum = 0;  // P_OneFreeEnd2OneFreeEnd
            for (int i = 1; i < length-1; i++) {
                sum += a->PDFarray[i]*b->PDFarray[i] *i*i;  // might over limit of double if accuracy changed
            }
            sum += a->PDFarray[length-1]*b->PDFarray[length-1]/2.0 *(length-1)*(length-1);  // devide by 2 because of the trapezoid rule
            sum *= 4*PI * Interval*Interval*Interval;

            this->destory();    // recycle the memory
            return sum/AVOGADRO;
        }

        double effc(squash_set* PPD) {  // calculate the effective concentration
            auto a = this->get_PDF_to_origin_with_a_squash(this->Des_receptor, this->origin, PPD);
            auto b = this->get_PDF_to_origin_with_a_squash(this->Des_ligand, this->origin, PPD);

            if (a->not_exist || b->not_exist)
                std::cout << "graph build error\n";

            a->normalize();
            b->normalize();
            /*freopen("testPDF.txt", "w", stdout);
            a->show_x_array();
            a->show_y_array();*/
            int length = (a->length() < b->length() ? a->length(): b->length());
            double sum = 0;  // P_OneFreeEnd2OneFreeEnd
            for (int i = 1; i < length-1; i++) {
                sum += a->PDFarray[i]*b->PDFarray[i] *i*i;  // might over limit of double if accuracy changed
            }
            sum += a->PDFarray[length-1]*b->PDFarray[length-1]/2.0 *(length-1)*(length-1);  // devide by 2 because of the trapezoid rule
            sum *= 4*PI * Interval*Interval*Interval;

            this->destory();    // recycle the memory
            return sum/AVOGADRO;
        }

        Node* BuildGraph (G_microstate* microstate1, G_microstate* microstate2, std::vector<G_ligand>& Ligand_set) {    // ligand set with receptor on [0], but no duplications in
            // from m1 add sth to get m2
            // find out if graph build is in need
            // and initiate some essence parameters in the class
            this->total_num_of_ligand = Ligand_set.size()-1;
            this->receptor_length = microstate1->holder.size();
            
            int origin_num = -1;
            int des_num = -1;
            for (int i = 0; i < microstate1->holder.size(); i++) {
                if (microstate1->holder[i].isVacant())
                if (!microstate2->holder[i].isVacant()) {
                    for (int j = 0; j < microstate1->holder.size(); j++) {
                        if (microstate2->holder[j].LigandType == microstate2->holder[i].LigandType && i != j) {   // find the origin
                            origin_num = j;
                            des_num = i;
                            break;  // need to build graph
                        }
                    }
                    if (origin_num == -1)
                        return NULL;    // no need to build graph
                    break;
                }
            }

            // build initial receptor and ligands net
            std::vector<Node*> heads;
            int mapping_from_ligand_type_to_heads_sub[100] = {0}; //----------------------------------could be done better, here might occurs segment exceptions
            heads.push_back(this->create_ori_net(0, Ligand_set));   // build up receptor net
            std::set<int> book;     // book in the first time
            std::set<int> book2;    // book the other time. with this two books, we can get rid of those free end ligands in microstate2
            for (int i = 0; i < microstate1->holder.size(); i++) {
                if (microstate2->holder[i].isVacant() || book2.find(microstate2->holder[i].LigandType) != book2.end()) {
                    continue;
                }
                if (book.find(microstate2->holder[i].LigandType) == book.end()) {   // didn't find
                    book.insert(microstate2->holder[i].LigandType);
                    continue;
                }
                else {  // find in book1 but not in book2
                    book2.insert(microstate2->holder[i].LigandType);
                    // build ligand net
                    mapping_from_ligand_type_to_heads_sub[microstate2->holder[i].LigandType] = heads.size();
                    heads.push_back(this->create_ori_net(this->get_true_type(microstate2->holder[i].LigandType), Ligand_set));
                }
            }

            // set origin and destiny
            this->Des_receptor = get_num_of_domain_through_0friend(heads[0], des_num);
            this->Des_ligand = get_num_of_domain_through_0friend(heads[mapping_from_ligand_type_to_heads_sub[ microstate2->holder[des_num].LigandType ]], microstate2->holder[des_num].DomainNumber-1);
            this->origin = get_num_of_domain_through_0friend(heads[0], origin_num);

            // connect each of ligands with receptor;
            for (int i = 0; i < microstate1->holder.size(); i++) {
                if (microstate1->holder[i].isVacant() || book2.find(microstate2->holder[i].LigandType) == book2.end()) // we don't connect them if vacant in microstate1 or didn't find in book2
                    continue;
                connect(get_num_of_domain_through_0friend(heads[0], i), get_num_of_domain_through_0friend(heads[mapping_from_ligand_type_to_heads_sub[ microstate2->holder[i].LigandType ]], microstate2->holder[i].DomainNumber-1));
            }

            return this->origin;
        }
    
    private:
        int total_num_of_ligand;
        int receptor_length;
        std::vector<Node*> recycle_bin_node;    // stop this program from memory leaking
        std::vector<PDF*> recycle_bin_PDF;
        Node* get_num_of_domain_through_0friend(Node* a, int step) {            // a very dangrous func, don't use it unless you read it thoroughly
            if (step <= 0)
                return a;
            Node* labor = a;
            labor = labor->friends[0];
            labor = labor->friends[1];
            for (int i = 1; i < step; i++) {
                labor = labor->friends[1];
                labor = labor->friends[1];  // again because of the linker between.
            }
            return labor;
        }
        Node* create_ori_net(int num, std::vector<G_ligand>& Ligand_set) {  // create a initial link
            Node* labor_new_node1 = new Node('d', Ligand_set[num].domain[0].diam, 0, num);  // the first domain
            this->recycle_bin_node.push_back(labor_new_node1);
            Node* labor_new_node2;
            Node* ans = labor_new_node1;
            for (int i = 1; i < Ligand_set[num].domain.size(); i++) {  
                labor_new_node2 = new Node('l', Ligand_set[num].linker[i-1].L_c, Ligand_set[num].linker[i-1].L_p, i-1, num);
                this->recycle_bin_node.push_back(labor_new_node2);
                connect(labor_new_node1, labor_new_node2);
                labor_new_node1 = new Node('d', Ligand_set[num].domain[i].diam, i, num);
                this->recycle_bin_node.push_back(labor_new_node1);
                connect(labor_new_node2, labor_new_node1);
            }
            return ans;
        }
        int get_true_type(int n) { // get the ture type, so that we can get its parameters from Ligand_set
            if (n <= this->total_num_of_ligand) {
                return n;
            }
            else {
                return (n-1-this->total_num_of_ligand)/(this->receptor_length-1) +1;
            }
        }
        
        PDF* get_PDF_to_origin(Node* now_node, Node* destiny) {  // recuirsion -> PDF. unchecked, and could be inproved with efficiency.
            //std::cout << now_node << "\t";
            if (now_node == destiny){
                //std::cout << "destiny\n";
                PDF* newP = new PDF('d', now_node->diam);
                this->recycle_bin_PDF.push_back(newP);
                return newP;    // destiny must be a domain
            }
            int avail_friends = 0;
            Node* my_friends[2];
            for (Node* i : now_node->friends) {
                if (!i->is_passed) {
                    my_friends[avail_friends] = i;
                    avail_friends++;
                }                    
            }
            if (avail_friends == 0) {
                //std::cout << "\n";
                PDF* newP = new PDF();
                this->recycle_bin_PDF.push_back(newP);
                return newP;
            }
            if (avail_friends == 1) {
                //std::cout << "->  " << my_friends[0] << "\n";
                now_node->is_passed = true;
                if (now_node->type == 'd') {
                    PDF* old_friend = get_PDF_to_origin(my_friends[0], destiny);
                    now_node->is_passed = false;
                    if (old_friend->not_exist)
                        return old_friend;
                    PDF* on_way_P = new PDF('d', now_node->diam);
                    PDF* newP = new PDF('c', old_friend, on_way_P);
                    this->recycle_bin_PDF.push_back(on_way_P);
                    this->recycle_bin_PDF.push_back(newP);
                    return newP;
                }
                if (now_node->type == 'l') {
                    PDF* old_friend = get_PDF_to_origin(my_friends[0], destiny);
                    now_node->is_passed = false;
                    if (old_friend->not_exist) 
                        return old_friend;
                    PDF* on_way_P = new PDF('l', now_node->L_c, now_node->L_p);
                    PDF* newP = new PDF('c', old_friend, on_way_P);
                    this->recycle_bin_PDF.push_back(on_way_P);
                    this->recycle_bin_PDF.push_back(newP);
                    return newP;
                }
            }
            if (avail_friends == 2) {
                //std::cout << "->  " << my_friends[0] << "    " << my_friends[1] << "\n";
                my_friends[0]->is_passed = true;     // cut down one way
                PDF* left_friend = get_PDF_to_origin(now_node, destiny);
                my_friends[0]->is_passed = false;    // recover
                my_friends[1]->is_passed = true;     // cut down the other
                PDF* right_friend = get_PDF_to_origin(now_node, destiny);
                my_friends[1]->is_passed = false;    // recover
                PDF* newP = new PDF('m', left_friend, right_friend);
                this->recycle_bin_PDF.push_back(newP);
                return newP;
            }
            else 
                std::cout << "\nmore than three connections?\n";
            PDF* newP = new PDF();
            this->recycle_bin_PDF.push_back(newP);
            return newP;
        }

        PDF* get_PDF_to_origin_with_a_squash(Node* now_node, Node* destiny, squash_set* PPD) {  // recuirsion -> PDF. with preprocessed data, namely PPD.
            //std::cout << now_node << "\t";
            if (now_node == destiny){
                //std::cout << "destiny\n";
                //PDF* newP = new PDF('d', now_node->diam);
                //this->recycle_bin_PDF.push_back(newP);
                //return newP;    // destiny must be a domain
                return PPD->set_of_squash[now_node->true_type].get_seed(now_node->label+1, now_node->label+1);
            }
            int avail_friends = 0;
            Node* my_friends[2];
            for (auto i : now_node->friends) {
                if (!i->is_passed) {
                    my_friends[avail_friends] = i;
                    avail_friends++;
                }                    
            }
            if (avail_friends == 0) {
                //std::cout << "\n";
                PDF* newP = new PDF();
                this->recycle_bin_PDF.push_back(newP);
                return newP;
            }
            if (avail_friends == 1) {

                if (now_node->type == my_friends[0]->type) {  // cross the binding site. (from ligand to receptor or reverse)
                    PDF* old_friend = get_PDF_to_origin_with_a_squash(my_friends[0], destiny, PPD);
                    if (old_friend->not_exist)
                        return old_friend;
                    PDF* newP = new PDF('c', 
                                        PPD->set_of_squash[now_node->true_type].get_seed( now_node->label+1, now_node->label+1 ), 
                                        old_friend);    // if two parameters in get_seed function equals to each other. it will return PDF* of the binding domain.
                    this->recycle_bin_PDF.push_back(newP);
                    return newP;
                }

                //std::cout << "->  " << my_friends[0] << "\n";
                Node* spider_man = now_node;    // our friendly neighbour spider man will help us to traverse the graph in a role of now_node
                std::vector<Node*> Pizza_delivery;  // spider man pizza delivery go back to the pizza store. this is used to recover is_passed term
                int neighbour = 0;
                while (true) {
                    spider_man->is_passed = true;
                    Pizza_delivery.push_back(spider_man);
                    neighbour = 0;
                    Node* Gotham[2];    // friends of spider man live in the Gotham City                    
                    for (auto i : spider_man->friends) {    // get dewellers in the Gotham city
                        if (!i->is_passed) {
                            Gotham[neighbour] = i;
                            neighbour++;
                        }
                    }                    
                    if (neighbour != 1 || (spider_man->type == Gotham[0]->type)) {
                        // use PPD to speed up linear convolution. (topology linear I mean)

                        PDF* old_friend = get_PDF_to_origin_with_a_squash(spider_man, destiny, PPD);
                        if (old_friend->not_exist) {
                            for (auto i : Pizza_delivery) {
                                i->is_passed = false;
                            }
                            return old_friend;
                        }
                            
                        PDF* newP = new PDF('c', 
                                            PPD->set_of_squash[spider_man->true_type].get_seed( now_node->label+1, spider_man->label+1 ), 
                                            old_friend);
                        this->recycle_bin_PDF.push_back(newP);
                        // unpass
                        for (auto i : Pizza_delivery) {
                            i->is_passed = false;
                        }
                        return newP;
                    }
                    spider_man = Gotham[0];
                }
            }
            if (avail_friends == 2) {
                //std::cout << "->  " << my_friends[0] << "    " << my_friends[1] << "\n";
                my_friends[0]->is_passed = true;     // cut down one way
                PDF* left_friend = get_PDF_to_origin_with_a_squash(now_node, destiny, PPD);
                my_friends[0]->is_passed = false;    // recover
                my_friends[1]->is_passed = true;     // cut down the other
                PDF* right_friend = get_PDF_to_origin_with_a_squash(now_node, destiny, PPD);
                my_friends[1]->is_passed = false;    // recover
                PDF* newP = new PDF('m', left_friend, right_friend);
                this->recycle_bin_PDF.push_back(newP);
                return newP;
            }
            else 
                std::cout << "\nmore than three connections?\n";
            PDF* newP = new PDF();
            this->recycle_bin_PDF.push_back(newP);
            return newP;
        }
};

#endif