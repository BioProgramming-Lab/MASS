// linux version

#include "GraphBuild.cpp"
#include "CellSurfaceEffc.cpp"
#include "PDFdef.cpp"
#include "Pumpkin.cpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <sys/stat.h>
#include <algorithm>
#include <filesystem>

int main (int argc, char* argv[]) { // file input all in this cpp
    std::string dir("!");
    std::string RLTfile_dir("/usr/src/CodePile/RLT_files/");
    std::string topology("!");
    int job_config = -1;
    bool silence = false;
    bool clean = false;
    bool out_relations = false;
    bool CellSurface_mode = false;
    bool PreProcess_switch = false;
    int PreProcess_job = -1;
    bool use_PPD_switch = false;
    int max_l = 2000;
    int max_s = 1000;
    
    for (int i = 1; i < argc; i++) {
        std::string arg_string(argv[i]);
        if (arg_string[0] == '-') {
            if (arg_string.compare("-d") == 0) {        // directory
                i++;
                std::string arg_string(argv[i]);
                dir = arg_string;
                continue;
            }
            else if (arg_string.compare("-t") == 0) {   // topology
                i++;
                std::string arg_string(argv[i]);
                topology = arg_string;
                continue;
            }
            else if (arg_string.compare("-s") == 0) {   // run in silence
                silence = true;
                continue;
            }
            else if (arg_string.compare("-j") == 0) {   // parallel mode, set single job to work with
                i++;
                std::string arg_string(argv[i]);
                job_config = std::stoi(arg_string, nullptr, 10);
                continue;
            }
            else if (arg_string.compare("-c") == 0) {   // clean temp files mode
                clean = true;
                continue;
            }
            else if (arg_string.compare("-R") == 0) {   // RLT file directory, default to be /usr/src/CodePile/RLT_files/
                i++;
                std::string arg_string(argv[i]);
                RLTfile_dir = arg_string;
                if (RLTfile_dir.back() != '/') {    // add '/' to the end of directory
                    RLTfile_dir += '/';
                }
                if (!std::filesystem::exists(RLTfile_dir)) {   // check if directory exists
                    std::cout << "ERROR: " << RLTfile_dir << " directory does not exist.\n";
                    exit(1);    // terminate with error, add do you wonna create directory? --------------------------------------------------------------------------------------
                }
                continue;
            }
            else if (arg_string.compare("--ham") == 0) {  // use preprocessed data in ./temp/
                use_PPD_switch = true;
                continue;
            }
            else if (arg_string.compare("--rlt") == 0) {  // output number of relations
                out_relations = true;
                continue;
            }
            else if (arg_string.compare("--cs") == 0) {  // cell surface mode
                CellSurface_mode = true;
                continue;
            }
            else if (arg_string.compare("--pp") == 0) {  // conduct preprocess and output preprocessed data
                PreProcess_switch = true;
                i++;
                std::string arg_string(argv[i]);
                if ( arg_string != "all" ) {
                    PreProcess_job = std::stoi(arg_string, nullptr, 10);
                }
                continue;
            }
            else if (arg_string.compare("--ml") == 0) {  // max length of possible contour length
                i++;
                std::string arg_string(argv[i]);
                max_l = std::stoi(arg_string, nullptr, 10);
                continue;
            }
            else if (arg_string.compare("--ms") == 0) {  // max splice in PDF
                i++;
                std::string arg_string(argv[i]);
                max_s = std::stoi(arg_string, nullptr, 10);
                continue;
            }
            else if (arg_string.compare("-v") == 0) {
                std::cout << "Geffc Version: Geffc_v3\n";  // set up standard version control---------------------------------------------------------------------------------------------------------
                return 0;
            }
            else if (arg_string.compare("-h") == 0) {   // help
                std::cout << "\n";
                std::cout << "usage: Geffc [-h][-v] ...\n";
                std::cout << "\n";
                std::cout << "necessary arguments:\t-d\t[dir]\tset direction of your work folder\n";
                std::cout << "                    \t-t\t[topology]\tset topology for your reaction system, e.g. 2-2\n";
                std::cout << "\n";
                std::cout << "optional arguments:\t-j\tcalculate effective concentration of a certain connection. (this is for parallel computing)\n";
                std::cout << "                   \t-s\trun in silence\n";
                std::cout << "                   \t-c\tremove temp files while integrating them into effcs.effc. (won't have any effect when [-j] is given)\n";
                std::cout << "                   \t-R\t[dir]\tset directory of RLT files. (default to be /usr/src/CodePile/RLT_files/)\n";
                std::cout << "                   \t--pp\tpre-process and store Pumpkin files in ./temp/\n";
                std::cout << "                   \t--ham\tuse pre-processed data store in ./temp/\n";
                std::cout << "                   \t--cs\tcell surface mode (each ligand on cell surface will be viewed as monomer without linkers between each other)\n";
                std::cout << "                   \t--ml\t[num]\talter the max permisive contour length of the whole structure, unit of angstrom, default to be 2000\n";
                std::cout << "                   \t--ms\t[num]\talter the splice in above \"ml\", default to be 1000. this setting is for accuracy controll\n";
                std::cout << "\n";
                std::cout << "commands:\t-h help\n";
                std::cout << "         \t-v version\n";
                std::cout << "         \t--rlt\tprint number of relations. (also require [-t] in this case)\n";
                std::cout << "\n";
                return 0;
            }
            else {
                std::cout << "unrecogonized argument [" << arg_string << "]\n";
                return 0;
            }
        }
    }

    if (out_relations) {    // output number of relations
        if (topology == "!") {
            std::cout << "[-t] argument is necessary in relations counting\n";
            return 0;
        }
        std::ifstream work_in_file (RLTfile_dir+topology+"relations.udot");
        if (not silence)
            std::cout << "number of relations in " << topology << "topology:\t" << std::count(std::istreambuf_iterator<char>(work_in_file), std::istreambuf_iterator<char>(), '\n') << "\n";
        else
            std::cout << std::count(std::istreambuf_iterator<char>(work_in_file), std::istreambuf_iterator<char>(), '\n');
        work_in_file.close();
        return 0;
    }

    if (topology == "!" || dir == "!") {    // necessary arguments missing
        std::cout << "[-d]&[-t] arguments are necessary\n";
        return 0;
        /*
        job_config = 73;
        dir = "multivalency/CodePile/julia_sim_pile/TrueSim/trivalentWS/";
        topology = "3-1-1-1";
        CellSurface_mode = true;
        */
    }

    if (job_config != -1 || PreProcess_switch) {         // parallel running on a certain job, whether if need to make a temp directory. We make temp folder also for pumpkins
        struct stat check_info;
        if (stat((dir+"/temp").c_str(), &check_info) != 0) {
            // mkdir((dir+"/temp").c_str());   // for windows
            mkdir((dir+"/temp").c_str(), 0777);  // for linux and mac
        }
    }

    PDFdef::const_alter(max_l, max_s);
    
    // microstates file read in
    std::ifstream states_file;
    states_file.open (RLTfile_dir+topology+"states.states");
    if (!states_file.is_open()) {
        std::cout << "failed to open states file\n";
        return 0;
    }
    int num_of_states = std::count(std::istreambuf_iterator<char>(states_file), std::istreambuf_iterator<char>(), '\n') - 5; // 5 lines of annotations
    G_microstate* microstates_set[num_of_states+2];   // microstates list
    for (int i = 0; i <= num_of_states+1; i++) {
        microstates_set[i] = NULL;
    }
    states_file.close();    // toy method to return to the beginning of the file
    states_file.open (RLTfile_dir+topology+"states.states");

    std::string line;
    int state_num;
    int lig_type;
    int lig_num;
    int number_of_states = 0;
    while (std::getline(states_file, line)) {
        if (line[0] == '#')
            continue;        
        std::stringstream labor(line);
        labor >> state_num;
        microstates_set[state_num] = new G_microstate;
        while (labor >> lig_type) {
            labor >> lig_num;
            microstates_set[state_num]->holder.push_back( G_ligandLocus(lig_type, lig_num) );
        }        
        number_of_states ++;
    }
    if (!silence)
        std::cout  << number_of_states << " states" << " confirmed\n";
    states_file.close();

    // read in ligand set file
    std::ifstream lig_file;
    std::vector<G_ligand> lig_set;
    double diam_lig;
    double L_c_lig, L_p_lig;
    lig_file.open (dir+"/ligand_set.lig");
    while (std::getline(lig_file, line)) {
        if (line[0] == '#')
            continue;
        G_ligand lig;
        std::stringstream labor(line);
        labor >> diam_lig;
        lig.domain.push_back(G_Domain(diam_lig));
        while (labor >> L_c_lig) {
            labor >> L_p_lig;
            labor >> diam_lig;
            lig.linker.push_back ( G_Linker(L_c_lig, L_p_lig) );
            lig.domain.push_back ( G_Domain(diam_lig) );
        }
        lig_set.push_back(lig);
    }
    if (!silence)
        std::cout << lig_set.size() << " ligands confirmed (receptor included)\n";
    lig_file.close();
    
    // read in cell surface file
    std::vector<double> Density_surfaceLig;
    if (CellSurface_mode) {     // cell surface mode on
        std::ifstream CS_file;
        CS_file.open (dir+"/CellSurfaceConsts.csc");
        double frog;    // frog has a huge mouth. they are good at swallowing and puking.
        while (std::getline(CS_file, line)) {
            if (line[0] == '#')
                continue;
            std::stringstream labor(line);
            while (labor >> frog) {
                Density_surfaceLig.push_back (frog / 1e20); // 1e20 is the conversion factor that converts the number of ligands per square meter to which per square angstrom
            }
            break;
        }
        if (!silence)
            std::cout << Density_surfaceLig.size() << " surface ligands confirmed (receptor not included)\n";
        CS_file.close();
    }
    
    squash_set PPD;
    if (PreProcess_switch) {
        if (PreProcess_job == -1) {
            PPD.squash_plant(lig_set);
            PPD.output(dir);
        }
        else {
            std::cout << "preprocessing ligand #" << PreProcess_job << "\n";
            pumpkin spiderman(lig_set[PreProcess_job].domain.size());
            spiderman.fill_pumpkin(lig_set[PreProcess_job]);
            std::cout << "\nsaving preprocessed data:\n";

            std::ofstream file_out (dir + "/temp/Pumpkin_" + std::to_string(PreProcess_job) + ".ppd");
            spiderman.output (file_out);
            file_out.close();
            
            std::cout << "finished\n";
        }        
        return 0;
    }
    if (use_PPD_switch) {
        PPD.input(dir, lig_set.size());
    }
    // read in working states relations and output the effcs    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (job_config == -1) {
        std::ifstream work_in_file (RLTfile_dir+topology+"relations.udot");
        std::ofstream work_out_file (dir+"/effcs.effc");
        int st_1;
        int st_2;
        int working_relation = 0;
        if (!silence)
            std::cout << "working on relation:  ";
        while (std::getline(work_in_file, line)) {
            if (line[0] == '#')
            continue;

            working_relation++;
            if (!silence) {
                std::cout << std::string(std::to_string(working_relation).length(), '\b');
                std::cout << working_relation << std::flush;
            }
            std::ifstream temp_file (dir+"/temp/temp_"+std::to_string(working_relation)+".effc");
            if (temp_file.is_open()) {  // work already done in temp file
                if (std::getline(temp_file, line)) {
                    work_out_file << line << "\n";
                    temp_file.close();
                    if (clean) {
                        remove((dir+"/temp/temp_"+std::to_string(working_relation)+".effc").c_str());
                    }
                    continue;
                }

            }
            std::stringstream labor(line);
            labor >> st_1;
            labor >> st_2;

            if (CellSurface_mode) {
                work_out_file << st_1 << "\t" << st_2 << "\t" << CellSurface_Effc(microstates_set[st_1], microstates_set[st_2], lig_set, Density_surfaceLig) << "\n";
            }
            else {
                NodesGraph my_graph;
                if (my_graph.BuildGraph(microstates_set[st_1], microstates_set[st_2], lig_set) == NULL) // no effc
                    work_out_file << st_1 << "\t" << st_2 << "\t" << -1 << "\n";
                else {
                    if (use_PPD_switch)
                        work_out_file << st_1 << "\t" << st_2 << "\t" << my_graph.effc(&PPD) << "\n";
                    else
                        work_out_file << st_1 << "\t" << st_2 << "\t" << my_graph.effc() << "\n";
                }
            }            
        }
        work_out_file.close();
        work_in_file.close();

        // free memory
        int labor_i = 1;
        while(microstates_set[labor_i] != NULL) {
            delete microstates_set[labor_i];
            labor_i++;
        }
        if (!silence)
            std::cout << "\nfinished\n";
    }
    else {  // work on specific job
        std::ifstream work_in_file (RLTfile_dir+topology+"relations.udot");
        std::ofstream work_out_file (dir+"/temp/temp_"+std::to_string(job_config)+".effc");
        int st_1;
        int st_2;
        int working_relation = 1;
        if (!silence)
            std::cout << "working on relation:  ";
        while (std::getline(work_in_file, line)) {
            if (line[0] == '#')
                continue;
            if (working_relation < job_config) {
                working_relation++;
                continue;
            }
            if (!silence)
                std::cout << working_relation << std::flush;
            
            std::stringstream labor(line);
            labor >> st_1;
            labor >> st_2;
            if (CellSurface_mode) {
                work_out_file << st_1 << "\t" << st_2 << "\t" << CellSurface_Effc(microstates_set[st_1], microstates_set[st_2], lig_set, Density_surfaceLig) << "\n";
            }
            else {
                NodesGraph my_graph;
                if (my_graph.BuildGraph(microstates_set[st_1], microstates_set[st_2], lig_set) == NULL) // no effc
                    work_out_file << st_1 << "\t" << st_2 << "\t" << -1 << "\n";
                else {
                    if (use_PPD_switch)
                        work_out_file << st_1 << "\t" << st_2 << "\t" << my_graph.effc(&PPD) << "\n";
                    else
                        work_out_file << st_1 << "\t" << st_2 << "\t" << my_graph.effc() << "\n";
                }
            }
            break;
        }
        work_out_file.close();
        work_in_file.close();

        // free memory
        int labor_i = 1;
        while(microstates_set[labor_i] != NULL) {
            delete microstates_set[labor_i];
            labor_i++;
        }
        if (!silence)
            std::cout << "\nfinished\n";
    }
}
