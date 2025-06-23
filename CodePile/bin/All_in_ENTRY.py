#!/usr/bin/python
# All in one entry

import sys
import os
import json
import subprocess
from colorama import init, Fore, Back, Style
import argparse

GRLT_path = "GRLT"
Sort_file_path = "/usr/src/CodePile/bin/sort_file.py"
Geffc_path = "Geffc"
Sim_path = "/usr/src/CodePile/sim_pile/titrate.py"
Merge_path = "/usr/src/CodePile/sim_pile/merge_results.py"

def ParseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--Directory', type=str, default="Null", help="directory of the simulation")
    parser.add_argument('--use-previous-effc', action='store_true', help="Use the previously calculated Effc file")
    parser.add_argument('--use-previous-sim', action='store_true', help="Use the previously calculated simulation results. Works only when the simulation runs in parallel.")

    args = parser.parse_args()
    # the two arguments are necessary
    if args.Directory == "Null":
        print("Please specify the directory of the simulation by -d")
        exit()
    # make the directory absolute path
    if args.Directory[0] != "/":
        args.Directory = os.getcwd() + "/" + args.Directory
    # add "/" to the end of the directory
    if args.Directory[-1] != "/":
        args.Directory += "/"

    return args

def ConfigFileCheck(config_data):
    if (not "geometrical_parameters" in config_data.keys()):    # check if "geometrical_parameters" exists
        print ("invalid config file! Missing term \"geometrical_parameters\"\n")
        exit()    
    if (not "model" in config_data.keys()):    # check if "model" exists
        print ("invalid config file! Missing term \"model\"\n")
        exit()
    if (config_data["model"] != "InSolution" and config_data["model"] != "OnSurface"):    # check if "model" is valid
        print ("invalid model. Model must be either \"InSolution\" or \"OnSurface\"\n")
        exit()
    if (not "RLT_type" in config_data.keys()):    # check if "RLT_type" exists
        print ("invalid config file! Missing term \"RLT_type\"\n")
        exit()
    if (config_data["RLT_type"] != "full" and config_data["RLT_type"] != "simplified"):    # check if "RLT_type" is valid
        print ("invalid RLT type. RLT type must be either \"full\" or \"simplified\"\n")
        exit()
    if (not "RLT_path" in config_data.keys()):    # check if "RLT_path" exists
        print ("invalid config file! Missing term \"RLT_path\"\n")
        exit()
    if (not "simulation_settings" in config_data.keys()):    # check if "simulation_settings" exists
        print ("invalid config file! Missing term \"simulation_settings\"\n")
        exit()
    else:
        if (not "binder_concentration" in config_data["simulation_settings"].keys()):
            print ("invalid config file! Missing term \"simulation_settings\"/\"binder_concentration\"\n")
            exit()
        if (not "target_concentration" in config_data["simulation_settings"].keys()):
            print ("invalid config file! Missing term \"simulation_settings\"/\"target_concentration\"\n")
            exit()
        if (not "association_time" in config_data["simulation_settings"].keys()):
            print ("invalid config file! Missing term \"simulation_settings\"/\"association_time\"\n")
            exit()
        if (not "dissociation_time" in config_data["simulation_settings"].keys()):
            print ("invalid config file! Missing term \"simulation_settings\"/\"dissociation_time\"\n")
            exit()
        

def ConfigFileDigestion():
    # read in config.json
    all_in_one_config_file = open(Args.Directory + "config.json", 'r')
    config_data = json.load(all_in_one_config_file)
    all_in_one_config_file.close()

    if (config_data["model"] == "OnSurface"):
        config_data["RLT_type"] = "full"  # OnSurface model only supports full RLT type
        config_data["geometrical_parameters"]["receptor"] = config_data["geometrical_parameters"]["binder"]  # receptor is the binder in OnSurface model
        topology_i = 1
        while f"antigen{topology_i}" in config_data["geometrical_parameters"].keys():
            config_data["geometrical_parameters"][f"ligand{topology_i}"] = config_data["geometrical_parameters"][f"antigen{topology_i}"]
            config_data["kinetic_config"][f"ligand{topology_i}"] = config_data["kinetic_config"][f"antigen{topology_i}"]  # OnSurface model uses antigen{topology_i} as ligand{topology_i}
            topology_i = topology_i + 1
        config_data["simulation_settings"]["target_concentration"] = config_data["simulation_settings"]["antigen_concentration"]  # target_concentration is the antigen_concentration in OnSurface model
        config_data["simulation_settings"]["target_density"] = config_data["simulation_settings"]["antigen_density"]  # target_density is the antigen_density in OnSurface model
    if (config_data["model"] == "InSolution"):
        config_data["geometrical_parameters"]["receptor"] = config_data["geometrical_parameters"]["target"]  # receptor is the target in InSolution model
        topology_i = 1
        while f"binder{topology_i}" in config_data["geometrical_parameters"].keys():
            config_data["geometrical_parameters"][f"ligand{topology_i}"] = config_data["geometrical_parameters"][f"binder{topology_i}"] # InSolution model uses binder{topology_i} as ligand{topology_i}
            config_data["kinetic_config"][f"ligand{topology_i}"] = config_data["kinetic_config"][f"binder{topology_i}"]  # InSolution model uses binder{topology_i} as ligand{topology_i}
            topology_i = topology_i + 1

    # check if the config file is valid
    ConfigFileCheck(config_data)

    print ("digesting config.json")

    # read topology
    receptor_length = (len(config_data["geometrical_parameters"]["receptor"]) + 2) / 3  # 3: {contour length, persistence length, radius}
    if receptor_length != int(receptor_length):  # check if receptor_length is an integer
        print ("invalid receptor geometrical parameters\n")
        exit()
    receptor_length = int(receptor_length)
    ligand_length = []
    topology_i = 1
    while f"ligand{topology_i}" in config_data["geometrical_parameters"].keys():
        ligand_length.append((len(config_data["geometrical_parameters"][f"ligand{topology_i}"]) + 2) / 3)
        topology_i = topology_i + 1
        if ligand_length[-1] != int(ligand_length[-1]):    # check if ligand_length is an integer
            print ("invalid ligand geometrical parameters\n")
            exit()
        ligand_length[-1] = int(ligand_length[-1])
    topology = f"{receptor_length}"
    for i in ligand_length:
        topology = topology + f"-{i}"
    config_data["topology"] = topology
    print (f"topology detected: {config_data['topology']}")

    # read in RLT_path
    if config_data["RLT_path"][-1] != "/":
        config_data["RLT_path"] = config_data["RLT_path"] + "/"

    # read in RLT_type
    print (f"RLT files directory: {config_data['RLT_path']}{config_data['RLT_type']}/")

    # write in ligand_set.lig
    ligand_set_file = open(Args.Directory + "ligand_set.lig", 'w')
    for i in config_data["geometrical_parameters"]["receptor"]:
        ligand_set_file.write(f"{i} ")
    for i in range(1, topology_i):
        ligand_set_file.write("\n")
        for j in config_data["geometrical_parameters"][f"ligand{i}"]:
            ligand_set_file.write(f"{j} ")
    ligand_set_file.close()

    # write in K_const_ligX.lig
    for i in range(1, topology_i):
        K_const_ligX_file = open(Args.Directory + f"K_const_lig{i}.K", 'w')
        for k in config_data["kinetic_config"][f"ligand{i}"][0]:
            K_const_ligX_file.write(f"{k} ")
        for j in range(1, len(config_data["kinetic_config"][f"ligand{i}"])):
            K_const_ligX_file.write("\n")
            for k in config_data["kinetic_config"][f"ligand{i}"][j]:
                K_const_ligX_file.write(f"{k} ")
        K_const_ligX_file.close()

    # write in titration_conc.tc and conc.conc
    titration_conc_file = open(Args.Directory + "titration_conc.tc", 'w')
    if (config_data["model"] == "OnSurface"):
        for i in config_data["simulation_settings"]["binder_concentration"]:
            titration_conc_file.write(f"{i} ")
    if (config_data["model"] == "InSolution"):
        for i in config_data["simulation_settings"]["target_concentration"]:
            titration_conc_file.write(f"{i} ")
    titration_conc_file.close()

    conc_file = open(Args.Directory + "conc.conc", 'w')
    if (config_data["model"] == "OnSurface"):
        conc_file.write(f"{config_data['simulation_settings']['binder_concentration'][0]}\n")
        for i in config_data['simulation_settings']['target_concentration']:
            conc_file.write(f"{i}\n")
    else:
        conc_file.write(f"{config_data['simulation_settings']['target_concentration'][0]}\n")
        conc_file.write(f"{config_data['simulation_settings']['binder_concentration'][0]}\n")
    conc_file.close()

    # write in CellSurfaceConsts.csc. This file is only used in OnSurface model
    if config_data["model"] == "OnSurface":
        cellsurf_file = open(Args.Directory + "CellSurfaceConsts.csc", 'w')
        for i in config_data["simulation_settings"]["target_density"]:
            cellsurf_file.write(f"{i} ")
        cellsurf_file.close()

    # read in simulation_settings
    # read in step_length
    if not "saveat" in config_data["simulation_settings"].keys():
        config_data["simulation_settings"]["saveat"] = 0.02
    if not "output" in config_data["simulation_settings"].keys():
        config_data["simulation_settings"]["output"] = ["final_response"]
    if not "relative_tolerance" in config_data["simulation_settings"].keys():
        config_data["simulation_settings"]["relative_tolerance"] = 1e-3
    if not "absolute_tolerance" in config_data["simulation_settings"].keys():
        config_data["simulation_settings"]["absolute_tolerance"] = 1e-10
    
    print ("digestion finished")
    return config_data

# for RLT_files
# if RLT_file not exsits in ../RLT_files, create it
def ensure_RLT_file(RLT_file : str):        # e.g. ensure_RLT_file("2-2")
    try:
        f = open(f"{Config['RLT_path']}{Config['RLT_type']}/{RLT_file}relations_sort.rlt", 'r')
        f.close()
        print (f"RLT file {RLT_file} found\n")
    except IOError:
        print (f"RLT file {RLT_file} not found, creating..." + Style.DIM)
        Thomas_train = RLT_file.split("-")
        RLT_args = ""
        for i in range(0, len(Thomas_train)):
            RLT_args = RLT_args + f" {Thomas_train[i]}"
        # execute GRLT
        os.system(f"{GRLT_path} {Config['RLT_path']}{Config['RLT_type']}/ {Config['RLT_type']} {RLT_args}")
        os.system(f"python {Sort_file_path} {Config['RLT_path']}{Config['RLT_type']}/ {RLT_file}")
        print (Style.RESET_ALL + f"RLT file {RLT_file} created\n")

# effc part
def ensure_effc_file():
    skip_effc = False
    if (Args.use_previous_effc):
        try:
            f = open(f"{Args.Directory}effcs.effc", 'r')
            f.close()
            skip_effc = True
            print ("effective concentrations found...")
            print ("skipping effective concentrations calculation...\n")
        except IOError:
            pass

    if (not skip_effc):
        print ("calculating effective concentrations: \n")
        if Config['model'] == "OnSurface":
            os.system(f"{Geffc_path} -R {Config['RLT_path']}{Config['RLT_type']}/ -d {Args.Directory[0:-1]} -t {Config['topology']} --cs -s")
        else:
            os.system(f"{Geffc_path} -R {Config['RLT_path']}{Config['RLT_type']}/ -d {Args.Directory[0:-1]} -t {Config['topology']} -s")
        # add checking point -------------------------------------------------------------------------
        print (Style.RESET_ALL + "effective concentrations calculated\n")

# sim part
def simulate():
    print ("ODE solving start:\n" + Style.DIM)

    output_type = ""
    for i in Config['simulation_settings']['output']:
        output_type = output_type + f"{i} "
    output_type = output_type.strip()

    if Config['simulation_settings']['parallel']:
        CURLY_BRACE = "{}"
        os.system(f"seq {len(Config['simulation_settings']['binder_concentration'])} | parallel python {Sim_path} -R {Config['RLT_path']}{Config['RLT_type']}/ -d {Args.Directory[0:-1]} -t {Config['topology']} -l {Config['simulation_settings']['saveat']} -o {output_type} -m {Config['model']} --rtol {Config['simulation_settings']['relative_tolerance']} --atol {Config['simulation_settings']['absolute_tolerance']} --association_time {Config['simulation_settings']['association_time']} --dissociation_time {Config['simulation_settings']['dissociation_time']} -j {CURLY_BRACE} --use-previous-sim")
        os.system(f"python {Merge_path} -d {Args.Directory} -n {len(Config['simulation_settings']['binder_concentration'])}")
    else:
        os.system(f"python {Sim_path} -R {Config['RLT_path']}{Config['RLT_type']}/ -d {Args.Directory[0:-1]} -t {Config['topology']} -l {Config['simulation_settings']['saveat']} -o {output_type} -m {Config['model']} --rtol {Config['simulation_settings']['relative_tolerance']} --atol {Config['simulation_settings']['absolute_tolerance']} --association_time {Config['simulation_settings']['association_time']} --dissociation_time {Config['simulation_settings']['dissociation_time']}")
    
    print (Style.RESET_ALL + "\nODE solving finished\n")
    # add checking point -----------------------------------------------------------------------------------------------------------------------------------------------------


init()  # colorama init

Args = ParseCommandLine()
Config = ConfigFileDigestion()

ensure_RLT_file(Config["topology"])
ensure_effc_file()
simulate()

print (Style.RESET_ALL + "All tasks finished\n")