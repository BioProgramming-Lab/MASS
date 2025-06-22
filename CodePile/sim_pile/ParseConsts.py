import argparse
import numpy as np

# Parse the command line arguments
def ParseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', type=str, default="CodePile/MVsim_test/first_test", help="directory of the simulation")
    parser.add_argument('-t', '--topology', type=str, default="3-3", help="binding topology")
    parser.add_argument('-R', '--RLTfiles', type=str, default="CodePile/full_RLT_files", help="RLT file directory")
    parser.add_argument('-l', '--StepLength', type=float, default=0.02, help="ODE solving save time step length (unit: sec)")
    parser.add_argument('--rtol', type=float, default=1e-3, help="ODE solving relative tolerance")
    parser.add_argument('--atol', type=float, default=1e-10, help="ODE solving absolute tolerance")
    parser.add_argument('-o', '--Output', nargs='+', default=["final_response"], help="output type")
    parser.add_argument('-m', '--Mode', type=str, default="InSolution", help="mode of the simulation, either InSolution or OnSurface")
    parser.add_argument('-j', '--ParallelJob', type=int, default="-1", help="the number of parallel job, -1 for no parallel")
    parser.add_argument('--association_time', type=float, default=180.0, help="association time (unit: sec)")
    parser.add_argument('--dissociation_time', type=float, default=180.0, help="dissociation time (unit: sec)")
    parser.add_argument('--use-previous-sim', action='store_true', help="Use the previously calculated simulation results. Works only when the simulation runs in parallel.")

    args = parser.parse_args()
    # the two arguments are necessary
    if args.directory == "undefined":
        print("Please specify the directory of the simulation by -d")
        exit()
    if args.topology == "undefined":
        print("Please specify the binding topology by -t")
        exit()
    # add "/" to the end of the directory
    if args.directory[-1] != "/":
        args.directory += "/"
    if args.RLTfiles[-1] != "/":
        args.RLTfiles += "/"

    return args

# get the true type of a microstate
# for the full RLT mode, we'll duplicate those ligands to genarate all of the microstates
def GetTrueType(n: int):
    if n <= NumOfLigands:
        return n
    else:
        return int(((n-1-NumOfLigands) / (ReceptorLength-1)) + 1)


# read the file of effective_concentrations
# key->"father,son"     value->effc (-1 for null)
def ReadEffc(dir: str):
    # use dictionary to store effective concentrations, the sparse matrix
    effc_const = {}
    effc_file = open(dir+"effcs.effc", 'r')
    for line in effc_file:
        a = line.split()
        if (a[0] == '#'):
            continue
        effc_const[f"{int(a[0])},{int(a[1])}"] = float(a[2])  # key->"father,son"
    effc_file.close()
    return effc_const

# read number of ligands
def ReadNumberOfLiagnds(dir: str, topology: str):
    num_file = open(dir + topology + "number_of_ligands.txt", 'r')
    num_ligands = int(num_file.readline())
    num_file.close()
    return num_ligands

# read the file of kinetic constants
# key->"Rec_domain;lig_type,lig_domain"     value->[kon, koff]
def ReadKconst(dir: str, number_of_ligands: int):
    K_const = {}
    for lig_num in range(1, number_of_ligands+1):
        Kc_file = open(dir+f"K_const_lig{lig_num}.K", "r")    # read in ligands config
        labor_i = 0
        for line in Kc_file:
            a = line.split() 
            if (a[0] == '#'):
                continue
            labor_i += 1       
            for i in range(1, int(len(a)/2+1)):
                K_const[f"{labor_i};{lig_num},{i}"] = [float(a[i*2-2]), float(a[i*2-1])]  # key->"Rec_domain;lig_type,lig_domain"     value->[kon, koff]
        del labor_i
        Kc_file.close()
    return K_const

# read microstates
def ReadMicrostates(dir: str, topology: str):
    state_file = open(dir+topology+"states_sort.states", "r")     # can only deal with sorted file------------------------------------
    states_num = 0
    states = [[]]   # sub->microstate_num
    for line in state_file:
        a = line.split()
        if (a[0] == '#'):
            continue
        states.append([])
        states_num += 1
        for j in a[1:]:
            states[int(a[0])].append(int(j))
    state_file.close()
    return states

# generate response of each microstate
# # response[i] = number of ligands bound to the receptor in microstate i
def StateResponse(microstates, receptor_length, num_of_ligands):
    response = []
    for i in range(len(microstates)):
        # use book array to record the binding of ligands, making each multivalent bound molecule of 1 unit response rather than multiple units
        book = np.zeros(receptor_length * num_of_ligands + 1)
        for j in range(len(microstates[i])):
            if j % 2 == 0 and microstates[i][j] != -1:
                book[microstates[i][j]] = 1
        response.append( sum(book) )
    return response

# generate consumption of each microstate
def StateConsumption(microstates, num_of_ligands):
    consumption = np.zeros((len(microstates), num_of_ligands + 1)) # this +1 is a result of the index of a list in python starting from 0
    for i in range(len(microstates)):
        for j in range(len(microstates[i])):
            if j % 2 == 0 and microstates[i][j] != -1:
                consumption[i, GetTrueType(microstates[i][j])] += 1
    return consumption

# read the relations of microstates
def ReadRelations(dir: str, topology:str):
    son_relation = [[]]
    father_relation = [[]]
    relation_file = open(dir+topology+"relations_sort.rlt", "r")
    for line in relation_file:
        a = line.split()
        if (a[0] == '#'):
            continue
        # one more microstate
        son_relation.append([])
        father_relation.append([])

        mode = 'N'  # N -> null, F -> father, S -> son
        for i in a:
            if i == 'F':
                mode = 'F'
                continue
            elif i == 'S':
                mode = 'S'
                continue

            if mode == 'N':
                continue
            if mode == 'F':
                father_relation[int(a[0])].append(int(i))
            elif mode == 'S':
                son_relation[int(a[0])].append(int(i))
    relation_file.close()
    return son_relation, father_relation

# read the configuration of concentrations
def ReadConcConfig(dir: str):
    conc_file = open(dir+"conc.conc", "r")
    conc_config = []
    for line in conc_file:
        a = line.split()
        if (a[0] == '#'):
            continue
        conc_config.append(float(a[0]))
    conc_file.close()
    return conc_config

# get the configuration of C
def get_C (m_1, m_2):    # father->m_1  son->m_2
    for i in range(0, len(Microstates[m_1]), 2):
        if (Microstates[m_1][i] == -1 and Microstates[m_2][i] != -1):
            return ConcConfig[GetTrueType(Microstates[m_2][i])]
    print("loading error2")
    quit()

def get_k_on (m_1, m_2):    # father->m_1  son->m_2
    for i in range(0, len(Microstates[m_1]), 2):
        if (Microstates[m_1][i] == -1 and Microstates[m_2][i] != -1):
            return KConst[f"{int(i/2)+1};{GetTrueType(Microstates[m_2][i])},{Microstates[m_2][i+1]}"][0]
    print("loading error2")
    quit()

def get_k_off (m_1, m_2):    # father->m_1  son->m_2
    for i in range(0, len(Microstates[m_1]), 2):
        if (Microstates[m_1][i] == -1 and Microstates[m_2][i] != -1):
            return KConst[f"{int(i/2)+1};{GetTrueType(Microstates[m_2][i])},{Microstates[m_2][i+1]}"][1]
    print("loading error2")
    quit()

def get_binding_ligand_type (m_1, m_2):    # father->m_1  son->m_2
    for i in range(0, len(Microstates[m_1]), 2):
        if (Microstates[m_1][i] == -1 and Microstates[m_2][i] != -1):
            return GetTrueType(Microstates[m_2][i])

# rearrange the data
def RearrangeData():
    # fill in k_on and k_off and C_dic dic
    k_on = {}
    k_off = {}
    C_dic = {}
    BindingLigand = {}
    for i in range(1, len(SonRelation)):     # traverse all microstates
        for j in SonRelation[i]:
            k_on[f"{i},{j}"] = get_k_on(i, j)
            k_on[f"{j},{i}"] = get_k_on(i, j)
            k_off[f"{i},{j}"] = get_k_off(i, j)
            k_off[f"{j},{i}"] = get_k_off(i, j)
            C_dic[f"{i},{j}"] = get_C(i, j)
            C_dic[f"{j},{i}"] = get_C(i, j)
            BindingLigand[f"{i},{j}"] = get_binding_ligand_type(i, j)
            BindingLigand[f"{j},{i}"] = get_binding_ligand_type(i, j)
    return k_on, k_off, C_dic, BindingLigand

# reload the configuration of concentrations
def ReloadConcConfig(ConcSeq):
    global ConcConfig
    global ReceptorConc
    global LigandConc
    ConcConfig = ConcSeq
    ReceptorConc = ConcConfig[0]
    C_dic = {}
    for i in range(1, len(SonRelation)):     # traverse all microstates
        for j in SonRelation[i]:
            C_dic[f"{i},{j}"] = get_C(i, j)
            C_dic[f"{j},{i}"] = get_C(i, j)
    LigandConc = C_dic

# set the concentration of titrate ligand (for in solution simulation)
# currently only works for the case NumOfLigands == 1 ------------------------------------------------------------------------------------------
def SetTitrateConc(Conc):
    TitrateConc = []
    TitrateConc.append(ReceptorConc)
    for i in range(NumOfLigands):
        TitrateConc.append(Conc)
    ReloadConcConfig(TitrateConc)

def SetDissociationConc():
    DissociationConc = []
    DissociationConc.append(ReceptorConc)
    for i in range(NumOfLigands):
        DissociationConc.append(0)
    ReloadConcConfig(DissociationConc)

Args = ParseCommandLine()
NumOfLigands = ReadNumberOfLiagnds(Args.RLTfiles, Args.topology)
Microstates = ReadMicrostates(Args.RLTfiles, Args.topology)
NumOfMicrostates = len(Microstates)-1
ReceptorLength = len(Microstates[1]) // 2
ResponseOfState = StateResponse(Microstates, ReceptorLength, NumOfLigands)
ConsumptionOfState = StateConsumption(Microstates, NumOfLigands)

EffcConst = ReadEffc(Args.directory)
KConst = ReadKconst(Args.directory, NumOfLigands)
SonRelation, FatherRelation = ReadRelations(Args.RLTfiles, Args.topology)
ConcConfig = ReadConcConfig(Args.directory)
ReceptorConc = ConcConfig[0]

Kon, Koff, LigandConc, BindingLigand = RearrangeData()
# Kon: ["father,son"] -> kon, ["father,son"] == ["son,father"]
# Koff: ["father,son"] -> koff, ["father,son"] == ["son,father"]
# LigandConc: ["father,son"] -> conc, ["father,son"] == ["son,father"]
# BindingLigand: ["father,son"] -> ligand type, ["father,son"] == ["son,father"]