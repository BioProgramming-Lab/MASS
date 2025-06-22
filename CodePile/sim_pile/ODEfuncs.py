import ParseConsts
import numpy as np

def GetOnConc (father, son):
    if ParseConsts.EffcConst[f"{father},{son}"] != -1:
        return ParseConsts.EffcConst[f"{father},{son}"]
    else:
        return ParseConsts.LigandConc[f"{father},{son}"]

def GetMergeConc ():
    merge_conc = {}
    for i in range(1, len(ParseConsts.SonRelation)):
        for j in ParseConsts.SonRelation[i]:
            merge_conc[f"{i},{j}"] = GetOnConc(i, j)
    return merge_conc

def ivp_simu_func (t, vars):
    ddt = np.zeros((len(ParseConsts.SonRelation)-1,))
    for i in range(1, len(ddt)+1):  # traverse all microstates
        for j in ParseConsts.SonRelation[i]:
            ddt[i-1] += vars[j-1] * ParseConsts.Koff[f"{i},{j}"]
            ddt[i-1] -= vars[i-1] * ParseConsts.Kon[f"{i},{j}"] * OnConcMerge[f"{i},{j}"]
        for j in ParseConsts.FatherRelation[i]:
            ddt[i-1] += vars[j-1] * ParseConsts.Kon[f"{i},{j}"] * OnConcMerge[f"{j},{i}"]
            ddt[i-1] -= vars[i-1] * ParseConsts.Koff[f"{i},{j}"]
    #print("var:    ", t, vars)
    #print("ddt:    ", t, ddt)
    return ddt

def cell_surf_simu_func(t, vars):
    ddt = np.zeros((ParseConsts.NumOfMicrostates + ParseConsts.NumOfLigands,))
    for i in range(1, len(ParseConsts.SonRelation)):
        for j in ParseConsts.SonRelation[i]:
            ddt[i - 1] += vars[j - 1] * ParseConsts.Koff[f"{i},{j}"]
            if (ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{i},{j}"]] != 0):
                ddt[i - 1] -= vars[i - 1] * ParseConsts.Kon[f"{i},{j}"] * OnConcMerge[f"{i},{j}"] * vars[ParseConsts.NumOfMicrostates + ParseConsts.BindingLigand[f"{i},{j}"]-1] / ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{i},{j}"]]
        for j in ParseConsts.FatherRelation[i]:
            if (ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{j},{i}"]] != 0):
                ddt[i - 1] += vars[j - 1] * ParseConsts.Kon[f"{i},{j}"] * OnConcMerge[f"{j},{i}"] * vars[ParseConsts.NumOfMicrostates + ParseConsts.BindingLigand[f"{j},{i}"]-1] / ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{j},{i}"]]
            ddt[i - 1] -= vars[i - 1] * ParseConsts.Koff[f"{i},{j}"]
        for ligand_num in range(1, ParseConsts.NumOfLigands + 1):
            ddt[ParseConsts.NumOfMicrostates + ligand_num - 1] -= ddt[i - 1] * ParseConsts.ConsumptionOfState[i][ligand_num]
    return ddt

def cell_surf_simu_dis_func(t, vars):
    ddt = np.zeros((ParseConsts.NumOfMicrostates + ParseConsts.NumOfLigands,))
    for i in range(2, len(ParseConsts.SonRelation)):
        for j in ParseConsts.SonRelation[i]:
            ddt[i - 1] += vars[j - 1] * ParseConsts.Koff[f"{i},{j}"]
            if (ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{i},{j}"]] != 0):
                ddt[i - 1] -= vars[i - 1] * ParseConsts.Kon[f"{i},{j}"] * OnConcMerge[f"{i},{j}"] * vars[ParseConsts.NumOfMicrostates + ParseConsts.BindingLigand[f"{i},{j}"]-1] / ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{i},{j}"]]
        for j in ParseConsts.FatherRelation[i]:
            if (ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{j},{i}"]] != 0):
                ddt[i - 1] += vars[j - 1] * ParseConsts.Kon[f"{i},{j}"] * OnConcMerge[f"{j},{i}"] * vars[ParseConsts.NumOfMicrostates + ParseConsts.BindingLigand[f"{j},{i}"]-1] / ParseConsts.ConcConfig[ParseConsts.BindingLigand[f"{j},{i}"]]
            ddt[i - 1] -= vars[i - 1] * ParseConsts.Koff[f"{i},{j}"]
        for ligand_num in range(1, ParseConsts.NumOfLigands + 1):
            ddt[ParseConsts.NumOfMicrostates + ligand_num - 1] -= ddt[i - 1] * ParseConsts.ConsumptionOfState[i][ligand_num]
    return ddt


OnConcMerge = GetMergeConc()

def ReloadOnConcMerge():
    global OnConcMerge
    OnConcMerge = GetMergeConc()