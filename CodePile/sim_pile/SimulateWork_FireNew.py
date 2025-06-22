#!/usr/bin/python
from scipy.integrate import solve_ivp
import numpy as np
import ParseConsts
import ODEfuncs
import matplotlib.pyplot as plt
import json
from json import JSONEncoder
import os

# subclass JSONEncoder
class ConfigEncoder(JSONEncoder):
        def default(self, o):
            return o.__dict__

# Time [timeStep] -> time value at the given time step
# Response [timeStep] -> Response value at the given time step
class ResponseData:
    def __init__ (self, t_list, y_dict):
        self.Time = t_list
        self.Output = y_dict

# Time [timeStep] -> time value at the given time step
# Solution [timeStep] [microstate_i] -> Concentration of microstate_i at the given time step
class ODEsolution:
    def __init__ (self, t_list, list_of_y_lists):
        self.Time = t_list
        self.Solution = list_of_y_lists
    def __getitem__ (self, index):
        return ODEsolution([self.Time[index]], [self.Solution[index]])

# DataPackage.data -> dictionary of ResponseData
class DataPackage:
    def __init__ (self):
        self.data = {}
    def DumpToJSON (self, filename):
        with open(ParseConsts.Args.directory + filename + ".json", "w") as f:
            json.dump(self.data, f, cls=ConfigEncoder)

def InSolution_FixTime_Sim (var0, StartTimePoint, EndTimePoint): # solve ODE in a fixed time range
    t_list = np.arange(StartTimePoint, EndTimePoint, ParseConsts.Args.StepLength)
    solution = solve_ivp(ODEfuncs.ivp_simu_func, [StartTimePoint, EndTimePoint], var0, method="Radau", t_eval=t_list, rtol=ParseConsts.Args.rtol, atol=ParseConsts.Args.atol)
    return ODEsolution (solution.t, solution.y.T)

# On Surface model simulation
def OnSurface_FixTime_Sim(var0, StartTimePoint, EndTimePoint):
    t_list = np.arange(StartTimePoint, EndTimePoint, ParseConsts.Args.StepLength)
    solution = solve_ivp(ODEfuncs.cell_surf_simu_func, [StartTimePoint, EndTimePoint], var0, method="Radau", t_eval=t_list, rtol=ParseConsts.Args.rtol, atol=ParseConsts.Args.atol)
    return ODEsolution(solution.t, solution.y.T)

# On Surface model simulation sepciailized for dissociation
def OnSurface_FixTimeDis_Sim(var0, StartTimePoint, EndTimePoint):
    t_list = np.arange(StartTimePoint, EndTimePoint, ParseConsts.Args.StepLength)
    solution = solve_ivp(ODEfuncs.cell_surf_simu_dis_func, [StartTimePoint, EndTimePoint], var0, method="Radau", t_eval=t_list, rtol=ParseConsts.Args.rtol, atol=ParseConsts.Args.atol)
    return ODEsolution(solution.t, solution.y.T)

# get the time line of a series of ODEsolutions
def GetTimeLine (list_of_ODEsolutions):
    TimeLength = 0
    for i in range(len(list_of_ODEsolutions)):
        TimeLength += len(list_of_ODEsolutions[i].Time)
    Time = np.zeros(TimeLength - len(list_of_ODEsolutions) + 1)
    current_time = 0
    for i in list_of_ODEsolutions:
        for j in range(len(i.Time)):
            Time[current_time] = i.Time[j]
            current_time += 1
        current_time -= 1
    return Time.tolist()

# response = sum of all microstates
def OnSurfaceSimToResponse(the_list_of_ODEsolutions, final_only = False):
    if final_only:  # only the last time point is needed
        list_of_ODEsolutions = [the_list_of_ODEsolutions[-1][-1]]
    else:
        list_of_ODEsolutions = the_list_of_ODEsolutions

    TimeLength = 0
    for i in range(len(list_of_ODEsolutions)):
        TimeLength += len(list_of_ODEsolutions[i].Time)
    Response = np.zeros(TimeLength - len(list_of_ODEsolutions) + 1)
    current_time = 0
    for i in list_of_ODEsolutions:
        for j in range(len(i.Time)):
            Response[current_time] = sum(i.Solution[j][1:ParseConsts.NumOfMicrostates])
            current_time += 1
        current_time -= 1
    return Response.tolist()

# response = sum (1/#bound_antigens)
# in bivalent binding BLI experiment
def BLISimToResponse(the_list_of_ODEsolutions, final_only = False):
    if final_only:  # only the last time point is needed
        list_of_ODEsolutions = [the_list_of_ODEsolutions[-1][-1]]
    else:
        list_of_ODEsolutions = the_list_of_ODEsolutions
    TimeLength = 0
    for i in range(len(list_of_ODEsolutions)):
        TimeLength += len(list_of_ODEsolutions[i].Time)
    Response = np.zeros(TimeLength - len(list_of_ODEsolutions) + 1)
    current_time = 0
    for i in list_of_ODEsolutions:
        for j in range(len(i.Time)):
            Response[current_time] = sum([i.Solution[j][k]/ParseConsts.ResponseOfState[k+1] for k in range(1, ParseConsts.NumOfMicrostates)])
            current_time += 1
        current_time -= 1
    return Response.tolist()

# a plugin for analysis of multivalent binding behavior
# return the each states of the last time point
def AllStates_OnSurfaceSimToResponse(list_of_ODEsolutions, final_only = False):
    if final_only:  # only the last time point is needed
        list_of_ODEsolutions = [list_of_ODEsolutions[-1][-1]]
    else:
        list_of_ODEsolutions = list_of_ODEsolutions
    
    Response = []
    current_time = 0
    for i in list_of_ODEsolutions:
        for j in range(len(i.Time)):
            if len(Response) > current_time:    # the current time point is already in the list, which means the time point is shared by two ODEsolutions
                Response[current_time] = i.Solution[j][1:ParseConsts.NumOfMicrostates].tolist()
            else:   # the current time point is not in the list
                Response.append(i.Solution[j][1:ParseConsts.NumOfMicrostates].tolist())
            current_time += 1
        current_time -= 1   # remove the last time point, which is the same as the first time point of the next ODEsolution
    return Response

# a plugin for analysis of multivalent binding behavior
# return concentrations of the each targets of the last time point
def Targets_OnSurfaceSimToResponse(list_of_ODEsolutions, final_only = False):
    if final_only:  # only the last time point is needed
        list_of_ODEsolutions = [list_of_ODEsolutions[-1][-1]]
    else:
        list_of_ODEsolutions = list_of_ODEsolutions
    
    Response = []
    current_time = 0
    for i in list_of_ODEsolutions:
        for j in range(len(i.Time)):
            if len(Response) > current_time:    # the current time point is already in the list, which means the time point is shared by two ODEsolutions
                Response[current_time] = i.Solution[j][ParseConsts.NumOfMicrostates:].tolist()
            else:   # the current time point is not in the list
                Response.append(i.Solution[j][ParseConsts.NumOfMicrostates:].tolist())
            current_time += 1
        current_time -= 1   # remove the last time point, which is the same as the first time point of the next ODEsolution
    return Response

def InSolutionSimToResponse (list_of_ODEsolutions, final_only = False):
    if final_only:  # only the last time point is needed
        list_of_ODEsolutions = [list_of_ODEsolutions[-1][-1]]
    else:
        list_of_ODEsolutions = list_of_ODEsolutions
    
    TimeLength = 0
    for i in range(len(list_of_ODEsolutions)):
        TimeLength += len(list_of_ODEsolutions[i].Time)
    Response = np.zeros(TimeLength - len(list_of_ODEsolutions) + 1)
    current_time = 0
    for i in list_of_ODEsolutions:
        for j in range(len(i.Time)):
            Response[current_time] = sum([ParseConsts.ResponseOfState[k]*i.Solution[j][k] for k in range(len(i.Solution[j]))])
            current_time += 1
        current_time -= 1 # remove the last time point, which is the same as the first time point of the next ODEsolution
    return Response.tolist()

# to integrate all required response output and put them into a dictionary
def OverAll_SimToResponse(list_of_ODEsolutions, response_func: dict, final_only: dict):
    all_final = True
    for key in final_only:
        all_final = all_final and final_only[key]

    output = {}
    if all_final:   # only the last time point is needed
        output['Time'] = GetTimeLine([list_of_ODEsolutions[-1][-1]])
    else:
        output['Time'] = GetTimeLine(list_of_ODEsolutions)
    output['Output'] = {}
    for key in response_func:
        output['Output'][key] = response_func[key](list_of_ODEsolutions, final_only[key])
    return output

# the function for external call
def OnSurfaceSim(association_time: float, dissociation_time: float, guest_concentration: float, final_only = {}, response_func = {}):
    var0 = np.zeros(ParseConsts.NumOfMicrostates + ParseConsts.NumOfLigands)
    var0[0] = guest_concentration
    for i in range(ParseConsts.NumOfLigands):
        var0[ParseConsts.NumOfMicrostates + i] = ParseConsts.ConcConfig[i+1]
    sov1 = OnSurface_FixTime_Sim(var0, 0, association_time)
    
    global TheDataPackage
    if (dissociation_time >= ParseConsts.Args.StepLength):
        sov2_start = sov1.Solution[-1]
        sov2_start[0] = 0.0
        sov2 = OnSurface_FixTimeDis_Sim(sov2_start, association_time, association_time + dissociation_time)
        TheDataPackage.data[guest_concentration] = OverAll_SimToResponse([sov1, sov2], response_func, final_only)
    else:
        TheDataPackage.data[guest_concentration] = OverAll_SimToResponse([sov1], response_func, final_only)

def InSolutionSim (association_time: float, dissociation_time: float, guest_concentration: float, final_only = {}, response_func = {}):
    ParseConsts.SetTitrateConc(guest_concentration)
    ODEfuncs.ReloadOnConcMerge()
    var0 = np.zeros(ParseConsts.NumOfMicrostates)
    var0[0] = ParseConsts.ReceptorConc
    sov1 = InSolution_FixTime_Sim(var0, 0, association_time)
    
    global TheDataPackage
    if (dissociation_time >= ParseConsts.Args.StepLength):
        sov2_start = sov1.Solution[-1]
        ParseConsts.SetDissociationConc()
        ODEfuncs.ReloadOnConcMerge()
        sov2 = InSolution_FixTime_Sim(sov2_start, association_time, association_time + dissociation_time)
        TheDataPackage.data[guest_concentration] = OverAll_SimToResponse([sov1, sov2], response_func, final_only)
    else:
        TheDataPackage.data[guest_concentration] = OverAll_SimToResponse([sov1], response_func, final_only)

TheDataPackage = DataPackage()