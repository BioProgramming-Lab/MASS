#!/usr/bin/python
import SimulateWork_FireNew
import ParseConsts

titration_conc_file = open(ParseConsts.Args.directory+"titration_conc.tc", "r")
titration_conc = []
for line in titration_conc_file:
    titration_conc.extend([float(i) for i in line.split()])
titration_conc_file.close()

WorkingConc = titration_conc

Final_only_dict = {"final_response": True, 
                   "final_states": True, 
                   "final_targets": True, 
                   "response": False, 
                   "states": False, 
                   "targets": False,
                   "BLI": False}

OnSurface_Rfunc = {"final_response": SimulateWork_FireNew.OnSurfaceSimToResponse,
                "final_states": SimulateWork_FireNew.AllStates_OnSurfaceSimToResponse,
                "final_targets": SimulateWork_FireNew.Targets_OnSurfaceSimToResponse,
                "response": SimulateWork_FireNew.OnSurfaceSimToResponse,
                "BLI": SimulateWork_FireNew.BLISimToResponse,
                "states": SimulateWork_FireNew.AllStates_OnSurfaceSimToResponse,
                "targets": SimulateWork_FireNew.Targets_OnSurfaceSimToResponse}

InSolution_Rfunc = {"final_response": SimulateWork_FireNew.InSolutionSimToResponse,
                "response": SimulateWork_FireNew.InSolutionSimToResponse}

rfunc_sim = {}
final_sim = {}
for i in ParseConsts.Args.Output:
    if ParseConsts.Args.Mode == "OnSurface":
        rfunc_sim[i] = OnSurface_Rfunc[i]
    else:
        rfunc_sim[i] = InSolution_Rfunc[i]
    final_sim[i] = Final_only_dict[i]

# should do as well to InSolution mode simulation (Rfunc)--------------------------------------------------------

if ParseConsts.Args.ParallelJob == -1:
    for i in WorkingConc:
        if ParseConsts.Args.Mode == "InSolution":
            SimulateWork_FireNew.InSolutionSim(ParseConsts.Args.association_time, ParseConsts.Args.dissociation_time, i, final_only=final_sim, response_func=rfunc_sim)
        elif ParseConsts.Args.Mode == "OnSurface":
            SimulateWork_FireNew.OnSurfaceSim(ParseConsts.Args.association_time, ParseConsts.Args.dissociation_time, i, final_only=final_sim, response_func=rfunc_sim)
    SimulateWork_FireNew.TheDataPackage.DumpToJSON("titration_result")
else:
    # if the job is not the first job, check if the previous job has been done
    if (ParseConsts.Args.use_previous_sim == True):
        try:
            temp_result_file = open(f"{ParseConsts.Args.directory}temp_titration_result_{ParseConsts.Args.ParallelJob}.json", "r")
            temp_result_file.close()
            print (f"Loading previous simulation result for job {ParseConsts.Args.ParallelJob}")
            exit(0)
        # if the previous job has not been done, do the simulation
        except:
            if ParseConsts.Args.Mode == "InSolution":
                SimulateWork_FireNew.InSolutionSim(ParseConsts.Args.association_time, ParseConsts.Args.dissociation_time, WorkingConc[ParseConsts.Args.ParallelJob-1], final_only=final_sim, response_func=rfunc_sim)
            elif ParseConsts.Args.Mode == "OnSurface":
                SimulateWork_FireNew.OnSurfaceSim(ParseConsts.Args.association_time, ParseConsts.Args.dissociation_time, WorkingConc[ParseConsts.Args.ParallelJob-1], final_only=final_sim, response_func=rfunc_sim)
            SimulateWork_FireNew.TheDataPackage.DumpToJSON(f"temp_titration_result_{ParseConsts.Args.ParallelJob}")