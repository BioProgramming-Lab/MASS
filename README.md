
# MASS: Multivalent Antigen Sensing Simulator

Multivalent Antigen Sensing Simulator (**MASS**) is a Dockerized simulation toolkit designed for modeling multivalent interactions. It facilitates flexible and customizable design of multivalent protein binders.
More details can be found in our reasearch paper associated with MASS at: LINK_TO_OUR_PAPER

## 🚀 Installation

### Using Docker (Recommended)

```sh
# Pull the prebuilt Docker image from DockerHub
docker image pull torreya/mass:latest

# Run the container with an interactive shell
docker run -it --name test_1 torreya/mass:latest
```

## 📦 Usage

### Quick Start with Docker

To get started using MASS inside the Docker container:

```sh
MASS -h
```

This will display help information and available options.

### Running a Simulation

1. Prepare a working directory containing your simulation configuration file: `config.json`.
2. Run the simulation using:

```sh
MASS -d RELATIVE_PATH_TO_WORKING_DIR
```

### Example Simulations

Example configuration files for both *in-solution* and *on-surface* modes are available in `/usr/src/tutorial`. To run them:

```sh
MASS -d tutorial/OnSurface
# or
MASS -d tutorial/InSolution
```

## ⚙️ Configuration File (`config.json`)

This file defines the parameters for your simulation. 

### On-Surface Simulation

Describes interactions between an $N$-valency binder and $M$ types of antigens on the cell surface.

**Required fields:**

1. `name` – Name of the simulation.
2. `RLT_path` – Absolute path to the directory containing RLT files, which specify microstate interaction networks. During the execution of one simulation task, MASS will search this directory for corresponding RLT files and create them automatically if they are not found.
3. `model` – Must be `"OnSurface"`.
4. `geometrical_parameters` – Defines the topology of the binder and antigens:
   - `binder`: Array of $3N-2$ floats \[diam, Lc, Lp, diam, ..., diam\], where diam is the diameter of the binding domain, Lc and Lp are the contour and persistence length of linkers, respectively. All units are in Å.
   - `antigen1` to `antigenM`: Diameter of each antigen in Å.
1. `simulation_settings`:
   - `binder_concentration`: Array of binder concentrations (unit: M), MASS will conduct simulation with each of them as if titrating.
   - `antigen_concentration`: Array of $M$ antigen concentrations (unit: M), it has a length of $M$.
   - `antigen_density`: Array of $M$ antigen surface densities (unit: m⁻²), it has a length of $M$.
   - `output`: List of desired outputs. Options:
     - `final_response`, as the concentration of all bound binders at the end of simulation
     - `final_targets`, as the list of concentrations of all unoccupied antigens left at the end of simulation
     - `final_states`, as the list of concentrations of all microstates at the end of simulation
     - `response`, as the concentration of all bound binders at each timepoint during simulation
     - `states`, as the list of concentrations of all microstates at each timepoint during simulation
     - `targets`, as the list of concentrations of all unoccupied antigens left at each timepoint during simulation
     - `BLI`, as the sum($C_i\times 1/n_i$) at each timepoint during simulation, where $C_i$ is the concentration of microstate $i$ of bound binder, $n_i$ is the number of bound antigens in microstate $i$
   - `association_time`, `dissociation_time` – duration of association and dissociation events in seconds.
   - `saveat` – Time interval for saving outputs (default: 0.02 s).
   - `relative_tolerance`, `absolute_tolerance` – ODE solver precision.
   - `parallel` – Boolean to enable parallel computation. Requires GNU parallel. 
7. `kinetic_config`: Defines binding kinetics for each binder-antigen pair. It has a dimension of $M$ rows, with each row containing $N$ pairs of on and off rates. See tutorial for an example. Format:
```
{
  "antigen1": [[kon1, koff1], ..., [konN, koffN]],
  ...
}
```

### In-Solution Simulation

This model simulates interactions that happen in solution between $N$-valency targets and $M$ types of binders. This model is modified from the MVsim toolkit.

**Required fields:**

1. `name` - Name of the simulation.
2. `RLT_type` - `"full"` or `"simplified"`. For the `"full"` type, we assume the number of target molecules is much smaller than that of binders, therefore a target molecule can bind to multiple binders (can be of the same type), while a binder can only bind to one target molecule. For the `"simplified"` type, we simplified the system by assuming that a target molecule can not bind to two binders of the same type. The `"simplified"` type can be used to simulate 1:1 binding, providing similar results while costing less computation resources in most cases.
3. `RLT_path` - Absolute path to the directory containing RLT files, which specify microstate interaction networks. During the execution of one simulation task, MASS will search this directory for corresponding RLT files and create them automatically if they are not found.
4. `model` – Must be `"InSolution"`.
5. `geometrical_parameters`:
   - `target`: Array of $3N-2$ floats \[diam, Lc, Lp, diam, ..., diam\], where diam is the diameter of the binding domain, Lc and Lp are the contour and persistence length of linkers, respectively. All units are in Å.
   - `binder1` to `binderM`: Each contains $3V_i - 2$ floats representing multivalent topology, where  $V_i$  is the valency of the binder type i.
6. `simulation_settings`:
   - `binder_concentration`: Array of binder concentrations (unit: M).
   - `target_concentration`: Array of target concentrations (unit: M), MASS will conduct simulation with each of them as if titrating.
   - `output`: List of desired outputs. Options:
     - `final_response`, as the concentration of all bound binders at the end of simulation
     - `response`, as the concentration of all bound binders at each timepoint during simulation
   - `association_time`, `dissociation_time` – duration of association and dissociation events in seconds.
   - `saveat` – Time interval for saving outputs (default: 0.02 s).
   - `relative_tolerance`, `absolute_tolerance` – ODE solver precision.
   - `parallel` – Boolean to enable parallel computation. Requires GNU parallel. 
7. `kinetic_config`:  Defines binding kinetics for each binder-target pair. It has a dimension of $M$ rows, with each row containing $N$ lists of on and off rates. See tutorial for an example. Format:
```
{
  "binder1": [[kon1_1, koff1_1, ..., kon1_V1, koff1_V1], ..., [konN_1, koffN_1, ..., konN_VN, koffN_VN]],
  ...
}
```
where `kon_ij, koff_ij` describes the binding between the i-th domain on the receptor and the j-th domain on the ligand.
