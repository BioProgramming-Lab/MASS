
# MASS: Multivalent Antigen Sensing Simulator

Multivalent Antigen Sensing Simulator (**MASS**) is a Dockerized simulation toolkit designed for modeling multivalent binding interactions, both in solution and on surfaces. It facilitates the flexible and customizable design of multivalent protein binders.

## 🚀 Installation

### Using Docker (Recommended)

```sh
# Pull the prebuilt Docker image from DockerHub
docker image pull torreya/mass:latest

# Run the container with an interactive shell
docker run -it --name test_1 torreya/mass:latest
```

### Manual Setup (Alternative)

To build MASS manually, ensure your environment meets the following requirement:

- `gcc` version 14.1.0

> **Note:** Full manual build instructions will be added in a future update.

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

This file defines the parameters for your simulation. Below are the necessary sections based on the model used.

### On-Surface Simulation

Describes interactions between an $N$-valency binder and $M$ antigen types on the cell surface.

**Required fields:**

1. `name` – Name of the simulation.
2. `RLT_path` – Absolute path to the directory containing RLT files. During the execution of one simulation task, MASS will search this directory for corresponding RLT files and create them automatically if they are not found.
3. `model` – Must be `"OnSurface"`.
4. `geometrical_parameters` – Defines the topology of the binder and antigens:
   - `receptor`: Array of $3N-2$ floats \[diam, Lc, Lp, diam, ..., diam\], where diam is the diameter of the binding domain, Lc and Lp are the contour and persistence length of linkers, respectively. All units are in Å. Here, receptor refers to the binder.
   - `ligand1` to `ligandM`: Diameter of each ligand binding domain in Å.
1. `simulation_settings`:
   - `binder_concentration`: Array of binder concentrations (unit: M).
   - `target_concentration`: Array of $M$ antigen concentrations (unit: M).
   - `target_density`: Array of $M$ antigen surface densities (unit: m⁻²).
   - `output`: List of desired outputs. Options:
     - `final_response`, `final_targets`, `final_states`
     - `response`, `states`, `targets`, `BLI`
   - `association_time`, `dissociation_time` – reaction duration of association and dissociation events in seconds.
   - `saveat` – Time interval for saving outputs (default: 0.02 s).
   - `relative_tolerance`, `absolute_tolerance` – ODE solver precision.
   - `parallel` – Boolean to enable parallel computation. Requires GNU parallel. 
7. `kinetic_config`: Defines binding kinetics for each binder-antigen pair. Format:

```json
{
  "ligand1": [[kon1, koff1], ..., [konN, koffN]],
  ...
}
```

### In-Solution Simulation

This model simulates interactions that happen in solution between $N$-valency targets and $M$ types of binders. This model is modified MVsim. We assume the number of receptor molecules is much less than the binders, meaning a target molecule can bind to multiple binders, while a binder can not bind to multiple receptors.

**Required fields:**

1. `name` - Name of the simulation.
2. `RLT_type` - `"full"` or `"simplified"` (whether (full) or not (simplified), there could be more than two analytes of the same type in a given microstate).
3. `RLT_path` - Absolute path to the directory containing RLT files.
4. `model` – Must be `"InSolution"`.
5. `geometrical_parameters`:
   - `receptor`: Array of $3N-2$ floats \[diam, Lc, Lp, diam, ..., diam\], where diam is the diameter of the binding domain, Lc and Lp are the contour and persistence length of linkers, respectively. All units are in Å. Here, receptor refers to the target molecule.
   - `ligand1` to `ligandM`: Each contains $3V_i - 2$ floats representing multivalent topology, where  $V_i$  is the valency of the binder type i.
1. `simulation_settings`:
   - `binder_concentration`: List of binder concentrations (unit: M).
   - `target_concentration`: Float (unit: M).
   - `output`: `final_response`, `response`
   - `association_time`, `dissociation_time`, `saveat`, `relative_tolerance`, `absolute_tolerance`, `parallel` as described above.
7. `kinetic_config`: Format:

```json
{
  "ligand1": [[kon1_1, koff1_1, ..., kon1_V1, koff1_V1], ..., [konN_1, koffN_1, ..., konN_VN, koffN_VN]],
  ...
}
```

Each `kon_ij, koff_ij` describes the binding between the i-th domain on the receptor and the j-th domain on the ligand.

## 📁 Repository Structure (Coming Soon)

Planned documentation updates will include:

- API references
- Build instructions for manual installations
- Extended tutorials and example analyses
