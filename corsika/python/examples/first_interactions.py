#!/usr/bin/python3

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import particle

from corsika8.io import Library

here = os.path.abspath(os.path.dirname(__file__))

parser = argparse.ArgumentParser()
parser.add_argument(
    "--input-dir", required=True, help="output directory of the CORSIKA 8 simulation"
)
parser.add_argument(
    "--output-dir",
    default=os.path.join(here, "example_plots"),
    help="output directory for plots",
)
args = parser.parse_args()

if not os.path.isdir(args.output_dir):
    print("Making directory", args.output_dir)
    os.makedirs(args.output_dir)

# Load the shower simulation output
lib = Library(args.input_dir)
lib_interactions = lib.get("interactions")
if lib_interactions is None:
    raise RuntimeError("Error getting interactions lib from " + str(args.input_dir))

# Load the primary particle information from the shower
interactions = lib_interactions.data
interaction_config = lib_interactions.config
projectiles = lib_interactions.projectiles  # type: ignore

shower_ids = np.unique(interactions["shower"])

ncols = 2
nrows = len(shower_ids)
fig, ax = plt.subplots(
    ncols=ncols,
    nrows=nrows,
    figsize=(ncols * 8, nrows * 8),
    gridspec_kw={"wspace": 0.2, "hspace": 0.3},
)
ax = np.atleast_2d(ax)
if ax.shape == (1, 2):
    ax = ax.T

for ish, sh_id in enumerate(shower_ids):
    # Find all of the secondary particles associated with this shower
    daughters = interactions[interactions["shower"] == sh_id]

    mother = projectiles[ish]

    print("Plotting shower", sh_id)
    print(mother)

    px, py, pz = mother.momentum

    mother_name = f"${particle.Particle.from_pdgid(mother.pdg).latex_name}$"
    ax[0, sh_id].plot([-px, 0], [-pz, 0], label=mother_name)
    ax[1, sh_id].plot([-py, 0], [-pz, 0], label=mother_name)

    # Get daughters, filter out zero-momentum particles
    p_tot = np.sqrt(daughters["px"] ** 2 + daughters["py"] ** 2 + daughters["pz"] ** 2)
    daughters = daughters[p_tot > 0]

    # keep track of largest values
    max_pz = abs(pz)
    max_pxy = max(abs(px), abs(py))

    # plot the daughter particles
    for i in range(len(daughters)):
        pz = daughters["pz"].iloc[i]
        px = daughters["px"].iloc[i]
        py = daughters["py"].iloc[i]

        pdg = int(daughters["pdg"].iloc[i])
        try:
            name = f"${particle.Particle.from_pdgid(pdg).latex_name}$"
        except Exception:
            name = "Unknown"

        ax[0, sh_id].plot([px, 0], [pz, 0], label=name)
        ax[1, sh_id].plot([py, 0], [pz, 0], label=name)

        max_pz = max(max_pz, abs(pz))
        max_pxy = max(max_pxy, max(abs(px), abs(py)))

    # set plotting style and names
    for i in range(2):
        ax[i, sh_id].legend()

        ymin, ymax = ax[i, sh_id].get_ylim()
        ax[i, sh_id].set_ylim(-max_pz, max_pz)
        ax[i, sh_id].set_xlim(-max_pxy, max_pxy)

        unit_energy_str = interaction_config["units"]["energy"]
        title = f"Mother: {mother.kinetic_energy:0.2e} {unit_energy_str}"
        title += ", " + mother_name
        title += f"\n{mother.n_secondaries} secondaries"
        ax[i, sh_id].set_title(title)
        ax[i, sh_id].set_ylabel(f"$P_z$ ({unit_energy_str})")

    ax[0, sh_id].set_xlabel(f"$P_x$ ({unit_energy_str})")
    ax[1, sh_id].set_xlabel(f"$P_y$ ({unit_energy_str})")

plot_path = os.path.join(args.output_dir, "interactions.png")
print("Saving", plot_path)
fig.savefig(plot_path, bbox_inches="tight")
