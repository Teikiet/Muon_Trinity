#!/usr/bin/python3

import argparse
import os

import matplotlib as mpl
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

# Load the primary particle information from the shower
lib_prim = lib.get("primary")
if lib_prim is None:
    raise RuntimeError("Could not get primary lib from " + str(args.input_dir))
primaries = lib_prim.data
primary = primaries[0]
primary_config = lib_prim.config

# Get the contents of the "particles" sub-directory
lib_particles = lib.get("particles")
if lib_particles is None:
    raise RuntimeError("Could not get particle lib from " + str(args.input_dir))
particles_config = lib_particles.config  # meta information
particles = lib_particles.astype("pandas")  # particle info

# Quick sanity check
if particles_config["plane"] != primary_config["plane"]:
    print("Primary", primary_config["plane"])
    print("Particles on observation plane", particles_config["plane"])
    msg = "The observation plane of the primary is not the same as that"
    msg += " of the particle output. This example will not work as intended"
    msg += " since they do not share the same coordinate system"
    raise RuntimeError(msg)

r = np.sqrt(particles["x"] ** 2 + particles["y"] ** 2)

if not len(r):
    print("No particles were found on the observation plane, quitting...")
    exit()

max_r = max(r)

######################
# Make a plot of the x-y positions of all particles on the observation plane
# Note that these are in the coordinate system of the plane
######################

n_bins = int(np.sqrt(len(r)) * 2)
bins = np.linspace(-max_r, max_r, n_bins)

length_unit_str = particles_config["units"]["length"]  # read in the printed units

title = f"Primary: {primary.name}," + r" E$_{\rm tot}$:"
title += f" {primary.total_energy:.2e} {primary_config['units']['energy']},"
title += f" Zen: {np.rad2deg(np.pi - np.arccos(primary.nz)):0.1f} deg"

fig, ax = plt.subplots(1, 1)
ax.set_aspect("equal")
ax.set_title(particles_config["type"] + "\n" + title)
ax.set_xlabel(f"observation plane x ({length_unit_str})")
ax.set_ylabel(f"observation plane y ({length_unit_str})")

ax.hist2d(
    particles["x"],
    particles["y"],
    weights=particles["weight"],
    bins=(bins, bins),
    norm=mpl.colors.LogNorm(),
)

plot_path = os.path.join(args.output_dir, "particle_dist.png")
print("Saving plot", plot_path)
fig.savefig(plot_path)

######################
# Make a plot of the lateral distribution of particles
# this is based off the approximate shower axis
######################

r_dot_n = (particles["x"] - primary.x) * primary.nx + (
    particles["y"] - primary.y
) * primary.ny
displacement_sq = (particles["x"] - primary.x) ** 2 + (particles["y"] - primary.y) ** 2

if "z" in particles.keys():  # z may not be written out
    r_dot_n += (particles["z"] - primary.z) * primary.nz
    displacement_sq += (particles["z"] - primary.z) ** 2

r_axis = np.sqrt(displacement_sq - r_dot_n**2)

bins = np.linspace(0, max(r_axis) * 1.05, int(n_bins / 2 + 1))

fig, ax = plt.subplots(1, 1)
ax.set_title(title)
ax.set_xlabel(f"distance to shower axis ({length_unit_str})")
ax.set_ylabel("number of particles")
ax.set_yscale("log")

pids = [22, 11, -11, 211, -211, 13, -13, 2212, 2112]  # List of common EAS particles

for pid in pids:
    name = particle.Particle.from_pdgid(pid).latex_name
    if pid not in particles["pdg"]:
        continue
    plt.hist(
        r_axis,
        weights=(particles["pdg"] == pid).astype(float),
        bins=bins,
        label=f"${name}$, {pid}",
        histtype="step",
    )

ax.legend()

plot_path = os.path.join(args.output_dir, "particle_lateral.png")
print("Saving plot", plot_path)
fig.savefig(plot_path)
