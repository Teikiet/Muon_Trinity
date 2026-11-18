#!/usr/bin/python3

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

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


def plot_avg_profile(dat, part, ax):
    # Plots the average profile for particle type `part`
    # averages over all showers in the simulation output dir
    nshower = len(dat.shower.unique())  # number of showers in output dir
    max_dX = np.max(dat["X"])
    h = np.histogram(
        dat.X,
        bins=np.linspace(0, max_dX, len(dat["X"]) - 1),
        weights=dat[part] * 1 / float(nshower),
    )
    ax.plot(h[1][:-1], h[0], label=part)


# Load the shower simulation output
lib = Library(args.input_dir)

# Load the primary particle information from the shower
primaries = lib.get("primary").data
primary = primaries[0]
primary_config = lib.get("primary").config

# Get the contents of the "profile" sub-directory
pr_config = lib.get("profile").config  # meta information
pr = lib.get("profile").astype("pandas")  # particle-number values

# Get the contents of the "production" sub-directory
pr_prod_config = lib.get("production_profile").config  # meta information
pr_prod = lib.get("production_profile").astype("pandas")  # particle-number values
pr_prod_max = lib.get("production_profile").xmumax

# Get the contents of the "energyloss" sub-directory
prdx_config = lib.get("energyloss").config  # meta information
prdx = lib.get("energyloss").astype("pandas")  # energy loss spectrum
prdx_max = lib.get("energyloss").xmax

title = f"Primary: {primary.name}," + r" E$_{\rm tot}$:"
title += f" {primary.total_energy:.2e} {primary_config['units']['energy']},"
title += f" Zen: {np.rad2deg(np.pi - np.arccos(primary.nz)):0.1f} deg"
title += f" Xmax: {prdx_max} g/cm²"


def draw_profiles(pr, pr_config):
    # Plots the number of particles as a function of slant depth
    # individual distributions are made for each particle type
    fig, ax = plt.subplots(1, 1)
    for part in pr.keys():
        # Skip shower id number and slant depth columns
        if "shower" == part or "X" == part:
            continue
        plot_avg_profile(pr, part, ax)
    ax.set_title(pr_config["type"] + "\n" + title)
    unit_grammage_str = pr_config["units"]["grammage"]  # get units of simulation output
    ax.set_xlabel(f"slant depth, X ({unit_grammage_str})")
    ax.set_ylabel("N(X)")
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlim(min(pr["X"]), max(pr["X"]))

    plot_path = os.path.join(args.output_dir, "shower_profile_long_profiles.png")
    print("Saving", plot_path)
    fig.savefig(plot_path)


def draw_production_profiles(pr, pr_config, pr_max):
    # Plots the number of muons that are produced as a function of slant depth
    # individual distributions are made for each particle type
    fig, ax = plt.subplots(1, 1)
    for part in pr.keys():
        # Skip shower id number and slant depth columns
        if "shower" == part or "X" == part:
            continue
        plot_avg_profile(pr, part, ax)
    ax.set_title(
        pr_config["type"]
        + "\n"
        + title
        + r"$\langle X_{\mu}^{\rm max}\rangle=%6.2f$" % np.mean(pr_max)
    )
    unit_grammage_str = pr_config["units"]["grammage"]  # get units of simulation output
    ax.set_xlabel(f"slant depth, X ({unit_grammage_str})")
    ax.set_ylabel("N(X)")
    ax.legend()
    ax.set_yscale("log")
    sel = pr["all"] > 0
    ax.set_xlim(min(pr["X"]), max(pr["X"][sel]))

    plot_path = os.path.join(args.output_dir, "shower_profile_prod_profiles.png")
    print("Saving", plot_path)
    fig.savefig(plot_path)


def draw_energyloss(prdx, prdx_config):
    # Plots the energy deposition as a function of slant depth
    fig, ax = plt.subplots(1, 1)
    ax.set_title(prdx_config["type"] + "\n" + title)
    plot_avg_profile(prdx, "total", ax)
    unit_energy_str = prdx_config["units"]["energy"]  # get units of simulation output
    unit_grammage_str = prdx_config["units"]["grammage"]
    ax.set_xlabel(f"slant depth, X ({unit_grammage_str})")
    ax.set_ylabel(f"dE/dX ({unit_energy_str} / ({unit_grammage_str}))")
    ax.set_xlim(min(prdx["X"]), max(prdx["X"]))

    plot_path = os.path.join(args.output_dir, "shower_profile_energy_loss.png")
    print("Saving", plot_path)
    fig.savefig(plot_path)


draw_profiles(pr, pr_config)
draw_production_profiles(pr_prod, pr_prod_config, pr_prod_max)
draw_energyloss(prdx, prdx_config)
