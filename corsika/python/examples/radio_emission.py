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

# Load the shower simulation output
lib = Library(args.input_dir)

# Load the primary particle information from the shower
lib_prim = lib.get("primary")
if lib_prim is None:
    raise RuntimeError("Could not get primary lib from " + str(args.input_dir))
primaries = lib_prim.data
n_showers = len(primaries)
primary_config = lib_prim.config

found_zhs = False
try:
    zsh_lib = lib.get("ZHS")
    if zsh_lib is None:
        raise RuntimeError("Could not get ZSH lib from " + str(args.input_dir))
    waveforms_config_zhs = zsh_lib.config  # meta information
    waveforms_zhs = zsh_lib.astype("pandas")
    observers_zhs = zsh_lib.get_observers()  # type: ignore
    found_zhs = True
except Exception:
    print("No ZHS waveforms in this directory")


found_coreas = False
try:
    coreas_lib = lib.get("CoREAS")
    if coreas_lib is None:
        raise RuntimeError("Could not get CoREAS lib from " + str(args.input_dir))
    waveforms_config_coreas = coreas_lib.config  # meta information
    waveforms_coreas = coreas_lib.astype("pandas")
    observers_coreas = coreas_lib.get_observers()  # type: ignore
    found_coreas = True
except Exception:
    print("No CoREAS waveforms in this directory")

if not found_zhs and not found_coreas:
    print("No electric field waveforms are in this file, quitting...")
    exit()

for ishower in range(n_showers):

    primary = primaries[ishower]

    title = f"Primary: {primary.name}," + r" E$_{\rm tot}$:"
    title += f" {primary.total_energy:.2e} {primary_config['units']['energy']}"

    # Get all of the unique observer names to match up the CoREAS/ZHS waveforms
    obs_names = []
    if found_zhs:
        obs_names += list(waveforms_zhs[str(ishower)].keys())
        zhs_dict = waveforms_zhs[str(ishower)]
    if found_coreas:
        obs_names += list(waveforms_coreas[str(ishower)].keys())
        coreas_dict = waveforms_coreas[str(ishower)]
    obs_names = np.unique(obs_names)

    ncols = 1
    nrows = len(obs_names)
    fig, ax = plt.subplots(
        ncols=ncols,
        nrows=nrows,
        figsize=(ncols * 8, nrows * 3),
        gridspec_kw={"wspace": 0.2, "hspace": 0.3},
    )

    if nrows == 1:
        ax = [ax]

    # Make a plot for each of the observer locations
    for iobs, obs_name in enumerate(obs_names):

        # Add blank entry to show colors in legend
        ax[iobs].fill_between([], [], [], color="k", label="Ex")
        ax[iobs].fill_between([], [], [], color="r", label="Ey")
        ax[iobs].fill_between([], [], [], color="b", label="Ez")

        if not iobs:
            ax[iobs].set_title(title)

        # Make the plots for the ZHS waveforms
        if found_zhs and obs_name in observers_zhs.keys():
            obs_position = np.array(
                [
                    observers_zhs[obs_name]["x"],
                    observers_zhs[obs_name]["y"],
                    observers_zhs[obs_name]["z"],
                ]
            ).flatten()
            times = np.array(zhs_dict[obs_name]["time"])

            ax[iobs].plot(
                times, zhs_dict[obs_name]["Ex"], color="k", linestyle="-", label="ZHS"
            )
            ax[iobs].plot(times, zhs_dict[obs_name]["Ey"], color="r", linestyle="-")
            ax[iobs].plot(times, zhs_dict[obs_name]["Ez"], color="b", linestyle="-")
            ax[iobs].set_xlabel(f"Time [{waveforms_config_zhs['units']['time']}]")
            ax[iobs].set_ylabel(
                f"Amp ({waveforms_config_zhs['units']['electric field']})"
            )

        # Make the plots for the CoREAS waveforms
        if found_coreas and obs_name in observers_coreas.keys():
            obs_position = np.array(
                [
                    observers_coreas[obs_name]["x"],
                    observers_coreas[obs_name]["y"],
                    observers_coreas[obs_name]["z"],
                ]
            ).flatten()
            times = np.array(coreas_dict[obs_name]["time"])

            ax[iobs].plot(
                times,
                coreas_dict[obs_name]["Ex"],
                color="k",
                linestyle="--",
                label="CoREAS",
            )
            ax[iobs].plot(times, coreas_dict[obs_name]["Ey"], color="r", linestyle="--")
            ax[iobs].plot(times, coreas_dict[obs_name]["Ez"], color="b", linestyle="--")
            ax[iobs].set_xlabel(f"Time ({waveforms_config_coreas['units']['time']})")
            ax[iobs].set_ylabel(
                f"Amp ({waveforms_config_coreas['units']['electric field']})"
            )

        ymin, ymax = ax[iobs].get_ylim()
        max_amp = max(abs(ymin), abs(ymax))
        ax[iobs].set_ylim(-max_amp, max_amp)
        text = f"Observer @ ({obs_position[0]:0.0f}, {obs_position[1]:0.0f},"
        text += f" {obs_position[2]:0.0f})"
        ax[iobs].text(times[-1], 2 * max_amp * 0.05 - max_amp, text, ha="right")
        ax[iobs].legend(loc="upper right")

    plot_path = os.path.join(args.output_dir, f"ObserverWaveforms_Sh{ishower}.png")
    print("Saving", plot_path)
    fig.savefig(plot_path, bbox_inches="tight")
