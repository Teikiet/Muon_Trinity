#!/usr/bin/env python3

import argparse
import json
import logging
import math
import os

import numpy as np


def to_1e(val):
    exp = int(math.log10(val))
    coeff = val / 10**exp
    return f"{coeff:g}e{exp}"


def generate_energy_strings():
    """os.environ["MCEQ_LOG_LEVEL"] = "50"
    logging.disable(logging.CRITICAL)
    import MCEq.config as config
    from MCEq.core import MCEqRun
    import crflux.models as pm

    config.kernel_config = "MKL"
    config.e_min = 100
    config.integrator = "euler"

    mceq = MCEqRun(
        interaction_model="SIBYLL23C",
        primary_model=(pm.GlobalSplineFitBeta, None),
        theta_deg=0.0,
        density_model=("CORSIKA", ("USStd", None)),
    )

    energies = mceq.e_grid"""
    #All energy range from MCEq
    energies = np.array([np.float64(0.8912509381337457), np.float64(1.1220184543019636), np.float64(1.4125375446227548), np.float64(1.7782794100389236), np.float64(2.238721138568341), np.float64(2.818382931264455), np.float64(3.548133892335755), np.float64(4.466835921509633), np.float64(5.623413251903493), np.float64(7.079457843841384), np.float64(8.91250938133746), np.float64(11.220184543019636), np.float64(14.125375446227547), np.float64(17.782794100389236), np.float64(22.38721138568341), np.float64(28.183829312644548), np.float64(35.48133892335755), np.float64(44.668359215096324), np.float64(56.23413251903494), np.float64(70.79457843841384), np.float64(89.12509381337459), np.float64(112.20184543019641), np.float64(141.25375446227557), np.float64(177.82794100389228), np.float64(223.872113856834), np.float64(281.83829312644554), np.float64(354.81338923357566), np.float64(446.6835921509635), np.float64(562.3413251903496), np.float64(707.9457843841387), np.float64(891.2509381337459), np.float64(1122.018454301964), np.float64(1412.5375446227554), np.float64(1778.2794100389247), np.float64(2238.721138568342), np.float64(2818.3829312644552), np.float64(3548.133892335757), np.float64(4466.835921509635), np.float64(5623.4132519034965), np.float64(7079.457843841388), np.float64(8912.509381337459), np.float64(11220.18454301964), np.float64(14125.375446227554), np.float64(17782.794100389245), np.float64(22387.211385683426), np.float64(28183.82931264455), np.float64(35481.338923357565), np.float64(44668.359215096345), np.float64(56234.13251903497), np.float64(70794.57843841387), np.float64(89125.09381337458), np.float64(112201.84543019641), np.float64(141253.75446227554), np.float64(177827.94100389248), np.float64(223872.11385683424), np.float64(281838.2931264455), np.float64(354813.38923357567), np.float64(446683.5921509635), np.float64(562341.3251903496), np.float64(707945.7843841388), np.float64(891250.9381337459), np.float64(1122018.4543019629), np.float64(1412537.5446227554), np.float64(1778279.4100389264), np.float64(2238721.138568342), np.float64(2818382.931264455), np.float64(3548133.8923357534), np.float64(4466835.921509635), np.float64(5623413.251903502), np.float64(7079457.843841388), np.float64(8912509.381337458), np.float64(11220184.54301963), np.float64(14125375.446227556), np.float64(17782794.100389265), np.float64(22387211.385683425), np.float64(28183829.31264455), np.float64(35481338.923357606), np.float64(44668359.21509644), np.float64(56234132.51903503), np.float64(70794578.43841387), np.float64(89125093.8133746), np.float64(112201845.43019652), np.float64(141253754.46227583), np.float64(177827941.00389266), np.float64(223872113.85683423), np.float64(281838293.12644553), np.float64(354813389.23357606), np.float64(446683592.1509644), np.float64(562341325.1903503), np.float64(707945784.3841388), np.float64(891250938.1337459), np.float64(1122018454.3019652), np.float64(1412537544.6227584), np.float64(1778279410.0389264), np.float64(2238721138.568342), np.float64(2818382931.264455), np.float64(3548133892.33576), np.float64(4466835921.5096445), np.float64(5623413251.903502), np.float64(7079457843.841388), np.float64(8912509381.33746), np.float64(11220184543.019653), np.float64(14125375446.227583), np.float64(17782794100.389263), np.float64(22387211385.683422), np.float64(28183829312.644547), np.float64(35481338923.357605), np.float64(44668359215.096436), np.float64(56234132519.03502), np.float64(70794578438.41388), np.float64(89125093813.37459)])
    #Reduce energy range:
    #energies_1 = energies[(energies >= 1e1) & (energies <= 1e2)]
    #energies_2 = energies[(energies >= 1e3) & (energies <= 1e4)]
    #energies = np.concatenate((energies_1, energies_2), axis = 0)
    energies = energies[(energies >= 1e4) & (energies <= 1e5)]
    return [to_1e(e) for e in energies]


def build_geometry_defaults():
    zeniths = [f"{z:.1f}" for z in np.arange(87.0, 90.0, 0.3)] #
    zeniths_full = [f"{z:.1f}" for z in np.arange(86.1, 92.0, 0.3)]
    if "89.9" not in zeniths:
        zeniths.append("89.9")
        zeniths_full.append("89.9")
    zeniths_extended = list(set(zeniths_full)-set(zeniths))

    azimuths = [f"{a:.1f}" for a in np.arange(267.0, 273.3, 0.3)] #["270.0"]#
    azimuths_full = [f"{a:.1f}" for a in np.arange(266.1, 274.3, 0.3)]
    azimuths_extended = list(set(azimuths_full) - set(azimuths)) 

    heights  = [str(int(h)) for h in np.linspace(5000, 50000, 100)] #["5000", "30000","100000"]#
    heigths_extended =  [str(int(h)) for h in np.logspace(np.log10(5000), np.log10(100000), 20)]
    height_full = heights + heigths_extended
    return {
        "zeniths_deg": zeniths_full, # zeniths_full zeniths_extended zeniths
        "azimuths_deg": azimuths_full, # azimuths_full azimuths_extended azimuths
        "heights_m": height_full, # height_full heigths_extended heights
        "tel_xs_m": ["0.3", "1.0", "2.0", "3.0", "4.0", "5.0"], #["0","-0.1", "-0.2", "-0.3", "-0.4", "-0.5", "-0.6", "-0.7", "-0.8", "-0.9","-1.0","-2.0","-5.0"],
        "tel_zs_m": ["0.3", "1.0", "2.0", "3.0", "4.0", "5.0"], #["0","-0.1", "-0.2", "-0.3", "-0.4", "-0.5", "-0.6", "-0.7", "-0.8", "-0.9","-1.0","-2.0","-5.0"],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate simulation inputs for Trinity runs."
    )
    parser.add_argument("--output", required=True, help="Path to output JSON")
    args = parser.parse_args()

    data = {
        "pdg": [13],
        "seeds": [1],
        "tel_radius": 5,
        "tel_y": 0,
        "energy_strings": generate_energy_strings(),
        "geometry": build_geometry_defaults(),
    }

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
