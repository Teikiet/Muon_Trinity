import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.offline import plot
from read_telescope import survive_photon, VAtmosAbsorption, get_photon_bunches, simulate_davies_cotton_optics

import numpy as np
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
from collections import defaultdict
from scipy.stats import linregress
def get_lambda_fit():
    # Data from user
    x = np.array([
        0.7644305772230889,
    1.5288611544461779,
    2.293291731669267,
    3.2761310452418098,
    4.149765990639626,
    5.6786271450858035,
    6.1154446177847115,
    7.316692667706708,
    6.9890795631825275,
    8.190327613104525,
    9.391575663026522,
    10.265210608424338,
    11.24804992199688,
    12.230889235569423,
    13.213728549141965,
    14.087363494539781

    ])

    y = np.array([
        6.674311926605505,
        5.8348623853211015,
        5.277522935779817,
        4.727064220183487,
        4.217889908256881,
        3.729357798165138,
        3.295871559633028,
        2.883027523,
        2.793577981651376,
        2.3738532110091746,
        1.926605504587156,
        1.4862385321100917,
        1.1972477064220184,
        0.805045872,
        0,
        -0.275229358

    ])

    # Keep only positive y
    x_fit = x
    y_fit = np.log(10**y)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x_fit, y_fit)

    # Compute lambda
    lambda_fit = -1 / slope
    return lambda_fit

def noise_function(H, lambda_fit = 0.8946330713407988):
    size = H[0].shape[1]
    noise_matrix = []
    for i in range(len(H)):
        noise = np.random.exponential(scale=1/lambda_fit, size=(size,size))
        noise_matrix.append(noise)
    return np.array(noise_matrix)

def camera_response_function(t, t_photon, rise=5, fall=15):
    signal = np.zeros_like(t)
    for i, ti in enumerate(t):
        if ti < t_photon:
            signal[i] = np.exp(-(t_photon - ti)**2 / (2 * rise**2))
        else:
            signal[i] = np.exp(-(ti - t_photon)**2 / (2 * fall**2))
    if np.max(signal) == 0:
        print("Warning: PE signal is zero, check the time range and t0 values.")
        print("t0:", t_photon)
        return signal
    else:
        signal /= np.max(signal)
        #signal *= 24.1 #1PE = 24.1 DC
        #signal *= 0.1 #1DC = 0.1 mV
        return signal #in mV

import numpy as np
import plotly.graph_objects as go


def get_bin_centers(edges):
    return (edges[:-1] + edges[1:]) / 2

def make_hexbin_figure(
    H,              # array of shape (n_times, ny, nx)
    Times,          # length‑n_times array
    x_edges,        # length‑(nx+1)
    y_edges,        # length‑(ny+1)
    bin_PE_threshold=20,
    scale="log",
    hexa_size=15,
    name=None, color_theme="tempo", PE_max=350
):
    h = []
    T = []
    for i in range(len(Times)):
        if (H[i]).max() >= bin_PE_threshold:
            h.append(H[i])
            T.append(Times[i])
    H = np.asarray(h)
    Times = np.asarray(T)
    print(f"Number of frames with PE ≥ {bin_PE_threshold}: {len(H)}")
    n_times, ny, nx = H.shape
    Times = Times - Times[0]
    assert len(Times) == n_times
    assert len(x_edges) == nx + 1 and len(y_edges) == ny + 1

    # bin centers + hive offsets
    x_centers = get_bin_centers(x_edges)
    y_centers = get_bin_centers(y_edges)
    dx = x_centers[1] - x_centers[0]

    hex_x = np.zeros((ny, nx))
    hex_y = np.zeros((ny, nx))
    for j in range(ny):
        offset = (j % 2) * (dx / 2)
        hex_x[j, :] = x_centers + offset
        hex_y[j, :] = y_centers[j]
    flat_x = hex_x.flatten()
    flat_y = hex_y.flatten()

    # global color scale
    if PE_max == None:
        all_counts = H.flatten()
        high = all_counts[all_counts >= bin_PE_threshold]
    else:
        high = np.array([PE_max]) #Saturate response
    if scale == "log" and high.size:
        gmin = np.log10(bin_PE_threshold)
        gmax = np.log10(max(high.max(), bin_PE_threshold))
    else:
        gmin = bin_PE_threshold
        gmax = max(high.max() if high.size else bin_PE_threshold, bin_PE_threshold)

    # build frames
    frames = []
    for i, t in enumerate(Times):
        counts = H[i].flatten()
        hi_mask = counts >= bin_PE_threshold
        if scale == "log":
            col = np.log10(counts[hi_mask])
        else:
            col = counts[hi_mask]

        frames.append(go.Frame(
            data=[
                # above threshold: use shared coloraxis
                go.Scatter(
                    x=flat_x[hi_mask],
                    y=flat_y[hi_mask],
                    mode='markers',
                    marker=dict(
                        symbol='hexagon',
                        size=hexa_size,
                        color=col,
                    ),
                    marker_coloraxis="coloraxis",
                    hovertemplate="X: %{x:.1f}<br>Y: %{y:.1f}<br>PE: %{marker.color}<extra></extra>"
                ),
                # below threshold
                go.Scatter(
                    x=flat_x[~hi_mask],
                    y=flat_y[~hi_mask],
                    mode='markers',
                    marker=dict(
                        symbol='hexagon',
                        size=hexa_size,
                        color='white',
                        line=dict(width=0)
                    ),
                    hovertemplate="X: %{x:.1f}<br>Y: %{y:.1f}<br>PE: %{col} below threshold<extra></extra>"
                )
            ],
            name=f"{t:.2f}"
        ))

    # initial traces
    fig = go.Figure(data=frames[0].data, frames=frames)

    # slider
    steps = [
        dict(
            method="animate",
            args=[[f"{t:.2f}"], {"mode":"immediate","frame":{"duration":0},"transition":{"duration":0}}],
            label=f"{t:.1f}"
        )
        for t in Times
    ]

    fig.update_layout(
        #template="plotly_dark",
        width=700, height=700,
        xaxis=dict(title="X (cm)", scaleanchor="y", scaleratio=1),
        yaxis=dict(title="Y (cm)"),
        sliders=[dict(active=0, pad={"t":50}, steps=steps, currentvalue={"prefix":"Time: "})],
        title=(f"{name}"
            f"Hexagonal‐Bin Photon Map (PE ≥ {bin_PE_threshold})<br>"
            f"Scale: {scale}, Hex size: {hexa_size}<br>"
            f"<br>"
        ),
        # shared coloraxis lives here and never disappears
        coloraxis=dict(
            colorscale=color_theme,
            cmin=gmin,
            cmax=gmax,
            colorbar=dict(
                title=f"{'log₁₀' if scale=='log' else ''}(PE ≥ {bin_PE_threshold})",
                x=1.02,
                len=0.8
            )
        )
    )

    return fig

def make_bin_figure(
    H,              # array of shape (n_times, ny, nx)
    Times,          # length‑n_times array
    x_edges,        # length‑(nx+1)
    y_edges,        # length‑(ny+1)
    bin_PE_threshold=20,
    scale="log",
    hexa_size=15,
    name=None, color_theme="tempo", PE_max=350
):
    h = []
    T = []
    for i in range(len(Times)):
        if (H[i]).max() >= bin_PE_threshold:
            h.append(H[i])
            T.append(Times[i])
    H = np.asarray(h)
    Times = np.asarray(T)
    print(f"Number of frames with PE ≥ {bin_PE_threshold}: {len(H)}")
    n_times, ny, nx = H.shape
    Times = Times - Times[0]
    assert len(Times) == n_times
    assert len(x_edges) == nx + 1 and len(y_edges) == ny + 1

    # bin centers + hive offsets
    # bin centers for regular grid
    x_centers = get_bin_centers(x_edges)
    y_centers = get_bin_centers(y_edges)

    # Create regular grid
    grid_x, grid_y = np.meshgrid(x_centers, y_centers)
    flat_x = grid_x.flatten()
    flat_y = grid_y.flatten()

    # global color scale
    if PE_max == None:
        all_counts = H.flatten()
        high = all_counts[all_counts >= bin_PE_threshold]
    else:
        high = np.array([PE_max]) #Saturate response
    if scale == "log" and high.size:
        gmin = np.log10(bin_PE_threshold)
        gmax = np.log10(max(high.max(), bin_PE_threshold))
    else:
        gmin = bin_PE_threshold
        gmax = max(high.max() if high.size else bin_PE_threshold, bin_PE_threshold)

    # build frames
    frames = []
    for i, t in enumerate(Times):
        counts = H[i].flatten()
        hi_mask = counts >= bin_PE_threshold
        if scale == "log":
            col = np.log10(counts[hi_mask])
        else:
            col = counts[hi_mask]

        frames.append(go.Frame(
            data=[
                # above threshold: use shared coloraxis
                go.Scatter(
                    x=flat_x[hi_mask],
                    y=flat_y[hi_mask],
                    mode='markers',
                    marker=dict(
                        symbol='square',
                        size=hexa_size,
                        color=col,
                    ),
                    marker_coloraxis="coloraxis",
                    hovertemplate="X: %{x:.1f}<br>Y: %{y:.1f}<br>PE: %{marker.color}<extra></extra>"
                ),
                # below threshold
                go.Scatter(
                    x=flat_x[~hi_mask],
                    y=flat_y[~hi_mask],
                    mode='markers',
                    marker=dict(
                        symbol='square',
                        size=hexa_size,
                        color='white',
                        line=dict(width=0)
                    ),
                    hovertemplate="X: %{x:.1f}<br>Y: %{y:.1f}<br>PE: %{col} below threshold<extra></extra>"
                )
            ],
            name=f"{t:.2f}"
        ))

    # initial traces
    fig = go.Figure(data=frames[0].data, frames=frames)

    # slider
    steps = [
        dict(
            method="animate",
            args=[[f"{t:.2f}"], {"mode":"immediate","frame":{"duration":0},"transition":{"duration":0}}],
            label=f"{t:.1f}"
        )
        for t in Times
    ]

    fig.update_layout(
        #template="plotly_dark",
        width=700, height=700,
        xaxis=dict(title="X (cm)", scaleanchor="y", scaleratio=1),
        yaxis=dict(title="Y (cm)"),
        sliders=[dict(active=0, pad={"t":50}, steps=steps, currentvalue={"prefix":"Time: "})],
        title=(f"{name}"
            f"Camera Photon Map (PE ≥ {bin_PE_threshold})<br>"
            f"Scale: {scale}"
            f"<br>"
        ),
        # shared coloraxis lives here and never disappears
        coloraxis=dict(
            colorscale=color_theme,
            cmin=gmin,
            cmax=gmax,
            colorbar=dict(
                title=f"{'log₁₀' if scale=='log' else ''}(PE ≥ {bin_PE_threshold})",
                x=1.02,
                len=0.8
            )
        )
    )

    return fig