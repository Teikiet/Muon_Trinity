#!/usr/bin/env python3
"""
particle_location_correction.py

Post-processing script for CORSIKA 8 Cherenkov output in eventio format.
Uses single-pass 1D histogram binning in X and Y separately to find the
Cherenkov photon hotspot, then filters by radius around that peak and
optionally recenters.

Algorithm:
    1) Extract all photon bunch positions (cm) and bunch weights.
    2) Choose one bin count for both X and Y histograms:
             - if --nbins is provided, use it directly;
             - otherwise compute N from area/500^2, where
                 area = pi * R_max^2 and R_max = max(sqrt(x^2 + y^2)).
    3) Build weighted 1D histograms for X and Y separately and take the
         highest-count bin center from each axis as (center_x, center_y).
    4) Filter with user telescope radius and optionally recenter.
"""

import argparse
import sys
import struct
import numpy as np

# ============================================================
# Eventio parsing
# ============================================================

def detect_sync_marker(raw):
    SYNC_BE = b"\xd4\x1f\x8a\x37"
    SYNC_LE = b"\x37\x8a\x1f\xd4"
    pos_be = raw.find(SYNC_BE)
    pos_le = raw.find(SYNC_LE)
    if pos_be >= 0 and (pos_le < 0 or pos_be <= pos_le):
        return SYNC_BE, "big-endian"
    elif pos_le >= 0:
        return SYNC_LE, "little-endian"
    return None, None

def parse_toplevel_header(raw, pos, sync_bytes):
    sync_len = len(sync_bytes)
    if pos + sync_len + 12 > len(raw):
        return None
    offset = pos + sync_len
    type_ver = struct.unpack_from("<I", raw, offset)[0]
    ident = struct.unpack_from("<i", raw, offset + 4)[0]
    length_field = struct.unpack_from("<I", raw, offset + 8)[0]

    block_type = type_ver & 0xFFFF
    version = (type_ver >> 20) & 0xFFF
    extended = bool(type_ver & (1 << 17))
    only_sub = bool(length_field & (1 << 30))

    header_size = sync_len + 12
    if extended:
        if pos + header_size + 4 > len(raw):
            return None
        ext = struct.unpack_from("<I", raw, offset + 12)[0]
        content_length = (length_field & 0x3FFFFFFF) | (ext << 30)
        header_size += 4
    else:
        content_length = length_field & 0x3FFFFFFF

    return {
        'type': block_type, 'version': version, 'id': ident,
        'only_sub': only_sub, 'header_size': header_size,
        'content_length': content_length,
        'content_start': pos + header_size,
        'block_end': pos + header_size + content_length,
        'pos': pos,
    }

def parse_subobject_header(raw, pos):
    if pos + 12 > len(raw):
        return None
    type_ver = struct.unpack_from("<I", raw, pos)[0]
    ident = struct.unpack_from("<i", raw, pos + 4)[0]
    length_field = struct.unpack_from("<I", raw, pos + 8)[0]

    block_type = type_ver & 0xFFFF
    version = (type_ver >> 20) & 0xFFF
    extended = bool(type_ver & (1 << 17))
    only_sub = bool(length_field & (1 << 30))

    header_size = 12
    if extended:
        if pos + 16 > len(raw):
            return None
        ext = struct.unpack_from("<I", raw, pos + 12)[0]
        content_length = (length_field & 0x3FFFFFFF) | (ext << 30)
        header_size = 16
    else:
        content_length = length_field & 0x3FFFFFFF

    return {
        'type': block_type, 'version': version, 'id': ident,
        'only_sub': only_sub, 'header_size': header_size,
        'content_length': content_length,
        'content_start': pos + header_size,
        'block_end': pos + header_size + content_length,
        'pos': pos,
    }

def extract_photon_bunches(raw, sub_info):
    content_start = sub_info['content_start']
    content = raw[content_start:sub_info['block_end']]
    if len(content) < 12:
        return None, False
    n_bunches_raw = struct.unpack_from("<i", content, 8)[0]
    is_compact = (n_bunches_raw < 0)

    if n_bunches_raw > 0:
        n = n_bunches_raw
        data_offset = 12
        if len(content) < data_offset + n * 32:
            return None, False
        arr = np.frombuffer(content[data_offset:data_offset + n * 32],
                            dtype=np.float32).reshape(-1, 8).copy()
        return arr, False
    elif n_bunches_raw < 0:
        n = abs(n_bunches_raw)
        data_offset = 12
        if len(content) < data_offset + n * 16:
            return None, True
        arr = np.frombuffer(content[data_offset:data_offset + n * 16],
                            dtype=np.int16).reshape(-1, 8).astype(np.float32).copy()
        arr[:, 0] *= 0.1
        arr[:, 1] *= 0.1
        return arr, True
    else:
        return np.empty((0, 8), dtype=np.float32), False

def rebuild_subobject_header(sub_info, new_content_length):
    need_ext = (new_content_length > 0x3FFFFFFF)
    type_ver = (sub_info['type'] & 0xFFFF) | (sub_info['version'] << 20)
    if need_ext:
        type_ver |= (1 << 17)
    header = struct.pack("<I", type_ver)
    header += struct.pack("<i", sub_info['id'])
    length_word = new_content_length & 0x3FFFFFFF
    if sub_info['only_sub']:
        length_word |= (1 << 30)
    header += struct.pack("<I", length_word)
    if need_ext:
        header += struct.pack("<I", (new_content_length >> 30) & 0xFFFFFFFF)
    return header

def rebuild_toplevel_header(raw, pos, sync_bytes, new_content_length):
    sync_len = len(sync_bytes)
    offset = pos + sync_len
    type_ver = struct.unpack_from("<I", raw, offset)[0]
    ident_bytes = raw[offset + 4:offset + 8]
    orig_length_field = struct.unpack_from("<I", raw, offset + 8)[0]
    orig_only_sub = bool(orig_length_field & (1 << 30))
    orig_extended = bool(type_ver & (1 << 17))
    need_ext = (new_content_length > 0x3FFFFFFF) or orig_extended

    header = bytearray(sync_bytes)
    if need_ext:
        type_ver |= (1 << 17)
    else:
        type_ver &= ~(1 << 17)
    header += struct.pack("<I", type_ver)
    header += ident_bytes
    length_word = new_content_length & 0x3FFFFFFF
    if orig_only_sub:
        length_word |= (1 << 30)
    header += struct.pack("<I", length_word)
    if need_ext:
        header += struct.pack("<I", (new_content_length >> 30) & 0xFFFFFFFF)
    return bytes(header)

def build_photons_content(bunches, is_compact=False):
    n_out = bunches.shape[0]
    total_photons = float(np.sum(np.abs(bunches[:, 6]))) if n_out > 0 else 0.0
    content = bytearray()
    content += struct.pack("<h", 0)
    content += struct.pack("<h", 0)
    content += struct.pack("<f", total_photons)
    if is_compact:
        content += struct.pack("<i", -n_out)
        if n_out > 0:
            c = bunches.copy()
            c[:, 0] /= 0.1
            c[:, 1] /= 0.1
            content += c.astype(np.int16).tobytes()
    else:
        content += struct.pack("<i", n_out)
        if n_out > 0:
            content += bunches.astype(np.float32).tobytes()
    return bytes(content)

# ============================================================
# Peak finding
# ============================================================

def compute_default_bins_from_area(px, py, bin_area_cm2=500.0 ** 2):
    """Compute default histogram bin count from max-radius area in cm^2."""
    r = np.sqrt(px * px + py * py)
    r_max = float(np.max(r)) if r.size > 0 else 0.0
    area_max = np.pi * r_max * r_max
    n_bins = max(1, int(area_max / bin_area_cm2))
    return n_bins, r_max, area_max


def find_peak_1d_axis(values, weights, n_bins, v_min, v_max):
    """Return weighted 1D histogram peak bin center and diagnostics."""
    if np.isclose(v_max, v_min):
        # Degenerate axis range: all values are equal.
        peak = float(v_min)
        peak_count = float(np.sum(weights))
        hist = np.array([peak_count], dtype=np.float64)
        edges = np.array([v_min, v_max if v_max > v_min else v_min + 1.0], dtype=np.float64)
        return peak, hist, edges, 0, peak_count

    hist, edges = np.histogram(
        values,
        bins=n_bins,
        range=(v_min, v_max),
        weights=weights,
    )
    idx = int(np.argmax(hist))
    peak = 0.5 * (edges[idx] + edges[idx + 1])
    peak_count = float(hist[idx])
    return peak, hist, edges, idx, peak_count

# ============================================================
# Main processing
# ============================================================

def find_toplevel_blocks(raw, sync_bytes):
    positions = []
    start = 0
    while True:
        p = raw.find(sync_bytes, start)
        if p == -1:
            break
        positions.append(p)
        start = p + len(sync_bytes)

    blocks = []
    for sp in positions:
        info = parse_toplevel_header(raw, sp, sync_bytes)
        if info:
            blocks.append(info)
    return blocks

def get_all_photon_bunches(raw, sync_bytes, blocks):
    """Extract all photon bunches from all TelescopeData blocks."""
    all_bunches = []
    for info in blocks:
        if info['type'] != 1204:
            continue
        sub_pos = info['content_start']
        content_end = info['block_end']
        while sub_pos < content_end:
            sub = parse_subobject_header(raw, sub_pos)
            if sub is None:
                break
            if sub['type'] == 1205:
                bunches, _ = extract_photon_bunches(raw, sub)
                if bunches is not None and bunches.shape[0] > 0:
                    all_bunches.append(bunches)
            sub_pos = sub['block_end']
    if all_bunches:
        return np.vstack(all_bunches)
    return np.empty((0, 8), dtype=np.float32)

def filter_and_write(raw, sync_bytes, blocks, output_path,
                     center_x_cm, center_y_cm, radius_cm, recenter,
                     tel_x_cm=0.0, tel_y_cm=0.0):

    """Filter photons by radius and write output."""
    total_in = 0
    total_out = 0
    output_parts = []
    last_end = 0

    for info in blocks:
        if info['pos'] > last_end:
            output_parts.append(raw[last_end:info['pos']])

        if info['type'] == 1204:
            new_block, n_in, n_out = _process_telescope_data(
            raw, info, sync_bytes,
            center_x_cm, center_y_cm, radius_cm, recenter,
            tel_x_cm, tel_y_cm)
            output_parts.append(new_block)
            total_in += n_in
            total_out += n_out
        else:
            output_parts.append(raw[info['pos']:info['block_end']])

        last_end = info['block_end']

    if last_end < len(raw):
        output_parts.append(raw[last_end:])

    output_data = b"".join(output_parts)
    print(f"\nWriting: {output_path}")
    with open(output_path, "wb") as f:
        f.write(output_data)
    print(f"  Output size: {len(output_data)} bytes")

    return total_in, total_out

def _process_telescope_data(raw, block_info, sync_bytes,
                            center_x_cm, center_y_cm, radius_cm, recenter,
                            tel_x_cm=0.0, tel_y_cm=0.0):

    content_start = block_info['content_start']
    content_end = block_info['block_end']

    total_in = 0
    total_out = 0
    new_subobjects = bytearray()
    sub_pos = content_start

    while sub_pos < content_end:
        sub = parse_subobject_header(raw, sub_pos)
        if sub is None:
            new_subobjects += raw[sub_pos:content_end]
            break

        if sub['type'] == 1205:
            bunches, is_compact = extract_photon_bunches(raw, sub)
            if bunches is not None and bunches.shape[0] > 0:
                n_in = bunches.shape[0]
                dx = bunches[:, 0] - center_x_cm
                dy = bunches[:, 1] - center_y_cm
                r_sq = dx * dx + dy * dy
                mask = r_sq <= radius_cm * radius_cm
                filtered = bunches[mask]

                if recenter and filtered.shape[0] > 0:
                    filtered[:, 0] -= (center_x_cm - tel_x_cm)
                    filtered[:, 1] -= (center_y_cm - tel_y_cm)


                n_out = filtered.shape[0]
                total_in += n_in
                total_out += n_out

                new_content = build_photons_content(filtered, is_compact=is_compact)
                new_sub_header = rebuild_subobject_header(sub, len(new_content))
                new_subobjects += new_sub_header + new_content
            else:
                new_subobjects += raw[sub_pos:sub['block_end']]
        else:
            new_subobjects += raw[sub_pos:sub['block_end']]

        sub_pos = sub['block_end']

    new_header = rebuild_toplevel_header(
        raw, block_info['pos'], sync_bytes, len(new_subobjects))

    return bytes(new_header) + bytes(new_subobjects), total_in, total_out

def main():
    parser = argparse.ArgumentParser(
        description="Filter CORSIKA 8 Cherenkov eventio using single-pass 1D X/Y "
                    "histogram peak-finding to locate the Cherenkov photon hotspot, then "
                    "filter by radius and optionally recenter.")

    parser.add_argument("--cherenkov-eventio", required=True,
                        help="Input Cherenkov eventio file")
    parser.add_argument("--telescope-radius", type=float, required=True,
                        help="Final filter radius in meters")
    parser.add_argument("--output", required=True,
                        help="Output filtered eventio file")
    parser.add_argument("--recenter", action="store_true", default=False,
                        help="Recenter photon x,y so peak moves to telescope position")
    parser.add_argument("--telescope-x", type=float, default=0.0,
                        help="Telescope X position in meters (recenter target)")
    parser.add_argument("--telescope-y", type=float, default=0.0,
                        help="Telescope Y position in meters (recenter target)")
    parser.add_argument("--nbins", type=int, default=None,
                        help="Bins per axis for 1D X/Y histograms (override default area/500^2 rule)")

    args = parser.parse_args()


    print("=" * 60)
    print("Cherenkov Photon Peak Finder & Filter (eventio)")
    print("=" * 60)

    # --- Step 1: Read eventio and extract all photon bunches ---
    print(f"\nReading: {args.cherenkov_eventio}")
    with open(args.cherenkov_eventio, "rb") as f:
        raw = f.read()
    print(f"  File size: {len(raw)} bytes")

    sync_bytes, endian_desc = detect_sync_marker(raw)
    if sync_bytes is None:
        print("ERROR: No sync markers found!", file=sys.stderr)
        sys.exit(1)
    print(f"  Sync: {' '.join(f'{b:02x}' for b in sync_bytes)} ({endian_desc})")

    blocks = find_toplevel_blocks(raw, sync_bytes)
    print(f"  Blocks: {' '.join(str(b['type']) for b in blocks)}")

    all_photons = get_all_photon_bunches(raw, sync_bytes, blocks)
    n_total = all_photons.shape[0]
    print(f"  Total photon bunches: {n_total}")

    if n_total == 0:
        print("WARNING: No photon bunches found! Copying file as-is.")
        with open(args.output, "wb") as f:
            f.write(raw)
        return

    # Photon positions and weights (eventio coordinates, cm)
    px = all_photons[:, 0]
    py = all_photons[:, 1]
    pw = np.abs(all_photons[:, 6])

    x_min_full = np.min(px)
    x_max_full = np.max(px)
    y_min_full = np.min(py)
    y_max_full = np.max(py)

    # Compute max radius from origin for default bin-count estimate.
    n_bins_default, r_max_cm, area_max_cm2 = compute_default_bins_from_area(px, py)
    if args.nbins is not None:
        n_bins = max(1, args.nbins)
        binning_mode = "user-defined"
    else:
        n_bins = n_bins_default
        binning_mode = "auto (area/500^2)"

    print(f"\n  Photon x range: [{x_min_full:.1f}, {x_max_full:.1f}] cm")
    print(f"  Photon y range: [{y_min_full:.1f}, {y_max_full:.1f}] cm")
    print(f"  Max hit radius from origin: {r_max_cm:.1f} cm = {r_max_cm/100.0:.2f} m")
    print(f"  Max hit area from origin:   {area_max_cm2:.1f} cm^2")
    print(f"  Histogram bins/axis mode:   {binning_mode}")
    print(f"  Histogram bins/axis used:   {n_bins}")

    # ================================================================
    # Single-pass 1D X/Y peak finding
    # ================================================================
    print("\n--- Peak finding: weighted 1D histograms in X and Y ---")

    peak_x_cm, x_hist, x_edges, x_peak_idx, x_peak_count = find_peak_1d_axis(
        px, pw, n_bins, x_min_full, x_max_full)
    peak_y_cm, y_hist, y_edges, y_peak_idx, y_peak_count = find_peak_1d_axis(
        py, pw, n_bins, y_min_full, y_max_full)

    bw_x = x_edges[1] - x_edges[0] if len(x_edges) > 1 else 0.0
    bw_y = y_edges[1] - y_edges[0] if len(y_edges) > 1 else 0.0

    print(f"  X bin size: {bw_x:.1f} cm = {bw_x/100.0:.2f} m")
    print(f"  Y bin size: {bw_y:.1f} cm = {bw_y/100.0:.2f} m")
    print(f"  X peak bin index: {x_peak_idx}")
    print(f"  Y peak bin index: {y_peak_idx}")
    print(f"  X peak weighted count: {x_peak_count:.1f}")
    print(f"  Y peak weighted count: {y_peak_count:.1f}")
    print(f"  Peak center: ({peak_x_cm:.1f}, {peak_y_cm:.1f}) cm "
          f"= ({peak_x_cm/100.0:.2f}, {peak_y_cm/100.0:.2f}) m")

    center_x_cm = peak_x_cm
    center_y_cm = peak_y_cm

    # ================================================================
    # Final filter with user-specified telescope radius
    # ================================================================
    radius_cm = args.telescope_radius * 100.0

    print(f"\n--- Final filter ---")
    print(f"  Center: ({center_x_cm:.1f}, {center_y_cm:.1f}) cm "
          f"= ({center_x_cm/100:.2f}, {center_y_cm/100:.2f}) m")
    print(f"  Radius: {radius_cm:.1f} cm = {args.telescope_radius:.1f} m")
    if args.recenter:
        print(f"  Recenter target (tel):   ({args.telescope_x:.2f}, {args.telescope_y:.2f}) m")


    # Preview
    dx_final = px - center_x_cm
    dy_final = py - center_y_cm
    r_sq_final = dx_final * dx_final + dy_final * dy_final
    n_keep = np.sum(r_sq_final <= radius_cm ** 2)
    print(f"  Photons within final radius: {n_keep}/{n_total} "
          f"({n_keep/n_total*100:.1f}%)")

    # --- Write filtered output ---
    total_in, total_out = filter_and_write(
        raw, sync_bytes, blocks, args.output,
        center_x_cm, center_y_cm, radius_cm, args.recenter,
        args.telescope_x * 100.0, args.telescope_y * 100.0)


    print()
    print("=" * 60)
    print("Summary:")
    print(f"  Input photon bunches:    {total_in}")
    print(f"  Filtered photon bunches: {total_out}")
    if total_in > 0:
        print(f"  Retention rate:          {total_out / total_in * 100:.1f}%")
    print(f"  Peak method:             weighted 1D X/Y histogram")
    print(f"  Bin count / axis:        {n_bins}")
    print(f"  X/Y bin size:            {bw_x/100.0:.2f} x {bw_y/100.0:.2f} m")
    print(f"  X peak index/count:      {x_peak_idx} / {x_peak_count:.1f}")
    print(f"  Y peak index/count:      {y_peak_idx} / {y_peak_count:.1f}")
    print(f"  Final center:            ({center_x_cm:.1f}, {center_y_cm:.1f}) cm "
          f"= ({center_x_cm/100:.2f}, {center_y_cm/100:.2f}) m")
    print(f"  Final filter radius:     {args.telescope_radius:.1f} m")
    print(f"  Recentered:              {args.recenter}")
    print(f"  Output:                  {args.output}")
    print("=" * 60)

if __name__ == "__main__":
    main()
