#!/usr/bin/env python3
"""
particle_location_correction.py

Post-processing script for CORSIKA 8 Cherenkov output in eventio format.
Uses two-pass 100x100 histogram binning to find the Cherenkov photon peak,
then filters by radius around that peak and optionally recenters.

Algorithm:
  Pass 1: Bin all photons into 100x100 grid over full range.
           Find peak bin. Use peak center + half the original photon radius
           as the coarse region.
  Pass 2: Bin photons within coarse region into 100x100 grid.
           Find refined peak bin. Use that as the final center with
           the user-specified telescope radius.
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
# Two-pass histogram peak finding
# ============================================================

def find_peak_histogram(x, y, weights, x_min, x_max, y_min, y_max, nbins=100):
    """
    Bin photons into nbins x nbins histogram and return the center
    of the bin with the highest weighted count.

    Returns:
        peak_x, peak_y: center of the peak bin
        bin_width_x, bin_width_y: size of each bin
        hist: the 2D histogram array
        xedges, yedges: bin edges
    """
    # Create histogram
    hist, xedges, yedges = np.histogram2d(
        x, y, bins=nbins,
        range=[[x_min, x_max], [y_min, y_max]],
        weights=weights)

    # Find peak bin
    peak_idx = np.unravel_index(np.argmax(hist), hist.shape)
    ix, iy = peak_idx

    # Center of peak bin
    peak_x = 0.5 * (xedges[ix] + xedges[ix + 1])
    peak_y = 0.5 * (yedges[iy] + yedges[iy + 1])

    bin_width_x = xedges[1] - xedges[0]
    bin_width_y = yedges[1] - yedges[0]

    return peak_x, peak_y, bin_width_x, bin_width_y, hist, xedges, yedges

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
                     center_x_cm, center_y_cm, radius_cm, recenter):
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
                center_x_cm, center_y_cm, radius_cm, recenter)
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
                            center_x_cm, center_y_cm, radius_cm, recenter):
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
                    filtered[:, 0] -= center_x_cm
                    filtered[:, 1] -= center_y_cm

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
        description="Filter CORSIKA 8 Cherenkov eventio using two-pass histogram "
                    "peak-finding to locate the Cherenkov photon hotspot, then "
                    "filter by radius and optionally recenter.")

    parser.add_argument("--cherenkov-eventio", required=True,
                        help="Input Cherenkov eventio file")
    parser.add_argument("--telescope-radius", type=float, required=True,
                        help="Final filter radius in meters (used in pass 2)")
    parser.add_argument("--output", required=True,
                        help="Output filtered eventio file")
    parser.add_argument("--recenter", action="store_true", default=False,
                        help="Recenter photon x,y so peak is at origin")
    parser.add_argument("--nbins", type=int, default=100,
                        help="Number of bins per axis for histogram (default: 100)")

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

    # Compute the max photon radius (half-span of the photon distribution)
    x_span_cm = x_max_full - x_min_full
    y_span_cm = y_max_full - y_min_full
    max_photon_radius_cm = 0.5 * max(x_span_cm, y_span_cm)
    max_photon_radius_m = max_photon_radius_cm / 100.0

    print(f"\n  Photon x range: [{x_min_full:.1f}, {x_max_full:.1f}] cm")
    print(f"  Photon y range: [{y_min_full:.1f}, {y_max_full:.1f}] cm")
    print(f"  Max photon radius (half-span): {max_photon_radius_cm:.1f} cm "
          f"= {max_photon_radius_m:.1f} m")

    # ================================================================
    # PASS 1: Coarse histogram over full photon range
    # ================================================================
    print(f"\n--- Pass 1: Coarse binning ({args.nbins}x{args.nbins}) over full range ---")

    peak1_x, peak1_y, bw1_x, bw1_y, hist1, xe1, ye1 = find_peak_histogram(
        px, py, pw,
        x_min_full, x_max_full, y_min_full, y_max_full,
        nbins=args.nbins)

    # Find peak bin indices for reporting
    peak1_idx = np.unravel_index(np.argmax(hist1), hist1.shape)
    peak1_count = hist1[peak1_idx]

    print(f"  Bin size: {bw1_x:.1f} x {bw1_y:.1f} cm "
          f"= {bw1_x/100:.2f} x {bw1_y/100:.2f} m")
    print(f"  Peak bin index: ({peak1_idx[0]}, {peak1_idx[1]})")
    print(f"  Peak bin weighted count: {peak1_count:.1f}")
    print(f"  Peak center: ({peak1_x:.1f}, {peak1_y:.1f}) cm "
          f"= ({peak1_x/100:.2f}, {peak1_y/100:.2f}) m")

    # Coarse cut radius = half the max photon radius
    coarse_radius_cm = max_photon_radius_cm / 2.0
    coarse_radius_m = coarse_radius_cm / 100.0
    print(f"  Coarse cut radius: {coarse_radius_cm:.1f} cm = {coarse_radius_m:.1f} m")

    # Select photons within coarse radius of pass-1 peak
    dx1 = px - peak1_x
    dy1 = py - peak1_y
    r_sq1 = dx1 * dx1 + dy1 * dy1
    coarse_mask = r_sq1 <= coarse_radius_cm ** 2
    n_coarse = np.sum(coarse_mask)
    print(f"  Photons within coarse radius: {n_coarse}/{n_total} "
          f"({n_coarse/n_total*100:.1f}%)")

    if n_coarse == 0:
        print("WARNING: No photons in coarse region! Using pass-1 peak as final center.")
        center_x_cm = peak1_x
        center_y_cm = peak1_y
    else:
        # ================================================================
        # PASS 2: Fine histogram over coarse region
        # ================================================================
        print(f"\n--- Pass 2: Fine binning ({args.nbins}x{args.nbins}) "
              f"within coarse region ---")

        px_coarse = px[coarse_mask]
        py_coarse = py[coarse_mask]
        pw_coarse = pw[coarse_mask]

        # Define the fine binning range: square centered on pass-1 peak
        fine_x_min = peak1_x - coarse_radius_cm
        fine_x_max = peak1_x + coarse_radius_cm
        fine_y_min = peak1_y - coarse_radius_cm
        fine_y_max = peak1_y + coarse_radius_cm

        peak2_x, peak2_y, bw2_x, bw2_y, hist2, xe2, ye2 = find_peak_histogram(
            px_coarse, py_coarse, pw_coarse,
            fine_x_min, fine_x_max, fine_y_min, fine_y_max,
            nbins=args.nbins)

        peak2_idx = np.unravel_index(np.argmax(hist2), hist2.shape)
        peak2_count = hist2[peak2_idx]

        print(f"  Bin size: {bw2_x:.1f} x {bw2_y:.1f} cm "
              f"= {bw2_x/100:.2f} x {bw2_y/100:.2f} m")
        print(f"  Peak bin index: ({peak2_idx[0]}, {peak2_idx[1]})")
        print(f"  Peak bin weighted count: {peak2_count:.1f}")
        print(f"  Refined peak center: ({peak2_x:.1f}, {peak2_y:.1f}) cm "
              f"= ({peak2_x/100:.2f}, {peak2_y/100:.2f}) m")

        # Shift from pass-1 to pass-2
        shift = np.sqrt((peak2_x - peak1_x)**2 + (peak2_y - peak1_y)**2)
        print(f"  Shift from pass-1 to pass-2: {shift:.1f} cm = {shift/100:.2f} m")

        center_x_cm = peak2_x
        center_y_cm = peak2_y

    # ================================================================
    # Final filter with user-specified telescope radius
    # ================================================================
    radius_cm = args.telescope_radius * 100.0

    print(f"\n--- Final filter ---")
    print(f"  Center: ({center_x_cm:.1f}, {center_y_cm:.1f}) cm "
          f"= ({center_x_cm/100:.2f}, {center_y_cm/100:.2f}) m")
    print(f"  Radius: {radius_cm:.1f} cm = {args.telescope_radius:.1f} m")
    if args.recenter:
        print(f"  Will recenter photons to origin")

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
        center_x_cm, center_y_cm, radius_cm, args.recenter)

    print()
    print("=" * 60)
    print("Summary:")
    print(f"  Input photon bunches:    {total_in}")
    print(f"  Filtered photon bunches: {total_out}")
    if total_in > 0:
        print(f"  Retention rate:          {total_out / total_in * 100:.1f}%")
    print(f"  Pass-1 peak:             ({peak1_x:.1f}, {peak1_y:.1f}) cm")
    print(f"  Pass-1 coarse radius:    {coarse_radius_m:.1f} m "
          f"(half of {max_photon_radius_m:.1f} m)")
    if n_coarse > 0:
        print(f"  Pass-2 refined peak:     ({center_x_cm:.1f}, {center_y_cm:.1f}) cm")
    print(f"  Final center:            ({center_x_cm:.1f}, {center_y_cm:.1f}) cm "
          f"= ({center_x_cm/100:.2f}, {center_y_cm/100:.2f}) m")
    print(f"  Final filter radius:     {args.telescope_radius:.1f} m")
    print(f"  Recentered:              {args.recenter}")
    print(f"  Output:                  {args.output}")
    print("=" * 60)

if __name__ == "__main__":
    main()
