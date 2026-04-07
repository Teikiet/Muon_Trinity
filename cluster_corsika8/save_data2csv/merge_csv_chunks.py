#!/usr/bin/env python3
# merge_csv_chunks.py — combine all chunk CSVs into one final CSV

import glob
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    chunk_files = sorted(glob.glob(f"{args.chunk_dir}/result_*.csv"))
    print(f"Found {len(chunk_files)} chunk files")

    total_rows = 0
    with open(args.output, "w") as out:
        for i, cf in enumerate(chunk_files):
            with open(cf) as inp:
                for line_num, line in enumerate(inp):
                    # Write header only from first file
                    if line_num == 0:
                        if i == 0:
                            out.write(line)
                        continue
                    out.write(line)
                    total_rows += 1

    print(f"Merged {total_rows} rows into {args.output}")

if __name__ == "__main__":
    main()
