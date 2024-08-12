import argparse

def trf_to_bed(trf_output, bed_output):
    with open(trf_output, 'r') as infile, open(bed_output, 'w') as outfile:
        for line in infile:
            if line.startswith('@'):  # Skip header lines
                chrom = line.strip().split()[0][1:]  # Extract chromosome name
                continue
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            start = int(fields[0]) - 1  # Convert to 0-based
            end = int(fields[1])
            outfile.write(f"{chrom}\t{start}\t{end}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert TRF tandem regions output to BED format.")
    parser.add_argument('trf_output', type=str, help='Input TRF output file')
    parser.add_argument('bed_output', type=str, help='Output BED file')

    args = parser.parse_args()

    trf_to_bed(args.trf_output, args.bed_output)

if __name__ == "__main__":
    main()