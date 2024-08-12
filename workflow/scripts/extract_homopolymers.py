import argparse

def extract_homopolymers(reference_file, min_length, output_file):
    with open(reference_file, 'r') as ref, open(output_file, 'w') as out:
        chrom = ""
        seq = ""
        
        def process_seq(sequence, chromosome):
            start = 0
            count = 1
            last_base = sequence[0]
            length = len(sequence)
            
            for i in range(1, length):
                base = sequence[i]
                if base != last_base:
                    if count >= min_length:
                        out.write(f"{chromosome}\t{start}\t{i}\n")
                    last_base = base
                    start = i
                    count = 1
                else:
                    count += 1
            
            if count >= min_length:
                out.write(f"{chromosome}\t{start}\t{length}\n")
        
        for line in ref:
            if line.startswith('>'):
                if seq and chrom:
                    process_seq(seq, chrom)
                chrom = line.strip().split()[0][1:]
                seq = ""
            else:
                seq += line.strip()
        
        if seq and chrom:
            process_seq(seq, chrom)

def main():
    parser = argparse.ArgumentParser(description="Extract homopolymer regions from a reference genome and output them in BED format.")
    parser.add_argument('reference_file', type=str, help='Input reference genome file in FASTA format')
    parser.add_argument('min_length', type=int, help='Minimum length of homopolymer regions to extract')
    parser.add_argument('output_file', type=str, help='Output BED file to store the homopolymer regions')

    args = parser.parse_args()

    extract_homopolymers(args.reference_file, args.min_length, args.output_file)

if __name__ == "__main__":
    main()