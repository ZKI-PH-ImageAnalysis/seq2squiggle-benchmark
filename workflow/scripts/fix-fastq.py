def modify_fastq_headers(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('@'):
                header_parts = line.strip().split('_')
                new_header = header_parts[0] + '\n'  # Keep everything before the first underscore
                f_out.write(new_header)
            else:
                f_out.write(line)

# Usage
input_file = 'results/zymo-human/passed.fastq'
output_file = 'results/zymo-human/passed.fastq'
modify_fastq_headers(input_file, output_file)
