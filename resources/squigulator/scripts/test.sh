#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

echo "Basic DNA"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -q a.fasta -n 10 --seed 1 --dwell-std 1.0 -r 20000 -t1 || die "Running the tool failed"
diff -q test/fasta.exp a.fasta || die "diff failed"
diff -q test/slow5.exp a.slow5 || die "diff failed"

echo "Basic RNA"
ex ./squigulator -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -q a.fastq -n 10 --seed 1 --prefix=yes --dwell-std 3.0 -t1  || die "Running the tool failed"
diff -q test/rna_slow5.exp a.slow5 || die "diff failed"

echo "--ideal"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal  -r 20000 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_slow5.exp a.slow5 || die "diff failed"

echo "--ideal-time"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal-time  -r 20000 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_time_slow5.exp a.slow5 || die "diff failed"

echo "--ideal-amp"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --ideal-amp  -r 20000  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_amp_slow5.exp a.slow5 || die "diff failed"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --amp-noise 0.0  -r 20000  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_ideal_amp_slow5.exp a.slow5 || die "diff failed"


echo "--prefix=yes"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 -n 2 --seed 1 --prefix=yes  -r 20000  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_prefix_slow5.exp a.slow5 || die "diff failed"

echo "--prefix=no"
ex ./squigulator -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o a.slow5 -n 2 --seed 1 --dwell-std 3.0 -t1 || die "Running the tool failed"
diff -q test/rna_prefixno_slow5.exp a.slow5 || die "diff failed"

echo "--full-contigs"
ex ./squigulator test/nCoV-2019.reference.fasta -o a.slow5 --seed 1 --full-contigs  --dwell-std 5.0 -t1 || die "Running the tool failed"
diff -q test/dna_full_contig.exp a.slow5 || die "diff failed"

echo "r10 PAF out"
ex ./squigulator -x dna-r10-prom -o a.slow5 -n 1 --seed 1 --dwell-std 4.0 -t1 test/nCoV-2019.reference.fasta -c a.paf -q a.fa
diff -q test/dna_r10_paf.exp a.slow5 || die "diff failed"
diff -q test/dna_r10_paf.paf.exp a.paf || die "diff failed"
diff -q test/dna_r10_paf.fa.exp a.fa || die "diff failed"

echo "r9 rna paf out and sam out"
ex ./squigulator -x rna-r9-prom -o a.slow5 -n 1 --seed 1 --dwell-std 3.0 -t1 -t1 test/rnasequin_sequences_2.4.fa -c a.paf -q a.fa -a a.sam
diff -q test/rna_paf.exp a.slow5 || die "diff failed"
diff -q test/rna_paf.paf.exp a.paf || die "diff failed"
diff -q test/rna_paf.sam.exp a.sam || die "diff failed"
diff -q test/rna_paf.fa.exp a.fa || die "diff failed"

# --paf-ref and samout
echo "r10 dna paf out --paf-ref"
ex ./squigulator -x dna-r10-prom -o a.slow5 -n 2 --seed 2 --dwell-std 4.0 -t1 test/nCoV-2019.reference.fasta -c a.paf --paf-ref -a a.sam
diff -q test/dna_r10_paf-ref.exp a.slow5 || die "diff failed"
diff -q test/dna_r10_paf-ref.paf.exp a.paf || die "diff failed"
diff -q test/dna_r10_paf-ref.sam.exp a.sam || die "diff failed"

# samout only
ex ./squigulator -x dna-r10-prom -o a.slow5 -n 2 --seed 2 --dwell-std 4.0 -t1 test/nCoV-2019.reference.fasta -c a.paf -a a.sam
diff -q test/dna_r10_paf-ref.exp a.slow5 || die "diff failed"
diff -q test/dna_r10_paf-ref.sam.exp a.sam || die "diff failed"


redundancy_check () {
    N=$(grep -v ^[@#] a.slow5 | cut -f ${1}  | sort | uniq -c | sort -nr -k1,1 | head -1 | awk '{print $1}')
    [ "$N" != "1" ] && die "failed thread test for column ${1}"
}

ex ./squigulator -x dna-r9-min test/nCoV-2019.reference.fasta -n 100 -t 8 -K 10 -o a.slow5 --seed 1
# read_id        read_group      digitisation    offset  range   sampling_rate   len_raw_signal  raw_signal
# 9channel_number  10median_before   11read_number     12start_mux        13start_time
redundancy_check 1
redundancy_check 4
redundancy_check 7
redundancy_check 8
redundancy_check 10
redundancy_check 11
redundancy_check 13

echo "Test passed"