import pyslow5
import matplotlib.pyplot as plt

squi_l = None
seq2s_l = None

s5 = pyslow5.Open('squigulator_reads_fixed.slow5','r')
for read in s5.seq_reads(pA=True, aux='all'):
    if read['read_id'] == "aea835ca-08a8-4ecc-8eed-74059874cde5":
        squi_l = read['signal']

s5 = pyslow5.Open('seq2squiggle_reads.slow5','r')
for read in s5.seq_reads(pA=True, aux='all'):
    if read['read_id'] == "aea835ca-08a8-4ecc-8eed-74059874cde5":
        seq2s_l = read['signal']

print(squi_l)
import statistics
mean_squi = (sum(squi_l) / len(squi_l))
std_squi = statistics.stdev(squi_l)
print(mean_squi, std_squi)
print(len(squi_l))
#squi_l = [1/std_squi * (i - mean_squi) for i in squi_l]

print(seq2s_l)
mean_seq2s = (sum(seq2s_l) / len(seq2s_l))
std_seq2s = statistics.stdev(seq2s_l)
print(mean_seq2s, std_seq2s)
print(len(seq2s_l))
#seq2s_l = [1/std_seq2s * (i - mean_seq2s) for i in seq2s_l]

plt.plot(squi_l[:1000], label='squigulator')
plt.plot(seq2s_l[10:990], label='seq2squiggle')

# Adding labels and title
plt.xlabel('Time')
plt.ylabel('Current (pA)')

# Adding legend
plt.legend()

plt.savefig("example.png", dpi=600)