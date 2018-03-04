import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from collections import Counter
from collections import OrderedDict

sns.set(color_codes=True)

file_name = sys.argv[1]
dictionary = {}
dictionary[0] = 0
dictionary[1] = 0
dictionary[2] = 0

with open(file_name, "r") as ins:
    for line in ins:
        line = line.rstrip()
        if not line:
            continue
        line = line.split('\t')
        gt1 = int(line[7])
        dictionary[gt1] += 1
        if line[5] != '.':
            gt2 = int(line[9])
            dictionary[gt2] += 1


fig, ax = plt.subplots()

dictionary2 = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))

ax.bar(range(len(dictionary2)), dictionary2.values(), align='center')
plt.xticks(range(len(dictionary2)), ['Hom', 'Het', 'Hom-alt'])

for i, v in enumerate(dictionary2.values()):
    ax.text(i, v + 500, str(v) + ": " + str(round(v*100/sum(dictionary2.values()), 2)) + "%",
            fontweight='bold', ha='center', fontsize=8)
plt.savefig("Class_distribution_bed.png", dpi=400)