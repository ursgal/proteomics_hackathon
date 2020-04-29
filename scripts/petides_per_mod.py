import csv
import sys
from collections import Counter
import matplotlib.pyplot as plt


def main(all_merged):

    with open(all_merged) as fin:
        counter = Counter()
        reader = csv.DictReader(fin)
        for i, line in enumerate(reader):
            # if i > 10_000:
            #     break
            if line['Modifications'] == '':
                continue
            mods = line['Modifications'].split(';')
            mods = set([m.rsplit(':', maxsplit=1)[0] for m in mods])
            for m in mods:
                m = m.strip()
                counter[m] += 1
    names, counts = zip(*counter.most_common(50))
    # plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.bar(names, counts)
    plt.title('Peptides per modification')
    plt.xlabel("Modifications")
    plt.xticks(rotation=90)
    plt.ylabel("Count")
    plt.show()
    plt.savefig('peptides_per_mod.png')


if __name__ == '__main__':
    main(sys.argv[1])
