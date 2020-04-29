import csv
import sys
from collections import Counter
import matplotlib.pyplot as plt


def main(all_merged):

    with open(all_merged) as fin:
        reader = csv.DictReader(fin)
        mass_diffs = []
        for i, line in enumerate(reader):
            # if i > 10_000:
            #     break
            if line['Mass Difference'] == '':
                continue
            md = line['Mass Difference'].split(';')
            md = set([m.rsplit(':', maxsplit=1)[0] for m in md])
            # print(md)
            for _md in md:
                mass_diffs.append(_md)
    # plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.hist(mass_diffs)
    plt.title('Mass diff histo')
    plt.xlabel("mass diff")
    plt.xticks(rotation=90)
    plt.ylabel("Count")
    plt.show()
    plt.savefig('mass_diff_histo.png')


if __name__ == '__main__':
    main(sys.argv[1])
