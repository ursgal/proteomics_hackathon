import csv
import sys
from collections import Counter
import matplotlib.pyplot as plt


def main(all_merged):

    with open(all_merged) as fin:
        reader = csv.DictReader(fin)
        mass_diffs = []
        seen_pep_mass_diff = set()
        for i, line in enumerate(reader):
            # if i > 10_000:
            #     break
            if line['Mass Difference'] == '':
                continue
            if len(line['Mass Difference'].split(':'))>1:
                if line['Mass Difference'].split(':')[1] == 'n':
                    continue
            md = line['Mass Difference'].split(';')
            md = set([m.rsplit(':', maxsplit=1)[0] for m in md])
            # print(md)
            for _md in md:
                pep_mass_diff = line['Sequence'] + '#' + _md
                if pep_mass_diff in seen_pep_mass_diff:
                    continue
                else:
                    mass_diffs.append(_md)
                    seen_pep_mass_diff.add(pep_mass_diff)
    # plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.hist(mass_diffs)
    plt.title('Mass diff histo')
    plt.xlabel("mass diff")
    plt.xticks(rotation=90)
    plt.ylabel("Count")
    # plt.show()
    plt.savefig('mass_diff_histo.png')


if __name__ == '__main__':
    main(sys.argv[1])
