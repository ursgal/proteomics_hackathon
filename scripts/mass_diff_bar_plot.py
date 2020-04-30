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
            if line['Mass Difference'] == '':
                continue
            if len(line['Mass Difference'].split(':'))>1:
                if line['Mass Difference'].split(':')[1] == 'n':
                    continue
            md = line['Mass Difference'].split(';')
            md = set([m.rsplit(':', maxsplit=1)[0] for m in md])
            # print(md)
            for _md in md:
                # _md = _md.replace('(', '').replace(')', '')
                if '(' in _md:
                    _md = _md.split('(')[0]
                counter[round(float(_md), 5)] += 1
    names, counts = zip(*counter.most_common(50))
    names = [str(n) for n in names]
    print(names, counts)
    # plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.bar(names, counts)
    plt.title('Binned mass differences')
    plt.xlabel("delta mass")
    plt.xticks(rotation=90)
    plt.ylabel("Count")
    plt.show()
    plt.savefig('mass_diff_bar_plot.png')


if __name__ == '__main__':
    main(sys.argv[1])
