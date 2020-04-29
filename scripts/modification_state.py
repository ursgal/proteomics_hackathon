import csv
import sys
import matplotlib.pyplot as plt


def main(all_merged):

    with open(all_merged) as fin:
        reader = csv.DictReader(fin)
        counter = {
            'modified': 0,
            'unmodified': 0,
            'MassShift': 0,
            'Glyco': 0,
        }
        for i, line in enumerate(reader):
            # if i > 10_000:
            #     break
            if line['Modifications'] == '':
                counter['unmodified'] += 1
            if line['Modifications'] != '':
                counter['modified'] += 1
            if line['Mass Difference'] != '':
                counter['MassShift'] += 1
            if line['Glycan'] != '':
                counter['Glyco'] += 1

    names, counts = zip(*counter.items())
    print(counter)
    # plt.tight_layout()
    # plt.gcf().subplots_adjust(bottom=0.20)
    plt.bar(names, counts)
    plt.xlabel("Modification state")
    # plt.xticks(rotation=90)
    plt.ylabel("Count")
    plt.show()
    # plt.savefig('modification_state.png')


if __name__ == '__main__':
    main(sys.argv[1])
