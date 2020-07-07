from phylogeny import Phylogeny
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import argparse

'''
Arguments:
    filename: name of fasta file to read
Returns:
    sequences: dictionary of outputs (string (sequence id) -> sequence (string))
    size: length of each sequence
'''
def read_data(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        sequences = {}
        output = ''
        size = 0
        curr = ''
        flag = False

        for line in lines[1:]:
            l = line.split()
            sequences[l[0]] = l[1]
            size = len(l[1])
        
    return sequences, size


def main():
    # Define phylogeny.
    rates = {
        'r1': {'rate': 0.2, 'prob': np.log(0.30)},
        'r2': {'rate': 0.3, 'prob': np.log(0.30)},
        'r3': {'rate': 0.4, 'prob': np.log(0.40)}
    }
    # rates = {
    #     'Low': {'rate': 0.2, 'prob': np.log(0.333)},
    #     'Medium':{'rate': 0.3, 'prob': np.log(0.333)},
    #     'High': {'rate': 0.4, 'prob': np.log(0.333)}
    # }

    # rates = {
    #     'Medium':{'rate': 0.3, 'prob': np.log(1)}
    # }

    nodes = ['human',
        'gorilla',
        'chimp',
        'gibbon',
        'golden',
        'orangutan',
        'green',
        'root',
        'interim1',
        'interim2',
        'interim3',
        'interim4',
        'interim5']

    tree = {
        'human': [],
        'gorilla': [],
        'chimp': [],
        'gibbon': [],
        'golden':[],
        'orangutan':[],
        'green':[],
        'root':['interim1', 'interim2'],
        'interim1': ['human', 'interim5'],
        'interim2': ['interim3', 'interim4'],
        'interim3': ['green', 'gorilla'],
        'interim4': ['golden', 'gibbon'],
        'interim5': ['chimp', 'orangutan']
    }

    root = 'root'

    data, seqlen = read_data('../data/dna_data.txt')
    # seqlen = 500

    branch_lengths = {
        'human': {},
        'gorilla': {},
        'chimp': {},
        'gibbon': {},
        'golden': {},
        'orangutan':{},
        'green': {},
        'root':{'interim1': 1.197025, 'interim2': 1.197025},
        'interim1': {'human': 0.00006, 'interim5': 2.87046},
        'interim2': {'interim3': 0.12230, 'interim4': 1.82491},
        'interim3': {'green': 3.09816, 'gorilla': 3.81994},
        'interim4': {'golden': 4.93107, 'gibbon': 1.22929},
        'interim5': {'chimp': 4.14481, 'orangutan': 3.75325}
    }

    # Compute likelihood and most likely rates.
    phylo = Phylogeny(rates, nodes, tree, root, data, seqlen, branch_lengths)
    rate_list, ll = phylo.likelihood()

    print('log likelihood=' + str(ll))
    print('most likely sitewise rates=' + str(rate_list))

    # Collect data and format for plotting.
    evo_rates = []
    seq_diffs = []
    rate_range = {rate:[] for rate in rates}

    prev_elem = rate_list[0]
    range_start = 0
    for site in range(seqlen):
        if not (rate_list[site] == prev_elem):
            rate_range[prev_elem].append((range_start, site))
            prev_elem = rate_list[site]
            range_start = site

        evo_rates.append(rates[rate_list[site]]['rate'])

        diffs = 0
        for species1 in data:
            for species2 in data:
                if data[species1][site] != data[species2][site]:
                    diffs += 1
        seq_diffs.append(diffs)
    rate_range[prev_elem].append((range_start, seqlen - 1))

    print(rate_range)

    # Plot results.
    fig, ax = plt.subplots()

    colormap = cm.get_cmap('Spectral')
    rate_idx = 0
    for rate in rates:
        range_idx = 0
        color = colormap(rate_idx/len(rates))
        for start, end in rate_range[rate]:
            if range_idx == 0:
                ax.axvspan(start, end, alpha=0.5, color=color, label=rate)
            else:
                ax.axvspan(start, end, alpha=0.5, color=color)
            range_idx += 1
        rate_idx += 1

    # ax.plot(range(seqlen), evo_rates)
    ax.plot(range(seqlen), seq_diffs, label='Sequence Differences')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()