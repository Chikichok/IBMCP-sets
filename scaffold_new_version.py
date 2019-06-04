import csv
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib.colors as colors
import matplotlib.markers as markers
import matplotlib.pyplot as plt
import numpy as np



COLORS = ['#006666', '#ff1a66', '#66ff33', '#ff33ff', '#996600',
          '#00ffff', '#76a2a2', '#ffff1a', '#ff8000', '#999966']

GENOTYPE_FPATH = Path('/home/marta/IBCMP/marta/T1_genotype_all_genes/genes')
SUMMARY_FNAME = "summary.csv"
BRIEF_SUMMARY = "/home/marta/IBCMP/marta/t1_resumen.csv"


def list_efficiencies_data_fpaths(data_dir, summary_fname=SUMMARY_FNAME):
    results_fpaths = []
    for gene in data_dir.iterdir():
        summary_fpath = gene / summary_fname
        if summary_fpath.exists():
            results_fpaths.append(summary_fpath)
        else:
            raise RuntimeError("file {} does not exist".format(str(summary_fpath)))
    return results_fpaths


def _get_gene_and_sample_from_label(label):
    label_items = label.split('-')
    gene_id = label_items[0]
    sample = '-'.join(label_items[1:])
    return gene_id, sample

def abreviate_name(gene_id):
    gene_id = gene_id.replace('SPL', '')
    return(gene_id)


def arrange_efficiencies_data(results_fpaths, abreviate_names=True):
    overall_data = defaultdict(dict)
    labels = []
    for results_fpath in results_fpaths:
        results_fhand = open(results_fpath)
        for result in csv.DictReader(results_fhand):
            gene_id, sample = _get_gene_and_sample_from_label(result['Label'])
            if abreviate_names:
                gene_id = abreviate_name(gene_id)
            if result['ICE'] != 'None':
                    overall_data[sample][gene_id] = int(result['ICE'])
    return overall_data


def obtain_marker_and_color_combinations_for_plot(use_hex_colors=False):
    markers_series = []
    for marker_symbol in markers.MarkerStyle.markers.keys():
        if markers.MarkerStyle.markers[marker_symbol] == 'nothing':
            continue
        markers_series.append(marker_symbol)
    colors_series = []
    if use_hex_colors:
        colors_series = COLORS
    else:
        for color_symbol in colors.BASE_COLORS.keys():
            colors_series.append(color_symbol)
    marker_and_color_combinations = []
    for marker in markers_series:
        for color in colors_series:
            if color == 'w':
                continue
            marker_and_color_combination = [marker, color]
            marker_and_color_combinations.append(marker_and_color_combination)
    print(marker_and_color_combinations)
    return marker_and_color_combinations


def get_selected_features(fhand):
    values = Counter()
    selected_genes = []
    for sample in csv.DictReader(fhand, delimiter="\t"):
        if not selected_genes:
            selected_genes = [abreviate_name(gene) for gene in sample.keys() if gene != "Sample"]
        for key, eff in sample.items():
            try:
                eff = int(eff)
                if eff >= 90:
                    values[sample["Sample"]] += 1
            except:
                continue
    return selected_genes,  [sample[0] for sample in values.most_common(10)]


def plot_data(overall_data, use_hex_colors=False):
    marker_and_color_combinations = obtain_marker_and_color_combinations_for_plot(use_hex_colors=use_hex_colors)
    plot_counter = 0
    plt.figure(figsize=(20, 10))
    plt.axes([0.3, 0.3, .55, .55])
    for sample, data in overall_data.items():
        names = list(data.keys())
        numbers = list(data.values())
        color = marker_and_color_combinations[plot_counter][1]
        marker = marker_and_color_combinations[plot_counter][0]
        plt.plot(names, numbers, color=color, marker=marker, label=sample) 
        plot_counter += 1
        #if plot_counter == 10:
        #    break
    plt.legend(loc='upper right', bbox_to_anchor=(-0.1, 1), fontsize=10, ncol=3)
    plt.xlabel('Gene and Guide')
    plt.ylabel('ICE value')
    plt.savefig('plot.svg', dpi=300, format='svg')

def get_filtered_data(fhand, selected_samples):
    filtered_data = defaultdict(dict)
    for sample in csv.DictReader(fhand, delimiter="\t"):
        if sample["Sample"] not in selected_samples:
            continue
        for key, value in sample.items():
            if key == "Sample":
                continue
            if not value:
                filtered_data[sample['Sample']][key] = None
            else:
                filtered_data[sample['Sample']][key] = int(value)
                
    return filtered_data



def main():
    data_fpaths = list_efficiencies_data_fpaths(GENOTYPE_FPATH)
    overall_data = arrange_efficiencies_data(data_fpaths)
    #brief_summary_fhand = open(BRIEF_SUMMARY)
    #selected_genes, selected_samples = get_selected_features(brief_summary_fhand)
    #brief_summary_fhand.seek(0)
    #filtered_data = get_filtered_data(brief_summary_fhand, selected_samples)
    plot_data(overall_data)


if __name__ == "__main__":
    main()