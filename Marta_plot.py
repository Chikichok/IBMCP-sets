import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.colors as colors
import matplotlib.markers as markers
import matplotlib.pyplot as plt
import numpy as np

ROOT_DIR = Path('/home/carlos/Marta_efficiency_plot_2019-6-2')
GENOTYPE_FPATH = ROOT_DIR / "T1_genotype"
SUMMARY_FNAME = "summary.csv"


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


def obtain_marker_and_color_combinations_for_plot():
	markers_series = []
	for marker_symbol in markers.MarkerStyle.markers.keys():
		if markers.MarkerStyle.markers[marker_symbol] == 'nothing':
			continue
		markers_series.append(marker_symbol)
	colors_series = []
	for color_symbol in colors.BASE_COLORS.keys():
		colors_series.append(color_symbol)
	marker_and_color_combinations = []
	for marker in markers_series:
		for color in colors_series:
			if color == 'w':
				continue
			marker_and_color_combination = [marker, color]
			marker_and_color_combinations.append(marker_and_color_combination)
	return marker_and_color_combinations


def plot_data(overall_data):
	marker_and_color_combinations = obtain_marker_and_color_combinations_for_plot()
	plot_counter = 0
	plt.figure(figsize=(20, 10))
	plt.axes([0.3, 0.3, .55, .55])
	for sample, data in overall_data.items():
		names = list(data.keys())
		numbers = list(data.values())
		color=marker_and_color_combinations[plot_counter][1]
		marker=marker_and_color_combinations[plot_counter][0]
		plt.plot(names, numbers, color=color, marker=marker, label=sample) 
		plot_counter += 1
		if plot_counter == 30:
			break
	plt.legend(loc='upper right', bbox_to_anchor=(-0.1, 1), fontsize=10, ncol=3)
	plt.xlabel('Gene and Guide')
	plt.ylabel('ICE value')
	plt.savefig('plot.svg', dpi=300, format='svg')



def main():
	data_fpaths = list_efficiencies_data_fpaths(GENOTYPE_FPATH)
	data_len = len(data_fpaths)
	overall_data = arrange_efficiencies_data(data_fpaths)
	plot_data(overall_data)


if __name__ == "__main__":
	main()