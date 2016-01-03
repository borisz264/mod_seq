__author__ = 'boris'
"""
I really want to use this mod-seq data to generate something resembling an x-ray exposure of a sequencing gel

inputs:
    chromosome: the chromosome to plot from
    start: the position to start plotting from
    stop: the position to stop plotting
    mutations.pkl: from parse_shapemapper_counts2.py, pickled dict of [chromosome][dataset][position] = mutation counts
    coverage.pkl: from parse_shapemapper_counts2.py, pickled dict of [chromosome][dataset][position] = coverage counts
    datasets: any number of names matching datsets in the mutations and coverage pickles

"""

import sys, mod_utils, math
import numpy
import scipy.stats as stats
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines

def gaussian(x, mean, variance, scale):
    return (float(scale)/math.sqrt(variance*2*math.pi)) * math.exp(-1*((float(x)-float(mean))**2)/(2*variance))


def generate_mutation_rates(chromosome, start, stop, mutations, coverage, dataset_names):
    mutation_rates = {}
    for dataset_name in dataset_names:

        mutation_array = [float(mutations[chromosome][dataset_name][position]) if position in mutations[chromosome][dataset_name] else 0.0 for position in range(start, stop+1)]
        coverage_array = [float(coverage[chromosome][dataset_name][position]) if position in coverage[chromosome][dataset_name] else 0.0 for position in range(start, stop+1)]
        assert len(mutation_array) == len(coverage_array)
        mutation_rates[dataset_name] = numpy.array(mutation_array)/numpy.array(coverage_array)
    return mutation_rates

def generate_gaussian_densities(mutation_rates, band_spacing, band_variance):
    """

    :param mutation_rates: output of generate_mutation_rates
    :param band_spacing: number of positions (pixels) between band centers
    :param band_variance: the variance of the guassian representing each band (in positions/pixels)
    :return:
    """
    band_densities = {}
    for dataset in mutation_rates:
        #center of first band will be at band_spacing/2 (zero-indexed, so wiht a spacing of 20, it will be position 10,
        #  the 11th entry in the array)
        #center of each subsequent band will be band_spacing+1 units from the previous 1
        summed_density_array = numpy.array([0.0]*(len(mutation_rates[dataset])*(band_spacing+1)))
        for position in range(len(mutation_rates[dataset])):
            band_position = (band_spacing/2) + (band_spacing+1)*position
            band_intensity = mutation_rates[dataset][position]
            #compute the contribution, at each position in the area of interest, from this band
            temp_array = numpy.array([gaussian(x, band_position, band_variance, band_intensity) for x in range(len(summed_density_array))])
            #temp_array = numpy.array([gaussian(x, band_position, math.sqrt(band_intensity)*band_variance, band_intensity) for x in range(len(summed_density_array))])

            summed_density_array = summed_density_array + temp_array
        band_densities[dataset] = summed_density_array
    return band_densities

def plot_density_lines(dataset_names, gaussian_densities, band_spacing, start, stop, chromosome):
    fig = plt.figure()
    plot = fig.add_subplot(111)
    color_index = 0
    for dataset in dataset_names:
        plot.plot(gaussian_densities[dataset], color = mod_utils.colors[color_index], label= dataset)
        color_index +=1
    lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlabel('nt position in %s' % chromosome)
    plot.set_ylabel('band intensity')
    plot.set_xticks([(band_spacing/2) + (band_spacing+1)*position for position in range(1+stop-start)])
    plot.set_xticklabels(range(start, stop+1))
    plt.xticks(rotation=90)
    plt.show()


def main():
    '''
    chromosome, start, stop, mutations, coverage = sys.argv[1:6]
    start = int(start)
    stop = int(stop)
    dataset_names = sys.argv[6:]
    '''

    chromosome, start, stop, mutations, coverage = 'S.c.25S__rRNA', 2740, 2800, '/Users/boris/Desktop/20151129_combined_counts/20151129_mutations.pkl',\
                                                   '/Users/boris/Desktop/20151129_combined_counts/20151129_depths.pkl'

    dataset_names = ['80S_DMSO_DMSO', '80S_DMSO_20mMDMS', '80S_CHX_20mMDMS' ]

    mutations = mod_utils.unPickle(mutations)
    coverage = mod_utils.unPickle(coverage)
    mutation_rates = generate_mutation_rates(chromosome, start, stop, mutations, coverage, dataset_names)
    band_spacing = 200
    #band_variance = 50000.
    band_variance = 2000.
    gaussian_densities = generate_gaussian_densities(mutation_rates, band_spacing, band_variance )
    plot_density_lines(dataset_names, gaussian_densities, band_spacing, start, stop, chromosome)
main()

