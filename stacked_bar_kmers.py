import matplotlib.pyplot as plt
import collections
import numpy
import math
import aColors
import tps_utils

def plot_stack(sp, bin_edges, sequence_mappings_passing_cutoff, protein, scale='linear', thresh='2sig'):
    assert scale in ['linear', 'log']
    sorted_mappings = sorted(sequence_mappings_passing_cutoff)
    enrichments = [sequence_mapping.enrichment for sequence_mapping in sorted_mappings]
    if thresh == '2sig':
        thresh = numpy.mean(enrichments) + 2 * numpy.std(enrichments)
    bottoms = []
    for bin_left, bin_right in tps_utils.pairwise(bin_edges):
        color2height = collections.Counter()
        width = bin_right - bin_left
        for sequence_mapping, enrichment in zip(sorted_mappings, enrichments):
            if enrichment > bin_right or enrichment < bin_left:
                continue
            color = aColors.protein_colors(sequence_mapping, protein, enrichment >= thresh)
            color2height[color] += 1
        if scale == 'linear':
            bottom = 0.
            for color in aColors.ordered_colors:
                if not color2height[color]:
                    continue
                sp.bar(bin_left, color2height[color], width=width, facecolor=color, edgecolor='none', bottom=bottom)
                bottom += color2height[color]
        elif scale == 'log':
            bottom = 0.
            for color in aColors.ordered_colors:
                if not color2height[color]:
                    continue
                sp.bar(bin_left, math.log(color2height[color] + bottom + 1, 2) - math.log(bottom + 1, 2),
                  width=width, facecolor=color,
                  edgecolor='none', bottom=math.log(bottom + 1, 2))
                bottom += color2height[color]
        bottoms.append(bottom)
    return max(bottoms)
