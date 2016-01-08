from collections import defaultdict
import matplotlib.pyplot as plt
import re
import scipy.stats
import subprocess
import os
import cPickle
import mod_utils
import numpy as np
from collections import Counter

import pysam
import bzUtils

class TPS_Lib:
    def __init__(self, experiment_settings, lib_settings):
        """
        Constructor for Library class
        """
        self.experiment_settings = experiment_settings
        self.lib_settings = lib_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir

        self.count_reads()

        self.enrichment_sorted_mappings = None


    def count_reads(self):
        """
        reads in own mapped read SAM file and computes for each rRNA:
        1)read 5' ends mapping to each position
        2)reads overlapping each position
        3)mismatches at each position

        4)Total read 5' ends at each nt type
        5)total reads overlapping at each nt type
        6)Counts for every mutation type
        """