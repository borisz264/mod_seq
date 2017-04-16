# mod_seq

to run the pipeline:
1)compress all of your fastq files with gzip and place them in one folder
2)fill out a settings file, pay attention to the output location, or you risk overwriting previous runs
3)from the mod_seq folder run:
    python mod_seq_main.py my_experiment.settings.json

optional parameters (add after settings file):
    --threads : max number of parrallel processes to use. Default is 8.

Python dependencies (install with pip):
    simplejson
    bokeh (optional, but very helpful)
    statsmodels
    matplotlib
    numpy
    scipy

This package requires installation of the ShapeMapper package from Kevin Weeks' lab.
Download ShapeMapper from http://www.chem.unc.edu/rna/software.html compile it, and add it to your PATH.

note: increasing maxProc in shapemapper.py in shapemapper package improves parrallelization when run locally with multiple
cores available.

cutadapt is required: http://cutadapt.readthedocs.org/en/stable/guide.html

fastX toolkit is required: http://hannonlab.cshl.edu/fastx_toolkit/

the bokeh interactive plotting package is optionally required:
http://bokeh.pydata.org/en/latest/docs/user_guide/quickstart.html#userguide-quickstart


explanation of settings files:
since I can't figure out how to add comments to the settings files, please read settings.explanation.txt