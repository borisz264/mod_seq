import os, sys, subprocess

folder, threads = sys.argv[1:]
for file_name in os.listdir(folder):
    full_path = os.path.join(folder, file_name)
    print full_path
    command_to_run = 'python mod_seq_main.py %s --threads %s' % (full_path, threads)
    subprocess.Popen(command_to_run, shell=True).wait()