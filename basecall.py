import numpy as np
import pyarrow as pa
from pyarrow import csv
import shutil
import os
import subprocess
import concurrent
import concurrent.futures
from pathlib import Path
from pod5.tools import pod5_view
from pod5.tools import pod5_subset

### I think this function can be removed, since I'm not going to use duplex for right now, and this function is only useful for duplex.
# def pod5_split(path_out, path_pod5, python_env, account, mail, mail_type = 'None', workers=5) :
#     """
#     Runs slurm jobs to split .pod5 files into separate folders based on sequencing channel. This is needed for running dorado duplex on the OSC with multiple workers. pod5 files are split and saved under path_out/split_pod5/{channel}, and the job scripts and .out logs are saved under path_out/logs_split_pod5. Currently requires a specified python environment with the pod5 package from ONT installed.
    
#     Args:
#         path_out (str) : directory to be used as the output folder. This will be created if it doesn't exist.
#         path_pod5 (str) : directory containing .pod5 files. The function will recursively find the files in this folder.
#         python_env (str) : name of your conda environment with pod5 installed.
#         account (str) : OSC account to be billed for the compute time.
#         mail (str) : email for any status messages from the jobs.
#         mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
#         workers (int) : number of jobs to run for the pod5 splitting. Defaults to 5. Rough recommendation is to do at most 2M reads per worker ie use >= 10 workers for 20M reads. Doing more just runs the risk of memory and time issues, and the OSC is hardly ever out of CPU space.
#     """
#     Path(path_out).mkdir(parents=True, exist_ok=True)
#     if not os.path.isfile(Path(path_out + "/view.txt")) :
#         print("creating view file")
#         pod5_view.view_pod5([Path(path_pod5)], Path(path_out + '/view.txt'), include = "read_id, channel", force_overwrite=True, threads=workers)
#     view = csv.read_csv(path_out + "/view.txt", parse_options=csv.ParseOptions(delimiter='\t'))
#     channels = view.column('channel').unique().to_pylist()
#     if Path(path_out + '/split_pod5/').is_dir() :
#         shutil.rmtree(Path(path_out + '/split_pod5/'))
#     if Path(path_out + '/logs_split_pod5/').is_dir() :
#         shutil.rmtree(Path(path_out + '/logs_split_pod5/'))
#     Path(path_out + '/split_pod5/').mkdir(parents=True, exist_ok=True)
#     Path(path_out + '/logs_split_pod5/').mkdir(parents=True, exist_ok=True)
#     for channel in channels :
#         Path(path_out + '/split_pod5/' + str(channel)).mkdir(parents=True, exist_ok=True)
#     files = [str(x) for x in Path(path_pod5).iterdir()]
#     file_chunks = np.array_split(files, workers )
#     print("done with prep work")
#     for i in range(workers) :
#         script = [
#             "#!/bin/bash\n",
#             "#SBATCH --account=" + account + "\n",
#             "#SBATCH --job-name=pod5_split_" + str(i) + "\n",
#             "#SBATCH --nodes=1 --ntasks-per-node=1\n",
#             "#SBATCH --cpus-per-task=8\n",
#             "#SBATCH --output=" + path_out + "/logs_split_pod5/" + str(i) + ".out" + "\n",
#             "#SBATCH --mail-type=" + mail_type + "\n",
#             "#SBATCH --time=4:00:00\n",
#             "#SBATCH --mail-user=" + mail + "\n",
#             "module load miniconda3/24.1.2-py310\n",
#             "conda activate " + python_env + "\n",
#             "pod5 subset " + ' '.join(file_chunks[i]) + " --table " + path_out + "/view.txt -o " + path_out + "/split_pod5/ --threads 16 --template '{channel}/{channel}_" + str(i) + ".pod5' -M --columns channel"
#         ]
#         with open(str(path_out + "/logs_split_pod5/" + str(i) + ".sh"), 'w') as handle:
#             handle.writelines(script)
#         args = ["sbatch", str(path_out + "/logs_split_pod5/" + str(i) + ".sh")]
#         subprocess.run(args)
#     return

# def check_pod5_jobs(user) :
#     """
#     Checks squeue to see if any pod5_split jobs are running.
    
#     Args :
#         user (str) : the OSC username used for the pod5_split jobs.
#     Returns :
#         bool, True if any jobs are found.
#     """
#     args = ['squeue', '-u', user, '-v', '-r', '--states=BF,CF,CG,DL,F,NF,OOM,PD,PR,R,RD,RS,RV,ST,S,TO', '--format="%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R"']
#     squeue = subprocess.run(args, capture_output=True)
#     stdout = str(squeue.stdout)
#     squeue_nl_split = stdout.split('\\n')
#     for item in squeue_nl_split :
#         if 'pod5_split' in item :
#             return True
#     return False

def move_pod5s(path_pod5, path_out, workers = 1, copy_files = True, overwrite = False) :
    """
    Moves (or copies) pod5 files or folders of pod5s into grouped folders for parallel running of Dorado on an HPC. The resulting folders are put under path_out.
        Can be used following pod5_split() if running Dorado duplex, or you can just point path_pod5 to a folder containing your original pod5 files, if not running duplex.
        This will also remove the original path_pod5 directory, unless copy_files=True.
    
    Args:
        path_pod5 (str) : directory containing pod5 files, or containing the folders output by pod5_split()
        path_out (str) : directory to be used as the output folder.
        workers (int) : number of folder to split the pod5s into. Should match the number of dorado jobs that you want to run.
        copy_files (bool) : whether or not to copy the data instead of moving it. Moving is typically MUCH faster, but this is set to False for safety.
        overwrite (bool) : whether or not to delete any contents of path_out.
    """
    pod5_items = [x for x in Path(path_pod5).iterdir()]
    if len(Path(path_out).iterdir()) > 0 :
        if overwrite :
            shutil.rmtree(Path(path_out))
            Path(path_out).mkdir(parents=True, exist_ok=True)
        else :
            raise FileExistsError("Error: there are already contents in ", path_out, ". Set overwrite=True if you'd like to delete the existing contents.")
    pod5_item_chuncks = np.array_split( pod5_items, workers )
    for i in range(len(pod5_item_chuncks)) :
        print("started moving group ", i)
        Path( path_out + '/' + str(i) ).mkdir(parents=True, exist_ok=False)
        for pod5_item in pod5_item_chuncks[i] :
            if copy_files :
                shutil.copy(pod5_item, Path( path_out + '/' + str(i) + '/' ))
            else :
                shutil.move(pod5_item, Path( path_out + '/'  + str(i) + '/' ))
        print("finished moving group ", i)
    if not copy_files :
        shutil.rmtree(path_pod5)
    return

def run_dorado_job(path_pod5, path_out, path_logs, path_dorado, account, mail, mail_type = 'None', i = 0) :
    """
    This is a utility script that runs a slurm job for dorado basecaller, sending the output sam file to path_out. Creates a script file in path_logs, which is also where the .out file from the job will be.
        Currently is hardcoded to request 1 GPU, 2 CPU cores, 40G memory, and 12 hours, and to run basecaller with the sup model, emit-sam, and no trimming. The script settings below can be edited as needed for your purposes.
    
    Args :
        path_pod5 (str) : directory containing pod5 files.
        path_out (str) : directory to be used as the output folder for the resulting sam file.
        path_logs (str) : directory to hold the .sh and .out files for the slurm job.
        path_dorado (str) : path to the dorado executable.
        account (str) : OSC account to be billed for the compute time.
        mail (str) : email for any status messages from the jobs.
        mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
        i (int) : the worker_ID to be used in the title of the job script, log, and sam output file.
    """
    script = [
        "#!/bin/bash\n",
        "#SBATCH --account=" + account + "\n",
        "#SBATCH --job-name=dorado_" + str(i) + "\n",
        "#SBATCH --nodes=1\n",
        "#SBATCH --ntasks=1\n",
        "#SBATCH --cpus-per-task=2\n",
        "#SBATCH --mem=40G\n",
        "#SBATCH --gres=gpu:1\n",
        "#SBATCH --time=12:00:00\n",
        "#SBATCH --output=" + path_logs + "/dorado_" + str(i) + ".out" + "\n",
        "#SBATCH --mail-type=" + mail_type + "\n",
        "#SBATCH --mail-user=" + mail + "\n",
        str(path_dorado) + " basecaller sup --emit-sam -r -v --no-trim " + str(path_pod5) + '/' + str(i) + "/ > " + str(path_out) + '/' + str(i) + ".sam"
    ]
    with open(path_logs + "/dorado_" + str(i) + '.sh', 'w') as handle:
        handle.writelines(script)
    args = ['sbatch', str(path_logs + '/dorado_' + str(i) + '.sh')]
    subprocess.run(args)
    return

def dorado_slurm(path_pod5, path_out, path_logs, path_dorado, account, mail, mail_type = 'None', python_env = None, user = None, workers = 1) :
    """
    Run Dorado basecalling through slurm jobs. Can use move_pod5s() to move/copy the split pod5s into a folder for each dorado worker under path_pod5_folders. 
        Then runs parallel dorado jobs requesting GPUs. Note that Dorado will only run on NVIDIA GPUs with Tensor cores, which are only in the Voltair line and newer.
        The resulting sam files are put under path_out.
    
    Args :
        path_pod5 (str) : the output folder of move_pod5s. This should contain folders of pod5 files.
        path_out (str) : directory to be used as the output folder for the resulting sam files.
        path_logs (str) : directory to hold the .sh and .out files for the slurm job.
        path_dorado (str) : path to the dorado executable.
        account (str) : OSC account to be billed for the compute time.
        mail (str) : email for any status messages from the jobs.
        mail_type (str) : type of messages to recieve for the jobs ie 'Start', 'All', 'Error'. Defaults to 'None'
        workers (int) : number of jobs to run for dorado. If you already ran move_pod5s, this should be the SAME number of workers as used for that. You'll either leave some reads un-basecalled or run unnecessary jobs if not.
    """
    if len(Path(path_pod5).iterdir()) == 0 :
        raise FileNotFoundError("Error: no contents found in path_pod5.")
    if len(Path(path_out).iterdir()) > 0 :
        raise FileExistsError("Error: the path_out directory has data in it. Please check the directory and remove unwanted data.")
    Path(path_out).mkdir(parents=True, exist_ok=True)
    Path(path_logs).mkdir(parents=True, exist_ok=True)
    for i in range(workers) :
        run_dorado_job(path_pod5, path_out, path_logs, path_dorado, account, mail, mail_type, i)
    return

def check_dorado_jobs(path_logs) :
    """
    Checks the .out logs in path_logs to make sure that all dorado jobs ran successfully. Checks for any lines in the .out with 'Basecalled' in them, which appears to only be present if dorado has completed basecalling.
    
    Args :
        path_logs (str) : directory used as the path_logs folder for dorado.
    """
    dorado_logs = [ x for x in Path(path_logs).iterdir() if x.suffix == '.out' ]
    broken_jobs = []
    for dorado_log in dorado_logs :
        broken = True
        with open(dorado_log, 'r') as handle :
            for line in handle.readlines() :
                if 'Basecalled' in line :
                    broken = False
                    break
        if broken :
            name = dorado_log.stem
            i = name.split('_')[1]
            broken_jobs.append(i)
    if len(broken_jobs) != 0 :
        print("Broken jobs were found : ", broken_jobs)
    else : 
        print("No broken jobs found! :D")
    return
