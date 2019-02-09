import argparse
import random
import string
import glob
from Bio import SeqIO
import itertools
import os
import time
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description=('Align contigs in a job array'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--genome_dir',
                        default=None,
                        type=str,
                        help='Directory containing genome files.')
    parser.add_argument('--genome_ext',
                        default=None,
                        type=str,
                        help='Extension for genome files (e.g., .fasta')
    parser.add_argument('--temp_script_dir',
                        default=None,
                        type=str,
                        help='Temporary directory for run scripts.')
    parser.add_argument('--slurm_output_dir',
                        default=None,
                        type=str,
                        help='Directory for slurm output files. Useful for debugging.')
    parser.add_argument('--alignment_dir',
                        default=None,
                        type=str,
                        help='Directory for alignments.')
    parser.add_argument('--cluster_username',
                        default=None,
                        type=str,
                        help='Username on cluster.')
    parser.add_argument('--partition_name',
                        default=None,
                        type=str,
                        help='Partition to use on cluster.')
    parser.add_argument('--mugsy_path',
                        default=None,
                        type=str,
                        help='Path to mugsy.')
    parser.add_argument('--utility_path',
                        default=None,
                        type=str,
                        help='Path to utilities.')
    parser.add_argument('--base_name',
                        default=None,
                        type=str,
                        help='base output file name')
    parser.add_argument('--max_jobs',
                        default=None,
                        type=int,
                        help='Maximum number of jobs to run simultaneously.')
    parser.add_argument('--max_submit',
                        default=None,
                        type=int,
                        help='Maximum number of jobs to submit simultaneously.')
    parser.add_argument('--final_output_dir',
                        default='./',
                        type=str,
                        help='Directory for final output.')

    args = parser.parse_args()

    check_inputs(args)

    contig_file_list = glob.glob('{contigdir}/*{extension}'.format(contigdir=args.genome_dir,
                                                                   extension=args.genome_ext))
    print('Renaming genome files')
    renamed_contig_files = rename_files(contig_file_list)

    print('Making shell files.')
    final_simfile_list, final_scriptfile_list = make_shell_files(renamed_contig_files,
                                                                 args.utility_path,
                                                                 args.alignment_dir,
                                                                 args.mugsy_path,
                                                                 args.temp_script_dir)

    print('Aligning contigs')
    align_array(args.temp_script_dir,
                args.slurm_output_dir,
                args.cluster_username,
                args.partition_name,
                final_scriptfile_list,
                maxjobs=args.max_jobs,
                jobsubmit=args.max_submit)

    missing_simfiles = []
    final_lines = []
    for simfile in final_simfile_list:
        try:
            with open(simfile, 'r') as simfile:
                final_lines.append(next(simfile))
        except FileNotFoundError:
            missing_simfiles.append(simfile)

    final_outfile_name = '{final_output_dir}/{base_name}_cat_ssd.txt'.format(base_name=args.base_name,
                                                                             final_output_dir=args.final_output_dir)
    with open(final_outfile_name, 'w') as outfile:
        outfile.write('\t'.join(['Strain 1',
                                 'Strain 2',
                                 'Initial divergence',
                                 'Alignment size',
                                 'Genome 1 size',
                                 'Genome 2 size',
                                 'Observed SSD',
                                 'SSD 95 CI low',
                                 'SSD 95 CI high\n']))
        outfile.writelines(final_lines)
    if missing_simfiles != []:
        print('WARNING: {num} files were expected but not generated.'.format(num=str(len(missing_simfiles))))
        print(missing_simfiles)


def check_inputs(args):
    # Check for the existence of the genome directory
    os.listdir(args.genome_dir)

    # Check for the existence of the utility path directory
    os.listdir(args.utility_path)

    # Check that contig files exist in the directory
    contig_list = glob.glob('{contigdir}/*{extension}'.format(contigdir=args.genome_dir,
                                                              extension=args.genome_ext))
    if len(contig_list) == 0:
        raise FileNotFoundError('Files with contig extension not found in directory.')

    # Check for temp script directory. Makes it if it isn't there.
    if not os.path.exists(args.temp_script_dir):
        print('Temporary script directory does not exist. Creating new directory.')
        os.makedirs(args.temp_script_dir)

    # Check for slurm output directory. Makes it if it isn't there.
    if not os.path.exists(args.slurm_output_dir):
        print('Slurm output directory does not exist. Creating new directory.')
        os.makedirs(args.slurm_output_dir)

    # Check for alignment directory. Makes it if it isn't there.
    if not os.path.exists(args.alignment_dir):
        print('Slurm output directory does not exist. Creating new directory.')
        os.makedirs(args.alignment_dir)

    # Checks mugsy path
    if not os.path.exists(args.mugsy_path):
        raise FileNotFoundError('Invalid mugsy path.')

    # Check utility path for scripts
    if not os.path.exists('{utility_path}/align_contig_array.py'.format(utility_path=args.utility_path)):
        raise FileNotFoundError('align_contig_array.py not in utility path')

    if not os.path.exists('{utility_path}/calculate_length_bias.py'.format(utility_path=args.utility_path)):
        raise FileNotFoundError('calculate_length_bias.py not in utility path')


def align_array(array_script_path,
                slurm_output_dir,
                username,
                partition_name,
                scriptfile_list,
                maxjobs=200,
                jobsubmit=100):
    '''
    Runs
    '''
    with open('master_align_script.sh', 'w') as outfile:
        outfile.write('#!/bin/bash\n')
        outfile.write('sh {scriptdir}/${{SLURM_ARRAY_TASK_ID}}.sh'.format(scriptdir=array_script_path))

    random_job_name = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(6))
    nums = []

    for script in sorted(scriptfile_list):
        nums.append(script.split('/')[-1].split('.')[0])
    for x in range(0, len(nums), jobsubmit):
        while not check_queue(random_job_name, username, maxjobs, jobsubmit):
            time.sleep(10)

        num_string = ','.join(sorted(nums[x: min([x + jobsubmit, len(nums)])], key=lambda x: int(x)))
        cmd = 'sbatch -N 1 -n 1 -p {partition_name} --time=20:00 -o {slurm_output_dir}/output_%A_%a.out --array={num_string} -J {jobname} master_align_script.sh'.format(
            slurm_output_dir=slurm_output_dir,
            num_string=num_string,
            jobname=random_job_name,
            partition_name=partition_name)
        os.system(cmd)


    # Waits until all jobs are finished
    while not check_queue(random_job_name, username, 0, 0):
        time.sleep(60)


def make_shell_files(contig_file_list,
                     utility_path,
                     align_directory,
                     mugsy_path,
                     array_script_path):
    all_sim_files = []
    all_script_files = []
    all_jobscripts = defaultdict(list)

    for i, contig_pair in enumerate(itertools.combinations(contig_file_list, 2)):

        contig1, contig2 = contig_pair
        strain1 = '.'.join(contig1.split('/')[-1].split('.')[0:-3])
        strain2 = '.'.join(contig2.split('/')[-1].split('.')[0:-3])
        align_name = '{strain1}_@_{strain2}.maf'.format(strain1=strain1, strain2=strain2)
        sim_name = align_name + '.sim.txt'
        all_sim_files.append(align_directory + sim_name)

        if not os.path.isfile(align_directory + '/' + align_name):
            align_command = 'python {utility_path}/align_contig_array.py {contig1} {contig2} {align_dir} {mugsy_path}\n'.format(
                             utility_path=utility_path,
                             contig1=contig1,
                             contig2=contig2,
                             align_dir=align_directory,
                             mugsy_path=mugsy_path
                             )
            all_jobscripts[i].append(align_command)

        if not os.path.isfile(align_directory + '/' + sim_name):
            simulation_command = 'python {utility_path}/calculate_length_bias.py -i {align_dir}/{infile} -o {align_dir}/{outfile} -g1 {genome1} -g2 {genome2}\n'.format(
                                  utility_path=utility_path,
                                  align_dir=align_directory,
                                  infile=align_name,
                                  outfile=sim_name,
                                  genome1=contig1,
                                  genome2=contig2
                                  )
            all_jobscripts[i].append(simulation_command)

    for script_num, joblist in all_jobscripts.items():
        outfile_name = '{array_script_path}/{num}.sh'.format(array_script_path=array_script_path,
                                                             num=str(script_num))
        with open(outfile_name, 'w') as outfile:
            outfile.writelines(joblist)
        all_script_files.append(outfile_name)

    return (all_sim_files, all_script_files)


def rename_files(genome_file_list):
    mugsy_inputs = []
    for genome in genome_file_list:

        strain_name = '.'.join(genome.split('/')[-1].split('.')[0:-1])

        outname = genome + '.renamed'
        seqs = []
        for i, s in enumerate(SeqIO.parse(genome, 'fasta')):
            s.id = '{id}_{contig_num}'.format(id=strain_name, contig_num=str(i))
            seqs.append(s)
        SeqIO.write(seqs, outname, 'fasta')

        # We want to remove all periods and colons from sequence input so that mugsy doesn't break
        mugsy_outname = genome + '.renamed.mugsy'
        mugsy_inputs.append(mugsy_outname)
        mugsy_s = []
        mugsy_name = strain_name.translate(({ord(c): '_' for c in """ !@#$%^&*()[]{};:,./<>?\|`"'~-=+"""}))
        for i, s in enumerate(SeqIO.parse(genome, 'fasta')):
            s.description = ''
            s.id = '{id}_{contig_num}'.format(id=mugsy_name, contig_num=str(i))
            mugsy_s.append(s)
        SeqIO.write(mugsy_s, mugsy_outname, 'fasta')
    return mugsy_inputs


def check_queue(jobname, username, maxjobs, jobsubmit):
    array_pending = False
    os.system('squeue -u {username} -r -n {jobname} > {jobname}.txt'.format(username=username,
                                                                            jobname=jobname))
    time.sleep(5)
    all_lines = [line for line in enumerate(open('{jobname}.txt'.format(jobname=jobname)))]
    num_jobs = len(all_lines) - 1
    if num_jobs + jobsubmit > maxjobs:
        suitable = False
    else:
        suitable = True
    os.system('rm {jobname}.txt'.format(jobname=jobname))
    return suitable


if __name__ == '__main__':
    main()
