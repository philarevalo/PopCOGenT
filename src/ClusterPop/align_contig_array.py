import argparse
import os
import random
import string

def main():
    parser = argparse.ArgumentParser(
        description=('Makes protein files from hmmer output'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('contig1', help='path to contig 1')
    parser.add_argument('contig2', help ='path to contig 2')
    parser.add_argument('alignment_dir' , help='alignment directory')
    parser.add_argument('mugsy_path', help = 'path to mugsy')
    args = parser.parse_args()

    # Assumes that files are named strain.contigextension.renamed.mugsy
    strain1 = '.'.join(args.contig1.split('/')[-1].split('.')[0:-3])
    strain2 = '.'.join(args.contig2.split('/')[-1].split('.')[0:-3])
    out_id_1 = 'TEMPOUT1337' + ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))
    out_id_2 = 'TEMPOUT1337' + ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))

    prefix = out_id_1 + out_id_2
    correct_name = '{strain1}_@_{strain2}.maf'.format(strain1 = strain1, strain2 = strain2) 

    os.system('cp {contig1} {randomcontigname1}'.format(contig1 = args.contig1, randomcontigname1 = out_id_1))
    os.system('cp {contig2} {randomcontigname2}'.format(contig2 = args.contig2, randomcontigname2 = out_id_2))

    os.system('{mugsypath} --directory {align_directory} --prefix {prefix} {randomcontigname1} {randomcontigname2}'.format(mugsypath = args.mugsy_path, align_directory = args.alignment_dir, prefix = prefix, randomcontigname1 = out_id_1, randomcontigname2 = out_id_2))

    os.system('rm {random_contig1}'.format(random_contig1 = out_id_1))
    os.system('rm {random_contig2}'.format(random_contig2 = out_id_2))
    os.system('rm {prefix}.mugsy.log'.format(prefix = prefix))
    os.system('mv {random_alignment_name} {correct_name}'.format(random_alignment_name = args.alignment_dir + '/' + prefix + '.maf', correct_name = args.alignment_dir + '/' + correct_name))

if __name__=='__main__':
    main()
