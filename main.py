# import relevant packages
import sys
import os
import os.path as op
# import appnope
import yaml

from singleton_session import SingletonSession


# define main function
def main():
    # load settings from yaml
    with open(op.join(os.getcwd(), 'experiment_settings.yml'), 'r') as f_in:
        params = yaml.safe_load(f_in)

    # take user input

    # define participant number and open json parameter file
    if len(sys.argv) < 2:
        raise NameError('Please add subject number (ex:1) '
                        'as 1st argument in the command line!')

    elif len(sys.argv) < 3:
        raise NameError('Please add session number (ex:1) '
                        'as 2nd argument in the command line!')

    sj_num = str(sys.argv[1]).zfill(3)  # subject number
    exp_num = str(sys.argv[2])  # run number

 # TODO: Check if practice always in same file?
    # task name dictionary
    exps = params['study']['no_exps']

    exp_num = ''
    while exp_num not in range(exps):
        exp_num = input(f"Which experiment to run ({exps.join('/ ')})?: ")

    print(f"Running experiment {exp_num} for subject {sj_num}")

    # make output dir
    if params['paths']['curr_dir'] == 'lab':
        base_dir = params['paths']['data_pth']['lab']
    else:
        base_dir = op.split(os.getcwd())[0]  # main path for all folders of project
    output_dir = op.join(base_dir, 'output', 'sourcedata', f"{sj_num}")

    # if output path doesn't exist, create it
    if not op.isdir(output_dir):
        os.makedirs(output_dir)
    print(f"Saving files in {output_dir}")

    # TODO: string for output data
    output_str = f"behavioural_data_mieke{sj_num}"

    # if file already exists
    # behav_file = op.join(output_dir, f"{output_str}_events.tsv")
    # if op.exists(behav_file):
    #     print('File already exists!')
    #
    #     overwrite = ''
    #     while overwrite not in ('y', 'yes', 'n', 'no'):
    #         overwrite = input(f"overwrite {behav_file} \n(y/yes/n/no)?:")
    #
    #     if overwrite in ['no', 'n']:
    #         raise NameError(f"Run {behav_file} already in directory\nstopping experiment!")

    # load appropriate class object to be run

    exp_sess = SingletonSession(output_str=output_str,
                                output_dir=output_dir,
                                eyetracker_on=True,
                                settings_file=params,
                                subject_number=sj_num,
                                exp_num=exp_num)

    exp_sess.run()


if __name__ == '__main__':
    main()
