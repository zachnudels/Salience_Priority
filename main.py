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

    sj_num = int(sys.argv[1])  # subject number
    exp_num = int(sys.argv[2])  # run number

    # TODO: Check if practice always in same file?
    # task name dictionary
    print(f"Running experiment {exp_num} for subject {sj_num}")


    # make output dir
    if params['paths']['curr_dir'] == 'lab':
        output_dir = op.join(params['paths']['data_pth']['lab'], 'output','sourcedata', 'sub-{sj}'.format(sj=sj_num))
    else:
        base_dir = op.split(os.getcwd())[0] # main path for all folders of project
        output_dir = op.join(base_dir,'output','sourcedata', 'sub-{sj}'.format(sj=sj_num))

    # if output path doesn't exist, create it
    if not op.isdir(output_dir): 
        os.makedirs(output_dir)
    print('saving files in %s'%output_dir)

    # string for output data
    output_str = 'sub-{sj}_ses-{ses}_task-{task}'.format(sj = sj_num, ses = ses_type, task = tasks[exp_type])

    # if file already exists
    behav_file = op.join(output_dir,'{behav}_events.tsv'.format(behav=output_str))
    if op.exists(behav_file): 
        print('file already exists!')

        overwrite = ''
        while overwrite not in ('y','yes','n','no'):
            overwrite = input('overwrite %s\n(y/yes/n/no)?: '%behav_file)

        if overwrite in ['no','n']:
            raise NameError('Run %s already in directory\nstopping experiment!'%behav_file)

    exp_sess = SingletonSession(output_str=f"behavioural_data_mieke{sj_num}",
                                output_dir=output_dir,
                                eyetracker_on=True,
                                settings_file='experiment_settings.yml',
                                subject_number=sj_num,
                                exp_num=exp_num)

    exp_sess.run()


if __name__ == '__main__':
    main()
