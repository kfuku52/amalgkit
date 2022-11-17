import os
import sys


def check_directory(args):
    work_dir = os.path.join(args.out_dir, 'config')
    path_config = os.path.join(work_dir, args.config_dir)

    if not os.path.exists(os.path.join(work_dir)):
        os.makedirs(os.path.join(work_dir))

    print('Checking config directory ...')
    if os.path.exists(path_config):
        print(path_config, ' already exists.')
        if args.overwrite:
            print('--overwrite is set! Any default config files in ', path_config, ' will be overwritten.')
            print('Continuing.')
        else:
            print('Exiting.')
            sys.exit()
    else:
        os.makedirs(path_config)

    return


def write_files(args, ext = '.config'):
    print('Creating files...')
    work_dir = os.path.join(args.out_dir, 'config')
    # control_term
    path_file = os.path.join(work_dir, args.config_dir, 'control_term' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (attribute)[TAB](control term)\n'
            f.writelines([l1,l2,l3])

    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # exclude_id
    path_file = os.path.join(work_dir, args.config_dir, 'exclude_id' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (arbitrary reason of exclusion)[TAB](Bioproject/Biosample/Experiment/Run/SRA_Primary/SRA_Study ID)\n'
            f.writelines([l1, l2, l3])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # exclude_keyword
    path_file = os.path.join(work_dir, args.config_dir, 'exclude_keyword' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (comma-delimited attributes to be searched)[TAB](arbitrary reason of exclusion)[TAB](bad keyword)\n'
            f.writelines([l1, l2, l3])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # give_value
    path_file = os.path.join(work_dir, args.config_dir, 'give_value' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (ID)[TAB](column)[TAB](value)\n'
            l4 = '# df.loc[df[one_of_ID_columns]==ID, column] = value\n'
            f.writelines([l1, l2, l3, l4])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # group_attribute
    path_file = os.path.join(work_dir, args.config_dir, 'group_attribute' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (aggregated to)[TAB](aggregated from)\n'
            f.writelines([l1, l2, l3])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # group_tissue
    path_file = os.path.join(work_dir, args.config_dir, 'group_tissue' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (aggregated to)[TAB](aggregated from)\n'
            f.writelines([l1, l2, l3])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # orthographical_variant
    path_file = os.path.join(work_dir, args.config_dir, 'orthographical_variant' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (comma-delimited attributes to be searched)[TAB](replaced to)[TAB](replaced from)\n'
            f.writelines([l1, l2, l3])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # replace_value
    path_file = os.path.join(work_dir, args.config_dir, 'replace_value' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# (column1)[TAB](column1_value)[TAB](column2)[TAB](column3)\n'
            l4 = '# df.loc[df[column1]==column1_value, column2] = df.loc[df[column1]==column1_value, column3]\n'
            f.writelines([l1, l2, l3, l4])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # rescue_id
    path_file = os.path.join(work_dir, args.config_dir, 'rescue_id' + ext)
    try:
        with open(path_file, 'w') as f:
            with open(path_file, 'w') as f:
                l1 = '# case-insensitive\n'
                l2 = '# regular expressions allowed\n'
                f.writelines([l1, l2])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # search_term_exclusion
    path_file = os.path.join(work_dir, args.config_dir, 'search_term_exclusion' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# Compatible with Entrez"s search field description\n'
            l4 = '# https://www.ncbi.nlm.nih.gov/books/NBK49540/\n'
            l5 = '# (Search Field)[TAB](Value)'
            f.writelines([l1, l2, l3, l4, l5])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # search_term_other
    path_file = os.path.join(work_dir, args.config_dir, 'search_term_other' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            l3 = '# Compatible with Entrez"s search field description\n'
            l4 = '# https://www.ncbi.nlm.nih.gov/books/NBK49540/\n'
            l5 = '# (Search Field)[TAB](Value)\n'
            f.writelines([l1, l2, l3, l4, l5])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # search_term_species
    path_file = os.path.join(work_dir, args.config_dir, 'search_term_species' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            f.writelines([l1, l2])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')

    # search_term_tissue
    path_file = os.path.join(work_dir, args.config_dir, 'search_term_tissue' + ext)
    try:
        with open(path_file, 'w') as f:
            l1 = '# case-insensitive\n'
            l2 = '# regular expressions allowed\n'
            f.writelines([l1, l2])
    except FileNotFoundError:
        print('The directory ', os.path.join(work_dir, args.config_dir), 'does not exist.')
    return


def check_files(args, config_file_list, ext = '.config'):
    work_dir = os.path.join(args.out_dir, 'config')
    print('Checking if all config files are present in ', os.path.join(work_dir, args.config_dir), '...')
    missing_files = []
    for config_file in config_file_list:
        path_config_file = os.path.join(work_dir, args.config_dir, config_file + ext)
        if os.path.exists(path_config_file):
            print('found ', config_file + ext, '...')
        else:
            print('could not find ', config_file + ext, 'in ', os.path.join(work_dir, args.config_dir))
            print('make sure amalgkit config ran correctly')
            missing_files.append(config_file + ext)
    if len(missing_files) > 0:
        print('some files are missing: ')
        for missing_file in missing_files:
            print(missing_file)
    else:
        print('All config files are present!')
    return


def config_main(args):
    ext = '.config'
    config_file_list = ['control_term', 'exclude_id', 'exclude_keyword', 'give_value', 'group_attribute',
                        'group_tissue',
                        'orthographical_variant', 'replace_value', 'rescue_id', 'search_term_exclusion',
                        'search_term_other',
                        'search_term_species', 'search_term_tissue']
    check_directory(args)
    write_files(args)
    check_files(args, config_file_list, ext)
