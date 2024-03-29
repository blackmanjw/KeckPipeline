from utilities.utilities import import_user_config
import tools

## This is the main executable to run the functions in tools.py
## Run the program as "python calib.py" with the appropriate flags below depending on what you want to do.
## Eg. To move files, use the -m flag.

if __name__ == "__main__":
    # Allow user to in CLI optionally overwrite the config values
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
                ***********************************************************************
                TestyMcTesterson
                \tfunkytowntesting
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                '''))

    parser.add_argument('-c', '--config', type=str, nargs=1, help='usage -c /path/to/config.json')
    parser.add_argument('-l', '--list', type=str, nargs=1, help='usage -l file.list')
    parser.add_argument('-r', '--rename', action="store_true")
    parser.add_argument('-s', '--source', type=str, nargs=2, help='usage -c /path/to/config.json')
    parser.add_argument('-fs', '--flatsubtract', action="store_true")
    parser.add_argument('-dl', '--darklist', action="store_true")
    parser.add_argument('-fl', '--flatlist', action="store_true")
    parser.add_argument('-ds', '--darksubtract', action="store_true")
    parser.add_argument('-dc', '--darkcombine', action="store_true")
    parser.add_argument('-fc', '--flatcombine', action="store_true")
    parser.add_argument('-sky', '--skycombine', action="store_true")
    parser.add_argument('-t', '--test', action="store_true")

    args = parser.parse_args()
    #print(args)
    if args.config is None:
        config = import_user_config('config/my_config.json')
    else:
        config = import_user_config(args.config[0])
    if args.darkcombine:
        config.update({'dark_combine': True})
    if args.darksubtract:
        config.update({'dark_subtract': True})
    if args.darklist:
        tools.darklist('Darks')
    if args.source is None:
        dir='Data'
    else:
        source_dir=args.source[0]
        dest_dir = args.source[1]
    if args.rename:
        tools.rename(source_dir,dest_dir)
    #else:
    #    tools.darklist(args.darklist[0])
    if args.flatlist:
        tools.flatlist('Flats')
    if args.darkcombine:
        tools.darklist('Darks')
        tools.darkcombine()
    if args.flatcombine:
        config.update({'flat_combine': True})
        print('Please combine flats!!!')
    if args.flatsubtract:
        config.update({'flat_subtract': True})
    if args.flatsubtract:
        config.update({'move_data': True})
        tools.movedata('2018-05-05',3)
    if args.list:
        config.update({'list': args.list[0]})
    if args.test:
        print('We Win!!!!!')

    #movedata('2018-05-05',3)
    #for key, value in config.items():
    #	print(f"Key: {key}\t\tValue: {value}")