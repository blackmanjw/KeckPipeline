

from utilities.utilities import import_user_config
from the_real_code.the_real_code import yeah_baby 

if __name__ == "__main__":
    # Allow user to in CLI optionally overwrite the config values
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
                ***********************************************************************
                Back-testing suite for algorithmic trading (development ongoing)
                The user can edit the included turtle_config.json, pass custom config
                following the same format), and/or optionally bypass a selection of
                config parameters:
                \tLong entry/exit
                \tShort entry/exit
                \tATR multiplier
                default usage:
                \tpython3 back_testing.py
                bypass example:
                \tpython3 back_testing.py -c my_turtles.json -l 22 10 -s 15 5 -a 2.5
                in addition to the above parameters, you can have the code select a
                random period between 2016-01-01 and the present date. To do this chose
                a desire number of months (say 10) and run:
                \tpython3 back_testing.py --random 10
                including any additional flags  from above that you want
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                '''))
    parser.add_argument('-c', '--config', type=str, nargs=1, help='usage -c /path/to/config.json')
    parser.add_argument('-l', '--longs', type=int, nargs=2, help='usage -l <long_entry> <long_exit>')
    parser.add_argument('-s', '--shorts', type=int, nargs=2, help='usage -s <short_entry> <short_exit>')
    parser.add_argument('-a', '--atrmult', type=float, nargs=1, help='usage -a <atr multiplier>')
    parser.add_argument('-p', '--pyramid', type=int, nargs=1, help='usage -p <pyramid_limit>')
    parser.add_argument('--random', type=int, nargs=1, help='usage --random <months_to_test>')

    args = parser.parse_args()
    if args.config is None:
        config = import_user_config('config_files/my_config.json')
    else:
        config = import_user_config(args.config[0])
    if args.longs is not None:
        config.update({'enter_long': args.longs[0]})
        config.update({'exit_long': args.longs[1]})
    if args.shorts is not None:
        config.update({'enter_short': args.shorts[0]})
        config.update({'exit_short': args.shorts[1]})
    if args.atrmult is not None:
        config.update({'atr_multiplier': args.atrmult[0]})
    if args.pyramid is not None:
        config.update({'pyramid_limit': args.pyramid[0]})
    if args.random is not None:
        config.update({"random_period": True})
        config.update({"period_length": 10})

    yeah_baby('it works')
    for key, value in config.items():
    	print(f"Key: {key}\t\tValue: {value}")