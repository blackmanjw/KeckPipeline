def import_user_config(filename):
    import os
    import json
    cwd = os.getcwd()
    with open(os.path.join(cwd, filename), 'r') as config_file:
        user_config = json.load(config_file)
    return user_config
