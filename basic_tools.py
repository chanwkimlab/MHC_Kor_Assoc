import subprocess

from path_configure import *

def run_subprocess(command,quiet=False,dry=False):
    print("------{}-----".format("RUN"))
    print(command)
    if dry:
        return
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if quiet==False:
        print("------{}-----".format("ERROR"))
        print(stderr.decode())
        print("------{}-----".format("OUTPUT"))
        print(stdout.decode())
    return stdout,stderr


def replace_characters(str_input,character_dict,catch_exception=True):
    for to_replace,value in character_dict.items():
        if catch_exception:
            try:
                str_input=str_input.replace(to_replace,value)
            except:
                pass
        else:
            str_input=str_input.replace(to_replace,value)
    return str_input  