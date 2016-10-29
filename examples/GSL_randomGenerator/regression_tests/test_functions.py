import subprocess as sp


def run_and_check_sim(command, file):
    command_redirect = command + " > " + file
    proc = sp.Popen(command_redirect, shell=True)
    proc.communicate()
    output = open(file).read()
    if "ERROR" in output:
        print output
        return False
    else:
        return True


def run(command1):
    proc = sp.Popen(command1, shell=True)
    proc.communicate()


def run_return_stdout(command1):
    proc = sp.Popen(command1, shell=True, stdout=sp.PIPE)
    return proc.communicate()[0]


def identical_stdout(command1, command2):
    value1 = run_return_stdout(command1) 
    value2 = run_return_stdout(command2) 
    if value1 == value2:
        return True
    else:
        #print command1, value1
        #print command2, value2
        return False


def equal_output(file1, file2, option):
        c1 = 'tft_extract.py ' + file1 + " " + option
        c2 = 'tft_extract.py ' + file2 + " " + option
        return identical_stdout(c1, c2)

