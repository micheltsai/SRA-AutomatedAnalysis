import subprocess


def run_cmd(cmd, stdout=False, stderr=False):
    if not stdout:
        stdout_arg = subprocess.DEVNULL
    else:
        stdout_arg = subprocess.PIPE
    if not stderr:
        stderr_arg = subprocess.DEVNULL
    else:
        stderr_arg = subprocess.PIPE
    p = subprocess.run(cmd, stdout=stdout_arg, stderr=stderr_arg, check=True, shell=True)
    return p
