#!/usr/bin/python3
# -*- coding: utf-8 -*-


import os
import signal
import subprocess
from threading import Timer


# default values
data_dir = ".." + os.path.sep + "datasets" + os.path.sep
is_shell = True
timeout = 36000


def call(_args, _timeout):
	p = subprocess.Popen(_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=is_shell, start_new_session=True)
	timer = Timer(_timeout, lambda process: os.killpg(process.pid, signal.SIGUSR1), [p])
	try:
		timer.start()
		stdout, stderr = p.communicate()
		return_code = p.returncode
		return str(stdout, encoding="utf8"), str(stderr, encoding="utf8"), return_code
	finally:
		timer.cancel()


def run_cmd(_cmd: str):
	print("run " + _cmd + " start")
	print(call(_cmd, timeout))
	print("run " + _cmd + " end")


run = ""
run_cmd(run)
