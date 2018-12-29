import sys

def do_error(msg):	# global error handler
    sys.stderr.write('ERROR: ' + msg + '\n')
    sys.exit(-1)
