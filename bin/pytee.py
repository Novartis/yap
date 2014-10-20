"""An implementation of Unix's 'tee' command in Python

Implementing 'tee' in Python allows the ease of writing output once, which will
be written to many file/file-like objects.

Benefits:
* Cross-platform (e.g. a pipe to the actual tee command won't work on Windows)
* The ability to write to more than just files. Any file-like object that
implements write() and flush() is acceptable.

Sample Usage:

import pytee

tee = pytee.create_tee([ '/tmp/tee-test-1', '/tmp/tee-test-2' ], mode='a')
print >>tee, "a" * 100000,
tee.close()

# Need to wait for the child process to finish writing to the file(s)
# before we can measure the amount of data written to the file(s).
os.wait()

for filename in files:
with open(filename, 'r') as fh:
chars = len(fh.read())
print "File '%s' has %d chars" % (filename, chars)
"""


import sys
import os
import subprocess
from netsa.util.shell import *
import os
from string import Template

__author__ = 'Brandon Sandrowicz <brandon@sandrowicz.org>'
__version__ = '0.1'

valid_modes = ['a', 'w']


def create_tee(files, mode, buffer_size=128):
    """Get a file object that will mirror writes across multiple files objs

Options:
files A list of files and/or file objects. All strings will be
treated as file paths and opened for writing. Everything
else is assumed to be a file-like object that implements
both the write() and flush() methods.

mode Which mode to use when opening new files. Valid values
are 'a' (append) and 'w' (overwrite).

buffer_size
Control the size of the buffer between writes to the
resulting file object and the list of files.
"""
    if mode not in valid_modes:
        raise IOError("Only valid modes to create_tee() are: %s" %
                      ', '.join(valid_modes))

    tee_list = []
    for file in files:
        if isinstance(file, str):
            fp = open(file, mode)
            tee_list.append(fp)
        else:
            tee_list.append(file)

    pipe_read, pipe_write = os.pipe()
    pid = os.fork()
    if pid == 0:
        # Child -- Read bytes from the pipe and write them to the specified
        # files.
        try:
            # Close parent's end of the pipe
            os.close(pipe_write)

            bytes = os.read(pipe_read, buffer_size)
            while(bytes):
                for file in tee_list:
                    file.write(bytes)
                    file.flush()
                    # TODO maybe add in fsync() here if the fileno() method
                    # exists on file

                bytes = os.read(pipe_read, buffer_size)
        except:
            pass
        finally:
            os._exit(255)
    else:
        # Parent -- Return a file object wrapper around the pipe to the
        # child.
        return os.fdopen(pipe_write, 'w')
