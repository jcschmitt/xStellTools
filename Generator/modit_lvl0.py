#!/usr/bin/env python3

import sys
import shutil

# Parse input parameters
filename = sys.argv[1]

# Copy original file to 'backup'
backup_filename = filename + '.bak'

shutil.copyfile(filename, backup_filename)

# Open original file
orig_file = open(backup_filename, 'r')

# Open destination file
new_file = open(filename, 'w')

# read lines from orig file and close it
orig_text = orig_file.readlines()
orig_file.close()

# Loop 
for this_line in orig_text:
   if ('NS_ARRAY' in this_line):
      new_line = '  NS_ARRAY =  32 51 81 101 128\n'
   elif ('FTOL_ARRAY' in this_line):
      new_line = '  FTOL_ARRAY = 5e-7 5e-8 5e-9 1e-12 1e-15\n'
   elif ('NITER_ARRAY' in this_line):
      new_line = '  NITER = 99999\n'
   elif ('MPOL' in this_line):
      new_line = '  MPOL = 10\n'
   elif ('NTOR' in this_line):
      new_line = '  NTOR = 12\n  NZETA = 50\n'
   else:
      new_line = this_line
   # writeline to dest file
   new_file.write(new_line)

# Close file
new_file.close()

