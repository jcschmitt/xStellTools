import os
import shutil
import numpy as np
import f90nml

gen = {}  # make a dictionary

gen['vmec_ext'] = 'aten_a3b25'
gen['vmec_input_in_p1'] =  'input/' + 'input.' + gen['vmec_ext']
gen['vmec_input_in_p2'] =  'input/' + 'input.' + gen['vmec_ext'] + '_p2'
gen['vmec_input_out'] = 'input.' + gen['vmec_ext']
gen['boozer_input_out'] = 'in_booz.' + gen['vmec_ext']

#gen['slurm_file_in'] = 'input/slurm_vmec.sh'
#gen['slurm_file_out'] = 'slurm_vmec.sh'
gen['slurm_file_in'] = 'input/slurm_vmec_2.sh'
gen['slurm_file_out'] = 'slurm_vmec_2.sh'

gen['mpol'] = (8, 10, 12, 14, 16, 18)
gen['ntor'] = (8, 12, 16, 18, 20, 24, 26)

gen['ns'] = ([16, 35, 51], [16, 35, 51], [16, 35, 51])#,
#             [16, 35, 51, 65], [16, 35, 51, 65], [16, 35, 51, 65],
#             [16, 35, 51, 65, 81], [16, 35, 51, 65, 81], [16, 35, 51, 65, 81],
#             [16, 35, 51, 65, 81, 101], [16, 35, 51, 65, 81, 101], [16, 35, 51, 65, 81, 101],
#             [16, 35, 51, 65, 81, 101, 128], [16, 35, 51, 65, 81, 101, 128], [16, 35, 51, 65, 81, 101, 128],
#             [16, 35, 51, 65, 81, 101, 128, 151], [16, 35, 51, 65, 81, 101, 128, 151], [16, 35, 51, 65, 81, 101, 128, 151])
gen['niter'] = ([20000, 20000, 20000], [20000, 20000, 20000], [20000, 20000, 20000])#,
#                [20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000],
#                [20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000],
#                [20000, 20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000, 20000],
#                [20000, 20000, 20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000, 20000, 20000],
#                [20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000], [20000, 20000, 20000, 20000, 20000, 20000, 20000, 20000])
gen['ftol'] =([1.0e-5, 5.0e-7, 5.0e-9], [1.0e-5, 2.0e-7, 2.0e-9], [1.0e-5, 1.0e-7, 1.0e-9])#,
#              [1.0e-5, 5.0e-7, 5.0e-9, 1.0e-9], [1.0e-5, 2.0e-7, 2.0e-9, 1.0e-9], [1.0e-5, 1.0e-7, 1.0e-9, 5.0e-10],
#              [1.0e-5, 5.0e-7, 5.0e-9, 1.0e-9], [1.0e-5, 2.0e-7, 2.0e-9, 1.0e-9], [1.0e-5, 1.0e-7, 1.0e-9, 5.0e-10],
#              [1.0e-5, 5.0e-7, 5.0e-9, 1.0e-9], [1.0e-5, 2.0e-7, 2.0e-9, 1.0e-9], [1.0e-5, 1.0e-7, 1.0e-9, 5.0e-10],
#              [1.0e-5, 5.0e-7, 5.0e-9, 1.0e-9], [1.0e-5, 2.0e-7, 2.0e-9, 1.0e-9], [1.0e-5, 1.0e-7, 1.0e-9, 5.0e-10],
#              [1.0e-5, 5.0e-7, 5.0e-9, 1.0e-9], [1.0e-5, 2.0e-7, 2.0e-9, 1.0e-9], [1.0e-5, 1.0e-7, 1.0e-9, 5.0e-10] )
 
gen['mboz'] = (20, 40, 56, 64, 80, 88, 104, 112, 132)
gen['nboz'] = (20, 40, 56, 64, 80, 88, 104, 112, 132)

# read in namelists
my_parser = f90nml.Parser()
my_parser.default_start_index = 0
nml_in = my_parser.read(gen['vmec_input_in_p1'])

# directory structure and name(s)
# aten_a3b25_mpol#_ntor#_ns#_ftol#/input.*, in-booz.*_mboz#_nboz#

# loop over target values, generate output directory names and filenames, write files, close files
for i_mpol, mpol in zip(range(0, len(gen['mpol'])), gen['mpol']):
  for i_ntor, ntor in zip(range(0, len(gen['ntor'])), gen['ntor']):
    for i_nsiterftol, ns, niter, ftol in zip(range(0, len(gen['ns'])), gen['ns'], gen['niter'], gen['ftol']):

      dir_out = ( 'output/' + gen['vmec_ext'] + '_mpol' + str(i_mpol) + '_ntor' + str(i_ntor) +
                  '_nsiterftol' + str(i_nsiterftol) ) 
      vmec_out = gen['vmec_input_out']

      # make directory
      print('<----Making ' + dir_out)
      os.makedirs(dir_out, exist_ok=True)

      this_new_file = f90nml.Namelist()

      # vmec
      # create nml data
      this_new_file['indata'] = {}
      for variable in nml_in['indata']:
        this_new_file['indata'][variable] = nml_in['indata'][variable]

      # find properties to adjust
      # adjust them
      this_new_file['indata']['mpol'] = mpol
      this_new_file['indata']['ntor'] = ntor
      this_new_file['indata']['ns_array'] = ns
      this_new_file['indata']['niter_array'] = niter
      this_new_file['indata']['ftol_array'] = ftol

      # make files
      # full (relative) path + filename
      file_out_name = os.path.join(dir_out, vmec_out)
      rbczbs_file_name = gen['vmec_input_in_p2']
      # write the data

      this_new_file.write(file_out_name, force=True)
      # end loop over target values

      # add on RBC, ZBS stuff (since f90nml read it in wrong)
      # read the file and output the text right back to the file, skipping the
      # 'trailing backslash'
      rbczbs_file = open(rbczbs_file_name, 'r')
      rbczbs_text = rbczbs_file.readlines()
      rbczbs_file.close()
      input_file = open(file_out_name, 'r')
      file_text = input_file.readlines()
      input_file.close()
      input_file = open(file_out_name, 'w')
      for this_line in file_text:
          if (this_line != '/\n'):
              # in most cases, just print the line
              input_file.write(this_line)
          else:
              # in the case of the trailing backslash, print a newline
              input_file.write('\n')

      # now do the append of the Initial Position info (axis and boundary)
      for this_line in rbczbs_text:
          input_file.write(this_line)

      input_file.write('/\n')


      for i_mboz, mboz in zip(range(0, len(gen['mboz'])), gen['mboz']):
        for i_nboz, nboz in zip(range(0, len(gen['nboz'])), gen['nboz']):
          booz_out_name = (gen['boozer_input_out'] + '_mboz' + str(i_mboz) +
                           '_nboz' + str(i_nboz))
          boozfile_out_name = os.path.join(dir_out, booz_out_name)
          theboozfile = open(boozfile_out_name, 'w')
          nextline = (str(mboz) + '   ' + str(nboz) + '\n')
          theboozfile.write(nextline)
          vmec_wout = gen['vmec_ext']
          nextline = (vmec_wout + '\n')
          theboozfile.write(nextline)
          nextline = (str(list(range(0,1+ns[-1]))).replace('[', '').replace(']', '').replace("'", ''))

          theboozfile.write(nextline)
          theboozfile.close()

      # copy the slurm file
      source_file = gen['slurm_file_in']
      destination_file = os.path.join(dir_out, gen['slurm_file_out'])
      shutil.copyfile(source_file, destination_file)

print("<---end generator.py")


