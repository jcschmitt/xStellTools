import os
import numpy as np
import f90nml

gen = {}  # make a dictionary

gen['vmec_ext'] = 'aten_a3b25'
gen['vmec_input_in'] =  'input/' + 'input.' + gen['vmec_ext']
gen['vmec_input_out'] = 'input.' + gen['vmec_ext']
gen['boozer_input_out'] = 'in_booz.' + gen['vmec_ext']

gen['slurm_file_in'] = 'input/slurm_vmec.sh'
gen['slurm_file_out'] = 'slurm_vmec.sh'

gen['mpol'] = (6, 8, 10, 12, 14, 16)
gen['ntor'] = (8, 10, 12, 14, 16, 18, 20, 24)

gen['ns'] = ('16, 35, 51', '16, 35, 51')
gen['niter'] = ('20000, 20000, 20000', '20000, 20000, 20000')
gen['ftol'] =('1.0e-5, 5.0e-7, 5.0e-9', '1.0e-5, 2.0e-7, 2.0e-9' )
 
gen['mboz'] = (8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120)
gen['nboz'] = (8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120)

# read in namelists
my_parser = f90nml.Parser()
my_parser.default_start_index = 0
nml_in = my_parser.read(gen['vmec_input_in'])

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
      # write the data

      this_new_file.write(file_out_name, force=True)
      # end loop over target values

      for i_mboz, mboz in zip(range(0, len(gen['mboz'])), gen['mboz']):
        for i_nboz, nboz in zip(range(0, len(gen['nboz'])), gen['nboz']):
          booz_out_name = (gen['boozer_input_out'] + '_mboz' + str(i_mboz) +
                           '_nboz' + str(i_nboz))
          boozfile_out_name = os.path.join(dir_out, booz_out_name)
          theboozfile = open(boozfile_out_name, 'w')
          nextline = (str(mboz) + '   ' + str(nboz) + '\n')
          theboozfile.write(nextline)
          nextline = (vmec_out + '\n')
          theboozfile.write(nextline)
          theboozfile.close()

      # copy the slurm file
      source_file = gen['slur_file_in']
      destination_file = os.path.join(dir_out, gen['slur_file_out'])
      shutil.copyfile(soure_file, destination_file)

print("<---end generator.py")


