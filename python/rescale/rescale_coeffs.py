data = open('input.QHS_Rstart_1_513_32polmodes_18x24_axis_v2').readlines()
wf = open('input.hsx_scaled_to_QHS46','w')

# scaling HSX to have the same average major radius as QHS46
scale = 6.025006 / 1.211364

for line in data:

#Do the first two lines, the axis lines
  if 'RAXIS_CC' in line or 'ZAXIS_CC' in line:
    d = line.split()
    wf.write('  ')
    wf.write(d[0]+' '+d[1]+' ')
    for k in xrange(2,len(d)):
        v = float(d[k])*scale
        wf.write(str(v))
        wf.write(' ')
    wf.write('\n')

  elif 'RBC' in line and 'ZBS' in line:
    wf.write('  ')
    d = line.split()
    skip = False
    for index, val in enumerate(d):
      if skip:
        skip = False
        continue
      wf.write(val)
      wf.write(' ')
      if val == '=':
        skip = True
        wf.write("{:.12e}".format((float(d[index+1])*scale)))
        wf.write('   ')
    wf.write('\n')
    
  else:
    wf.write(line)
            
