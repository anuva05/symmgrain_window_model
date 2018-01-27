# parses the simulation result file to do convergence analysis
import sys
datafilename = sys.argv[1]
strainerr = 'STRAIN FIELD ERROR'

stresserr = 'STRESS FIELD ERROR'
#datafile = file('slurm-2233782.out')
datafile = file(datafilename)
found = False
ctr1=0
ctr2=0
strain_f = []
stress_f = []
for line in datafile:
        if strainerr in line:
            found = True
            l = line.split('=')
            strain_f.append( l[1])
            ctr1 = ctr1 + 1
        if stresserr in line:
            found = True
            l = line.split('=')
            stress_f.append( l[1] )
            ctr2 = ctr2+1 

print "\n".join(strain_f)
print "stress:"
#print stress_f
print "\n".join(stress_f)
datafile.close()

