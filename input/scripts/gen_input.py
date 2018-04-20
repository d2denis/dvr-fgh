#!/home/panzj/anaconda2/bin/python
from math import sqrt, exp, pi, cos
import json
import numpy as np
import sys
sys.path.append('/data2/panzj/trie-deom/scripts')
import armadillo as arma
import syst
import bath
#from unitconverter import cm2au, fs2au, kt2au#, au2ar, am2au
#Exciton parameters
cm2au=fs2au=kt2au=1.0
E0=21.70#12400
E1=24.31#12540
V=-1.98#106
U=0
lamd=0.65
gama=0.99
lamb1=0.65*cm2au
gama1=6.22*cm2au
lamb2=0.65*cm2au
gama2=6.22*cm2au
temp=1.0*kt2au

def argv(i):
    try:
	    sys.argv[i]
    except :
	return 0;
    else:
	return 1;

if __name__ == '__main__':

    with open('default.json') as f:
        ini = json.load(f)

    
    # syst
    hams = np.zeros((4,4),dtype=complex)
    hams[2,1] = hams[1,2] = V
    hams[1,1] = E0
    hams[2,2] = E1
    hams[3,3] = E0+E1+U
    qmds = np.zeros((2,4,4),dtype=complex)
    qmds[0,1,1] =qmds[0,3,3] =  1.0
    qmds[1,2,2] =qmds[1,3,3] =  1.0
    syst.init (ini['syst'],hams*cm2au,qmds)

    nmod=qmds.shape[0]
   
    # bath
    ini['bath']['temp'] = temp
#    ini['bath']	['jomg']=[{"jdru":[(lamb1,gama1)]},{"jdru":[(lamb2,gama2)]}]
    ini['bath']	['jomg']=[{"jdru":[(lamd,gama)]} for i in xrange(nmod)]
    ini['bath']['nmod'] = nmod
    ini['bath']['pade'] = 1
    ini['bath']['npsd'] = 2
    ini['bath']['amod'] = ['dru' for i in xrange(nmod)]
    ini['bath']['fmod'] = [0.0/lamd for i in xrange(nmod)]
    bath.init (ini['bath'])
    
    # hidx
    ini['hidx']['lmax'] = 5
    ini['hidx']['nmax'] = 30000
    ini['hidx']['ferr'] = 2.0e-12

    jsonInit = {"deom":ini,
           }
    i=1
    while(argv(i) > 0): 
            string = sys.argv[i]
            print "you want to calulate " ,string 
    	    filename = string+".json"
            with open(filename) as f:
                dic = json.load(f)
    	    jsonInit.update(dic)
            i+=1
    
    with open('input.json','w') as f:
        json.dump(jsonInit,f,indent=4) 
