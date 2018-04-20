#!/usr/bin/python
import re
import numpy as np
import json
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.sparse import coo_matrix
from math import sqrt, exp, pi, cos
#now use gen_input
import sys
sys.path.append('scripts')

#import armadillo as arma
import syst
import bath
from unitconverter import cm2au, fs2au, kt2au, au2ar, am2au
#Python3 use range instead of xrange
try:
    xrange
except NameError:
    xrange = range

def encode (i_list,cum_prod,nlen):
    i_rank = 0
    for i in xrange(nlen-1):
        i_rank += i_list[i]*cum_prod[i]
    i_rank += i_list[-1]
    return i_rank

def decode (i_rank,cum_prod,nlen):
    i_list = []
    for i in xrange(nlen-1):
        i_list.append(i_rank/cum_prod[nlen-i-2])
        i_rank %= cum_prod[nlen-i-2]
    i_list.append(i_rank)
    return i_list

def savelog(dic):
    for i in dic:
        with open(i,'wb') as f:
            np.savetxt(f,dic[i])

def main ():

    with open('morse.json','r') as f:
        d = json.load(f)
    ne = d['ne']
    ie = d['ie']
    modes = d['modes']
    gxy = d['gxy']
    for m in modes:
        m['m0'] *= am2au
        m['ri'] /= au2ar
        m['rf'] /= au2ar
        m['re'] /= au2ar
    vmn = [[v*cm2au for v in colv] for colv in gxy]

    n_list = np.array([m['nr'] for m in modes])
    ndim = np.prod(n_list)
    nmod = len(modes)
    cum_prod = np.zeros(nmod,dtype=int)
    cum_prod[0] = modes[nmod-1]['nr']
    for i in xrange(1,nmod):
        cum_prod[i] = cum_prod[i-1]*modes[nmod-i-1]['nr']

    with open('wavetrun.dat','r') as f:
        wavefunc=np.loadtxt(f)
    with open('energr.dat','r') as f:
        energies=np.loadtxt(f)

    qmod = []
    hami = []
    dipo = []
    dip2 = []
    for imod in xrange(nmod):
        m = modes[imod]
        q = np.empty(ndim,dtype=float)
        d = np.empty(ndim,dtype=float)
        d2= np.empty(ndim,dtype=float)
        for i in xrange(ndim):
            i_list = decode(i,cum_prod,nmod)
            d[i] = m['ri']+i_list[imod]*(m['rf']-m['ri'])/m['nr']
            d2[i]= (d[i]-m['re'])*(d[i]-m['re'])
            q[i] = m['q2']*(d[i]-m['re'])*(d[i]-m['re'])/2+m['q1']*(d[i]-m['re']) 
        qwave = np.array([[q[i]*wavefunc[i,j] for j in xrange(ne)] for i in xrange(ndim)])
        dwave = np.array([[d[i]*wavefunc[i,j] for j in xrange(ne)] for i in xrange(ndim)])
        d2wave= np.array([[d2[i]*wavefunc[i,j] for j in xrange(ne)] for i in xrange(ndim)])
        qmod.append(np.dot(wavefunc.T,qwave) )
        dipo.append(np.dot(wavefunc.T,dwave) )
        dip2.append(np.dot(wavefunc.T,d2wave))
#test new geninput
        ini={
    	"syst": {
    	    "hamsFile": "inp_hams.mat",
    	    "qmdsFile": "inp_qmds.mat",
    	    "dipsFile": "inp_dips.mat"
    	},
    	"bath": {
            "npsd": 6,
            "pade": 2,
            "temp": 300.0,
	    "modeFile": "inp_mode.mat",
    	    "etalFile": "inp_etal.mat",
    	    "etarFile": "inp_etar.mat",
    	    "etaaFile": "inp_etaa.mat",
    	    "expnFile": "inp_expn.mat",
    	    "bdipFile": "inp_bdip.mat",
    	    "delrFile": "inp_delr.mat"
        },
    	"hidx": {
    		"lmax": 5,
    		"nmax": 5000,
    		"ferr": 2.0e-7,
	    "expnFile": "inp_expn.mat"
    	}
    }

    qmds = np.array(qmod)*(1.0+0.0J)
    hams = np.diag(np.array(energies))*(1.0+0.0J)
    syst.init(ini['syst'],hams,qmds)
    mutot = np.empty(ne,dtype=float)
    pitot = np.empty(ne,dtype=float)
    qutot = np.empty(ne,dtype=float)
    for imod in xrange(nmod):
	    mutot=mutot+modes[imod]['d1']*dipo[imod]+modes[imod]['d2']*dip2[imod]

    # bath
    print ("temperature is K\n",ini['bath']['temp'])
    ini['bath']['temp'] *= kt2au
    ini['bath']['jomg'] = []
    for m in modes:
        lamd = m['jdru'][0]['lamd']*cm2au
        gamd = m['jdru'][0]['gamd']*cm2au
        drupara ={
                'jdru': [(lamd,gamd)],
                 }
        ini['bath']['jomg'].append(drupara)
    ini['bath']['amod'] = ['dru' for i in xrange(nmod)]
    ini['bath']['fmod'] = [0.0 for i in xrange(nmod)]
    ini['bath']['nmod'] = qmds.shape[0] 
    bath.init (ini['bath'])
    jsonInit = {"deom":ini,
    "rhot": {
        "dt": 0.005, 
        "nk": 10, 
        "nt": 2000, 
        "inistate": 0
    },
           }
    with open('input.json','w') as f:
        json.dump(jsonInit,f,indent=4) 
        
if __name__ == '__main__':
    main()
