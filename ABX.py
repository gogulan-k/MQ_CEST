#!/usr/bin/python3

# This code is designed for fitting multi-quantum chemical exchange saturation
# transfer data (MQ-CEST) see Karunanithy et al., J Phys. Chem. Lett., 2020
# (https://doi.org/10.1021/acs.jpclett.0c01322).

# ABX Functions
# Functions for 3 spin ABX systems
# GK UCL
# 29/10/19


import numpy as np
from scipy import linalg, kron, eye, transpose, sparse, stats
import sys, time
from math import log, ceil
import time, itertools, sys

e_mat = np.array([[1.0,0.0],[0.0,1.0]])
x_mat = np.array([[0.0,0.5],[0.5,0.0]])
y_mat = np.array([[0.0,-0.5j],[0.5j,0.0]])
z_mat = np.array([[0.5,0.0],[0.0,-0.5]])

ax = ['E']
a1 = ['Ix','Iy','Iz']
a2 = ['Rx','Ry','Rz']
a3 = ['Sx','Sy','Sz']

mats = [x_mat,y_mat,z_mat]

def buildProdOperators_IRS(ax,a1,a2,a3,mats):
    # build up 64 product operators sequentially for 3 spin system
    # use itertools product to get all combinations and as a helpful counter!
    # this returns a dictionary of 64x64 Product operator matrices and the associated
    # list of product operators

    prod = {}
    prodList = []

    # single spin first
    prodList.append('E')
    prod['E'] = np.kron(e_mat,kron(e_mat,e_mat))
    for i,vals in enumerate(itertools.product(ax,a1)): # for I spin
        prodList.append(vals[1])
        prod[vals[1]] = kron(mats[i],kron(e_mat,e_mat))

    for i,vals in enumerate(itertools.product(ax,a2)): # for R spin
        prodList.append(vals[1])
        prod[vals[1]] = kron(e_mat,kron(mats[i],e_mat))

    for i,vals in enumerate(itertools.product(ax,a3)): # for S spin
        prodList.append(vals[1])
        prod[vals[1]] = kron(e_mat,kron(e_mat,mats[i]))

    # double spin
    jj= [j for j in itertools.product([0,1,2],[0,1,2])]
    for i,vals in enumerate(itertools.product(a1,a2)): # for IR spins
        prodList.append(vals[0]+vals[1])
        prod[vals[0]+vals[1]] = 2.0*kron(mats[jj[i][0]],kron(mats[jj[i][1]],e_mat))

    for i,vals in enumerate(itertools.product(a1,a3)): # for IS spins
        prodList.append(vals[0]+vals[1])
        prod[vals[0]+vals[1]] = 2.0*kron(mats[jj[i][0]],kron(e_mat,mats[jj[i][1]]))

    for i,vals in enumerate(itertools.product(a2,a3)): # for RS spins
        prodList.append(vals[0]+vals[1])
        prod[vals[0]+vals[1]] = 2.0*kron(e_mat,kron(mats[jj[i][0]],mats[jj[i][1]]))

    # triple spin
    jjj= [j for j in itertools.product([0,1,2],[0,1,2],[0,1,2])]
    for i,vals in enumerate(itertools.product(a1,a2,a3)): # for IR spins
        prodList.append(vals[0]+vals[1]+vals[2])
        prod[vals[0]+vals[1]+vals[2]] = 4.0*kron(mats[jjj[i][0]],kron(mats[jjj[i][1]],mats[jjj[i][2]]))

    return prod, prodList

def rel_matrix(prodList):
     # create a dictionary of relaxation matrices for each
     # product operator
     # considers auto-relaxation only (no cross relaxation)
     # relies on product operators being in the right order
    rel = {}

    for i,prod in enumerate(prodList):
        rel[prod] = sparse.lil_matrix((127,127))
        rel[prod][i,i]+=1.0
        rel[prod][i+63,i+63]+=1.0
    return rel

def pop_rel(rel,rate_dic):
    # input dicitonary of sparse matrices (rel) and rate dictionary (rate_dic)
    # outputs a relaxation rate popuated 127x127 matrix
    rate_mat = np.zeros((127,127), dtype=complex)

    for prod in rel:
        if prod in rate_dic: # we only put in rates that are populated in the rate dicitonary
            rate_mat += rate_dic[prod]*rel[prod]
    return rate_mat


def CS_matrix(prod,prodList): # returns a matrix of 1s can be multiplied by I, R or S offsets
    # written for 127x127 matrix
    # we will write this for the current case (symmetric exchange between equally populated sites)
    # but could be easily generalised
    H_small_orig = 1.0*prod['Iz'] + 2.0*prod['Rz'] + 3.0*prod['Sz']
    H_small_ex = 1.0*prod['Iz'] + 3.0*prod['Rz'] + 2.0*prod['Sz']

    L_orig = kron(H_small_orig,prod['E']) - kron(prod['E'], transpose(H_small_orig))
    L_ex = kron(H_small_ex,prod['E']) - kron(prod['E'], transpose(H_small_ex))

    I_cs = np.zeros((127,127), dtype=complex)
    R_cs = np.zeros((127,127), dtype=complex)
    S_cs = np.zeros((127,127), dtype=complex)

    I_cs_ex = np.zeros((127,127), dtype=complex)
    R_cs_ex = np.zeros((127,127), dtype=complex)
    S_cs_ex = np.zeros((127,127), dtype=complex)
    # if we were to extended to unsymmetrical exchange case then we would need two more matrices here
    # I_cs_ex and S_cs_ex not needed in this case
    for i in range(1,64):
        for j in range(1,64):
            val_in = np.transpose(np.conjugate(prod[prodList[i]]).flatten())
            val_orig = np.real(1.0j*np.dot(val_in,np.dot(L_orig,(prod[prodList[j]].flatten())))) # this deals with evolution
            val_ex = np.real(1.0j*np.dot(val_in,np.dot(L_ex,(prod[prodList[j]].flatten()))))

            # NOTE: THERE'S A FUNNY NORMALISATION THING GOING ON HERE HENCE
            # VALUES ARE MULTIPLIED BY 2. EVERYTHING SEEMS TO BE BASICALLY
            # WORKING THOUGH

            if int(val_orig) == 2:
                I_cs[i,j]+=1.0
            elif int(val_orig) == -2:
                I_cs[i,j]-=1.0

            elif int(val_orig) == 4:
                R_cs[i,j]+=1.0
            elif int(val_orig) == -4:
                R_cs[i,j]-=1.0

            elif int(val_orig) == 6:
                S_cs[i,j]+=1.0
            elif int(val_orig) == -6:
                S_cs[i,j]-=1.0


            if int(val_ex) == 2:
                I_cs_ex[i+63,j+63]+=1.0
            elif int(val_ex) == -2:
                I_cs_ex[i+63,j+63]-=1.0

            elif int(val_ex) == 4:
                S_cs_ex[i+63,j+63]+=1.0
            elif int(val_ex) == -4:
                S_cs_ex[i+63,j+63]-=1.0

            elif int(val_ex) == 6:
                R_cs_ex[i+63,j+63]+=1.0
            elif int(val_ex) == -6:
                R_cs_ex[i+63,j+63]-=1.0

    return I_cs, R_cs, S_cs, I_cs_ex, R_cs_ex,S_cs_ex


def Pulse_matrix(prod_op): # work out where to add in pulses
    # prod op = 'I','R' or 'S' depending on which spin we want the pulse
    # to act on
    w1 = 1.0
    ham_rf_x = w1*(prod[prod_op+'x']) # assume pulse hits I and S here
    ham_rf_y =  w1*(prod[prod_op+'y']) # assume pulse hits I and S here

    x_rf = np.zeros((127,127), dtype = complex)
    y_rf = np.zeros((127,127), dtype = complex)

    L_rf_x = kron(ham_rf_x,prod['E']) - kron(prod['E'], transpose(ham_rf_x))
    L_rf_y = kron(ham_rf_y,prod['E']) - kron(prod['E'], transpose(ham_rf_y))
    for i in range(1,64):
        for j in range(1,64):
            inVal = np.transpose(np.conjugate(prod[prodList[i]]).flatten())
            val_x = np.real(1.0j*np.dot(inVal,np.dot(L_rf_x,(prod[prodList[j]].flatten()))))
            val_y = np.real(1.0j*np.dot(inVal,np.dot(L_rf_y,(prod[prodList[j]].flatten()))))

            if int(val_x) == 2:
                x_rf[i,j]+=1.0
                x_rf[i+63,j+63]+=1.0
            elif int(val_x) == -2:
                x_rf[i,j] -= 1.0
                x_rf[i+63,j+63]-=1.0
            if int(val_y) == 2:
                y_rf[i,j]+=1.0
                y_rf[i+63,j+63]+=1.0
            elif int(val_y) == -2:
                y_rf[i,j] -= 1.0
                y_rf[i+63,j+63]-=1.0

    return x_rf, y_rf

def J_coup_matrix(prod, spin1 = 'I', spin2 = 'R'):
    # prod is dictionary of product operator matrices
    # returns a matrix with 1s in the correct places for
    # adding J couplings between spin1 and spin2
    ham_J = 1.0*(prod[spin1+'z'+spin2+'z'])
    J_mat = np.zeros((127,127),dtype=complex)

    L_J = kron(ham_J,prod['E']) - kron(prod['E'],ham_J)

    for i in range(1,64):
        for j in range(1,64):
            val = np.real(1.0j*np.dot(np.transpose(np.conjugate(prod[prodList[i]]).flatten()),np.dot(L_J,(prod[prodList[j]].flatten()))))
            if int(val)==2:
                J_mat[i,j]+=1.0
                J_mat[i+63,j+63]+=1.0
            elif int(val)==-2:
                J_mat[i,j]-=1.0
                J_mat[i+63,j+63]-=1.0

    return J_mat

def add_J_coup(L_curr, J, J_mat):
    # populate the J matrix with couplings
    Jval = np.pi*J
    a = np.zeros((127,127), dtype = complex)
    a += L_curr
    a += Jval*J_mat

    return a


def L_basic(rate_mat, rate_dic, pb=0.5,kex = 300.0):
    # populate the Liouvillian with relaxation and chemical exchange terms
    a = np.zeros((127,127), dtype=complex)
    # relaxation bits
    pa = 1.0 - pb
    a += rate_mat

    for i, prod in enumerate(prodList):
        if prod=='Iz': # these are for the return of magnetisation to equilibrium
            a[i,0]-=rate_dic['Iz']*pa
            a[i+63,0]-=rate_dic['Iz']*pb
        elif prod=='Rz':
            a[i,0]-=rate_dic['Rz']*pa
            a[i+63,0]-=rate_dic['Rz']*pb
        elif prod=='Sz':
            a[i,0]-=rate_dic['Sz']*pa
            a[i+63,0]-=rate_dic['Sz']*pb

    # add exchange
    k_ge = kex*pb
    k_eg = kex*pa
    k_mat = np.array([[-k_ge, k_eg],[k_ge, -k_eg]])
    k_mat = kron(k_mat, np.eye(63))
    a[1:127,1:127]-=k_mat # signs are important!

    return a


def L_addElements_simp(L_basic, deltaI = 0.0, deltaI_ex = 0.0, deltaR=100.0,
                        deltaS=300.0,deltaR_ex=300, deltaS_ex=100):
    # add chemical shift terms to the Liouvillian
    # Requires that we've already calculated I_cs, I_cs_ex, R_cs ...
    deltaI = 2.0*np.pi*deltaI
    deltaR = 2.0*np.pi*deltaR
    deltaS = 2.0*np.pi*deltaS

    deltaI_ex = 2.0*np.pi*deltaI_ex
    deltaR_ex = 2.0*np.pi*deltaR_ex
    deltaS_ex = 2.0*np.pi*deltaS_ex

    a = np.zeros((127,127), dtype=complex)
    a += deltaI*I_cs
    a += deltaR*R_cs
    a += deltaS*S_cs
    a += deltaI_ex*I_cs_ex
    a += deltaR_ex*R_cs_ex
    a += deltaS_ex*S_cs_ex

    a += L_basic

    return a

def L_add_rf(L_basic, x_rf, y_rf, w1, phase = 0.0):
    # add pulses x_rf and y_rf need to be for the appropriate
    # nuclei. w1 is in Hz and needs to be entered as a numpy array
    # otherwise will crash and burn
    num = len(w1)
    a = np.zeros((num,127,127), dtype=complex)

    w1_x = 2.0*np.pi*w1*np.cos(phase)
    w1_y = 2.0*np.pi*w1*np.sin(phase)

    for i in range(num):
        a[i,:,:] += w1_x[i]*x_rf
        a[i,:,:] += w1_y[i]*y_rf

    a += L_basic

    return a


def tay_val(relaxL,x_rf, y_rf,b1_vals,timey,offsets,deltaN1, deltaN2):

    # let's work out the power and val we need for taylor expansion
    deltaN1_max = np.max(deltaN1-offsets)
    deltaN2_max = np.max(deltaN2-offsets)

    Ltemp = L_addElements_simp(relaxL, deltaR=deltaN1_max, deltaS=deltaN2_max, deltaR_ex = deltaN2_max, deltaS_ex = deltaN1_max)
    CESTp = L_add_rf(Ltemp, x_rf, y_rf, np.array([np.max(b1_vals)]), phase = 0.0)

    curr_L = CESTp[0,:,:]
    normy = linalg.norm(curr_L, ord = 2) # get 2 norm
    scale = 1./normy

    if scale > 1.0:
        time_step = 1.0
        pow = 1.0
        scale = 1.0
    else:
        pow = int(ceil(log(normy,2)))

    scale = 1./(2**pow)
    time_step = timey*scale
    val = -1.0*time_step

    return val, pow


def propagate_fast(init,L,val,pow):
    return np.dot([taylorProp_fast(L[i,:,:], val, pow) for i in range(L.shape[0])],  init)


def taylorProp_fast(L,val,pow):
    x = np.eye(L.shape[0])+val*L + (np.linalg.matrix_power(L,2)*val**2.0)/2.0+(np.linalg.matrix_power(L,3)*val**3.0)/6.0
    xx = np.linalg.matrix_power(x,2**pow)

    return xx


def taylorProp(L, time):
    # Use a Taylor expansion to find the matrix exponential - is slightly
    # less expensive. The time is scaled so that the Taylor series expansion is
    # valid. Can only take one Liouvillian at a time
    normy = linalg.norm(L, ord = 2) # get 2 norm
    scale = 1./normy

    if scale > 1.0:
        time_step = 1.0
        pow = 1.0
        scale = 1.0
    else:
        pow = int(ceil(log(normy,2)))
        scale = 1./(2**pow)
        time_step = time*scale

    val = -1.0*time_step
    x = np.eye(L.shape[0])+val*L + (np.linalg.matrix_power(L,2)*val**2.0)/2.0 +(np.linalg.matrix_power(L,3)*val**3.0)/6.0
    xx = np.linalg.matrix_power(x,2**pow)

    return xx


def propagate(init, L, time):
    # propagate the state either with direct matrix exponential or using
    # Taylor expansion method
    return np.dot([taylorProp(L[i,:,:], time) for i in range(L.shape[0])],  init)
    # return np.dot([linalg.expm(-L[i,:,:]*time) for i in range(L.shape[0])],init)


def propagate_cpd(state_init,cest1,cest2,cpdField,time):
    # propagate with composite decoupling: 90x 240y 90x scheme
    pw90 = 1./(4*cpdField)
    pw240 = pw90*8./3.
    one_lap = 2.0*pw90+pw240

    num_laps = int(ceil(time/one_lap))

    prop1 = np.array([taylorProp(cest1[i,:,:],pw90) for i in range(cest1.shape[0])])
    prop2 = np.array([taylorProp(cest2[i,:,:],pw240) for i in range(cest1.shape[0])])

    prop_main = np.array([np.dot(prop1[i,:,:], np.dot(prop2[i,:,:], prop1[i,:,:])) for i in range(prop1.shape[0])])

    prop_time = np.array([np.linalg.matrix_power(prop_main[i,:,:], num_laps) for i in range(prop_main.shape[0])])

    return np.dot([prop_time[i,:,:] for i in range(prop_main.shape[0])],  state_init)


def propagate_tpit(state_init,cest1,cest2, pulsedField, tau_time, time):
    # propagate CEST with tau-180-tau decoupling
    pw90 = (1./(4*pulsedField))*1.0 # can add in miscalibration here
    one_lap = 2*tau_time
    num_laps = int(ceil(time/(2.0*tau_time)))

    prop1 = np.array([taylorProp(cest1[i,:,:],2.0*pw90) for i in range(cest1.shape[0])]) # 180 deg pulse
    prop2 = np.array([taylorProp(cest2[i,:,:],tau_time-pw90) for i in range(cest1.shape[0])]) # tau time

    prop_main = np.array([np.dot(prop2[i,:,:], np.dot(prop1[i,:,:], prop2[i,:,:])) for i in range(prop1.shape[0])])
    prop_time = np.array([np.linalg.matrix_power(prop_main[i,:,:], num_laps) for i in range(prop_main.shape[0])])

    return np.dot([prop_time[i,:,:] for i in range(prop_main.shape[0])],  state_init)


def getB1_vals(B1_av, B1_std, B1_num):
    # get a Gaussian distribution of B1 values with their associated probabilities
    B1_vals = np.linspace(B1_av - 2.0*B1_std, B1_av + 2.0*B1_std, B1_num)
    B1_probs = []
    for vals in B1_vals:
        B1_probs.append(stats.norm(B1_av, B1_std).pdf(vals))
    # self.B1_probs = numpy.array(self.B1_probs)*(((self.B1_av + 2.0*self.B1_std) - (self.B1_av - 2.0*self.B1_std))/self.B1_num)
    B1_probs = B1_probs/np.sum(B1_probs)

    return np.array(B1_vals), np.array(B1_probs)

### Setup things we need from

prod, prodList = buildProdOperators_IRS(ax,a1,a2,a3,mats)
rel = rel_matrix(prodList)
I_cs, R_cs, S_cs, I_cs_ex, R_cs_ex,S_cs_ex = CS_matrix(prod,prodList)

Ix_rf, Iy_rf = Pulse_matrix('I')
Rx_rf, Ry_rf = Pulse_matrix('R')
Sx_rf, Sy_rf = Pulse_matrix('S')

RS_x_rf = Rx_rf + Sx_rf
RS_y_rf = Ry_rf + Sy_rf

J_IR_mat = J_coup_matrix(prod,spin1 = 'I',spin2 = 'R')
J_IS_mat = J_coup_matrix(prod,spin1 = 'I',spin2 = 'S')
