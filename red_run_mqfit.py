from ABX_reducedBasis import *
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report
import matplotlib.pyplot as plt
from matplotlib import cm
import os,sys

dRateMat = {'r_Cz':0.5, 'r_Nz':1.65, 'r_Nxy':16.0}
dGlobParams = {'J':20.0, 'kex':50.0, 'pb':0.5, 'deltaO':1.} # these do not change with field
# note that in this case w0 is the DQ chemical shift in ppm and deltaO is in pppm
# these can be fitting parameters or not we state what parameters we want to be
# fit below

dFitParams = ['kex','w0','deltaO']
dFitParams_rel = ['r_Nz','r_Nxy']
dExpParams = {'cest_time':0.25,'inhom_num':10,'phase':0.0}

# inhom_num = 5 # number of points to account for inhomogeneity

class fitGroup(object):
    def __init__(self, inputFile, outfile, fitParams=dFitParams, fitParams_rel = dFitParams_rel,
            dExpParams = dExpParams, dGlobParams = dGlobParams, dRateMat = dRateMat):
    # fields is list of 15N carrier frequencies; infiles are corresponding CEST
    # v1 are the B1 field strength in Hz we can have multiple B1 field strengths
    # for each static field but we should repeat fields when this is the case.
    # v1std is the corresponding standard deviation in v1 (Gaussian distribution)
        self.dGlobParams = dGlobParams
        self.dExpParams = dExpParams
        self.dRateMat = dRateMat
        self.readInput(inputFile)
        self.start = time.time()
        self.offsets = {} # dictionary of offsets: 1 key per field/v1
        self.intens = {} # dictionary of cest intensities: 1 key per field/v1
        self.intens_err = {}

        self.v1 = {}
        self.v1_probs = {}

        self.pathname = outfile+'/'+self.name
        if os.path.exists(self.pathname)==0:
            os.mkdir(self.pathname)
        else:
            print('folder {} already exists. Contents may be overwritten'.format(self.pathname))

        for i,lb in enumerate(self.label):
            self.offsets[lb],self.intens[lb],self.intens_err[lb]=self.readInfile(self.infiles[i],self.minerr[lb])
            # populate our offset and intensity dictionaries

        self.fitParams = fitParams
        self.fitParams_rel = fitParams_rel

        for i,lb in enumerate(self.label):
            self.v1[lb], self.v1_probs[lb] = getB1_vals(self.v1_list[i],self.v1_std_list[i],self.inhom_num)

        self.doFit()
        self.end = time.time()
        print('\nCompleted fit for ', self.name)
        print('time taken: ', self.end-self.start, '\n')

    def readInfile(self, infile, minerr): # get real/simulated CEST curve
        offset = []
        cesti = []
        cesterr = []
        with open(infile,'r') as inny:
            for line in inny:
                line = line.split()
                if len(line)>1:
                    if line[0]!='#':
                        offset.append(float(line[0]))
                        cesti.append(float(line[1]))
                        # cesterr.append(float(line[2]))
                        cesterrVal = np.fabs(float(line[2]))
                        if cesterrVal < minerr:
                            cesterr.append(minerr)
                        else:
                            cesterr.append(cesterrVal)

        return np.array(offset), np.array(cesti), np.array(cesterr) # return offsets and cesti

    def readInput(self,infile): # read residue specic input file
        allDat = {}
        self.fields = []
        self.v1_list = []
        self.v1_std_list = []
        self.minerr = {}
        self.w0 = {}
        self.carrier = {}
        self.rel = {}
        self.glob = {}

        with open(infile,'r') as inny:
            for line in inny:
                if not line.startswith('#'):
                    if len(line)>0:
                        line = line.split()
                        if str.lower(line[0])=='resi':
                            allDat['resi']=line[1]
                        elif str.lower(line[0])=='datafiles':
                            allDat['datafiles']=list(line[1:])
                        else:
                            allDat[line[0]]=list(map(float,line[1:]))

        try:
            self.name = allDat['resi']
            print('loading in parameters for ', self.name)
        except:
            print('need a residue name aborting now...')
            sys.exit()


        try:
            self.infiles = allDat['datafiles']
        except:
            print('did not find any cest data files aborting now...')
            sys.exit()

        try:
            self.fields = allDat['field']
            print('15N spectrometer fields are (MHz) ', self.fields)
        except:
            print('did not find any fields in input file aborting now...')
            sys.exit()

        try:
            self.v1_list = allDat['v1']
            print('CEST B1 field strengths are (Hz)', self.v1_list)
        except:
            print('did not find any CEST fields in input file aborting now...')
            sys.exit()

        try:
            self.v1_std_list = allDat['v1_std']

        except:
            print('did not find any SD for CEST field strengths going to default values (10% SD)')
            self.v1_std_list = [self.v1_list[i]*0.1 for i in range(len(self.v1_std_list))]
        print('CEST B1 field std are (Hz) ', self.v1_std_list)

        self.label = list(zip(self.fields,self.v1_list))

        assert len(self.infiles)==len(self.fields)==len(self.v1_list)==len(self.v1_std_list), \
                "Require number of data files, spectrometer fields, B1 fields and B1 field SD to be equal, aborting now..."
        self.num = len(self.fields)

        try:
            # print('minerr values are ',end=' ')
            for i,label in enumerate(self.label):
                self.minerr[label] = (allDat['minerr'])[i]
                # print(self.minerr[label],end=' ')
            print('')
        except:
            for i,label in enumerate(self.label):
                self.minerr[label] = 0.0
            print('not using minerr (setting to 0)')


        try:
            # print('carrier frequencies are (ppm)',end=' ')
            for i,label in enumerate(self.label):
                self.carrier[label] = (allDat['carrier'])[i]
                # print(self.carrier[label],end=' ')
            print('')
        except:
            print('did not find right number carrier frequencies in input file aborting now...')
            sys.exit()

        assert len(self.carrier)==self.num, "Require carrier frequency for every spectrometer field..."

        try:
            self.glob['w0'] = allDat['w0'][0]
            print('15N SQ resonance frequency is (ppm)', self.glob['w0'])
        except:
            print('did not find 15N resonance frequency in input file aborting now...')
            sys.exit()

        try:
            self.glob['J'] = allDat['J'][0]
            print('J (Hz) = ', self.glob['J'])
        except:
            print('did not find J coupling constant in input file using default value of ', self.dGlobParams['J'])
            self.glob['J'] = self.dGlobParams['J']

        try:
            self.glob['kex'] = allDat['kex'][0]
            print('kex initial (s^-1) = ', self.glob['kex'])
        except:
            print('did not find kex value in input file using default value of ', self.dGlobParams['kex'])
            self.glob['kex'] = self.dGlobParams['kex']

        try:
            self.glob['pb'] = allDat['pb'][0]
            print('pb = ', self.glob['pb'])
        except:
            print('did not find pb value in input file using default value of ', self.dGlobParams['pb'])
            self.glob['pb'] = self.dGlobParams['pb']

        try:
            self.glob['deltaO'] = allDat['deltaO'][0]
            print('deltaO initial = ', self.glob['deltaO'])
        except:
            print('did not find deltaO value in input file using default value of ', self.dGlobParams['deltaO'])
            self.glob['deltaO'] = self.dGlobParams['deltaO']

        try:
            self.cest_time = allDat['cest_time'][0]
            print('cest saturation time (s) = ', self.cest_time)
        except:
            print('did not find cest_time value in input file using default value of ', self.dExpParams['cest_time'])
            self.cest_time = self.dExpParams['cest_time']

        try:
            self.phase = allDat['phase'][0]
            print('cest pulse phase (in degrees) ', self.phase)
        except:
            print('did not find cest pulse phase in input file using default value of ', self.dExpParams['phase'])
            self.phase = self.dExpParams['phase']

        try:
            self.inhom_num = int(allDat['inhom'][0])
            print('number of B1 inhomogeneity points to use in calculation ', self.inhom_num)
        except:
            print('did not find number of inhomogeneity B1 points to trial in calculation using default value of ', self.dExpParams['inhom_num'])
            self.inhom_num = self.dExpParams['inhom_num']

        for key in self.dRateMat:
            for i,field in enumerate(self.fields):
                rfield = str(field).split('.')[0]
                try:
                    # self.rel[key+'_'+str(lb)+'_'+str(self.v1_list[i])] = allDat[key][i]
                    self.rel[key+'_'+rfield] = allDat[key][i]
                except:
                    # self.rel[key+'_'+str(lb)+'_'+str(self.v1_list[i])] = dRateMat[key]
                    self.rel[key+'_'+rfield] = self.dRateMat[key]


    def doFit(self):

        params_all = Parameters()
        meth = 'leastsq'
        for param in self.glob:
            # add in Global parameters
            params_all.add(param, value = self.glob[param], min = 0.0, vary = False)

        for param in self.rel:
            # add in relaxation terms
            params_all.add(str(param), value = self.rel[param], min = 0.0, vary = False)

        for param in self.fitParams:
            params_all[param].vary = True

        for param in self.fitParams_rel:
            for field in self.fields:
                rfield = str(field).split('.')[0]
                params_all[param+'_'+rfield].vary = True

        print('\nstarting minimisation for ', self.name, '...\n')

        minner = Minimizer(self.chiFunc, params_all)
        result = minner.minimize(method = meth)

        final_rep = fit_report(result)
        print(final_rep)
        with open(self.pathname+'/'+self.name+'_fitReport.out','w') as outy:
                outy.write('%s' %(final_rep))

        for lb in self.label:
            final_res = self.sim_cest(result.params,lb)
            ppm = self.offsets[lb]/lb[0]+self.carrier[lb]
            plt.plot(ppm, self.intens[lb], 'o')
            rfield = str(lb[0]).split('.')[0]
            rB1 = str(lb[1]).split('.')[0]
            plt.plot(ppm, final_res, '-')
            plt.xlim([np.max(ppm),np.min(ppm)])
            plt.xlabel(r'$^{15}$N (ppm)')
            plt.ylabel(r'I/I$_0$')
            plt.savefig(self.pathname+'/'+self.name+'_'+str(rfield)+'MHz_'+rB1+'Hz.eps')
            plt.clf()

            with open(self.pathname+'/'+self.name+'_'+str(rfield)+'MHz_'+rB1+'Hz.out','w') as outy:
                outy.write('# ppm\tExp\tFit\n')
                for i,val in enumerate(ppm):
                    outy.write('%f\t%f\t%f\n' %(val,self.intens[lb][i],final_res[i]))
            # plt.show()

    def chiFunc(self, params):

        for i,lb in enumerate(self.label):
            sim_res = self.sim_cest(params,lb)
            if i == 0:
                res = ((self.intens[lb] - sim_res)/self.intens_err[lb])
            else:
                res_curr = ((self.intens[lb] - sim_res)/self.intens_err[lb])
                res = np.concatenate((res,res_curr))

        return res


    def sim_cest(self,params,label):
        kex = params['kex'].value
        pb = params['pb'].value
        deltaO = params['deltaO'].value
        J = params['J'].value
        w0 = params['w0'].value

        field = label[0]
        rfield = str(field).split('.')[0]

        v1 = self.v1[label]
        v1_probs = self.v1_probs[label]
        offsets = self.offsets[label]
        carrier = self.carrier[label]


        deltaN1 = (w0 + 0.5*deltaO - carrier)*field # convert to Hz
        deltaN2 = (w0 - 0.5*deltaO - carrier)*field # convert to Hz
        cest_fin = np.zeros(len(offsets))
        state = np.zeros((1+2*rowLen),dtype="complex64")

        rate_dic = {}
        cz = params['r_Cz'+'_'+rfield].value
        nz = params['r_Nz'+'_'+rfield].value
        nxy = params['r_Nxy'+'_'+rfield].value
        rate_dic['Iz'] = cz
        rate_dic['Rz'] = nz
        rate_dic['Sz'] = nz
        rate_dic['RxSx'] = 2.0*nxy
        rate_dic['RxSy'] = 2.0*nxy
        rate_dic['RySx'] = 2.0*nxy
        rate_dic['RySy'] = 2.0*nxy
        rate_dic['RxSz'] = nxy + nz
        rate_dic['RySz'] = nxy + nz
        rate_dic['RzSx'] = nxy + nz
        rate_dic['RzSy'] = nxy + nz
        rate_dic['RzSz'] = 2.0*nz
        rate_dic['IzRxSx'] = cz+(2.0*nxy)
        rate_dic['IzRxSy'] = cz+(2.0*nxy)
        rate_dic['IzRySx'] = cz+(2.0*nxy)
        rate_dic['IzRySy'] = cz+(2.0*nxy)
        rate_dic['IzRxSz'] = cz+nxy+nz
        rate_dic['IzRySz'] = cz+nxy+nz
        rate_dic['IzRzSx'] = cz+nxy+nz
        rate_dic['IzRzSy'] = cz+nxy+nz
        rate_dic['IzRzSz'] = (2.0*nz)+cz

        state[0] += 1.0
        state[rowLen] += 0.5 # start on IzSz ground
        state[2*rowLen] += 0.5 # split between ground and excited evenly
        rate_mat = pop_rel_red(rel,rate_dic,rowLen)

        relaxL = L_basic_red(rate_mat,rate_dic,rowList,rowLen, pb=pb, kex = kex)
        relaxL = add_J_coup_red(relaxL, J, J_IR_mat,rowLen) # add IR J-coupling
        relaxL = add_J_coup_red(relaxL, J, J_IS_mat,rowLen)
        val, pow = tay_val_red(relaxL,rowLen,RS_x_rf, RS_y_rf,v1,self.cest_time,offsets,deltaN1,deltaN2)

        for k,offset in enumerate(offsets):
            Ltemp = L_addElements_simp_red(relaxL,rowLen, deltaR=deltaN1-offset, deltaS=deltaN2-offset, deltaR_ex = deltaN2-offset, deltaS_ex = deltaN1 - offset)
            CESTp = L_add_rf_red(Ltemp, rowLen, RS_x_rf, RS_y_rf, v1, phase = 0.0)
            state_curr = propagate_fast(state,CESTp,val,pow)
            #state_curr = propagate(state, CESTp, self.cest_time)
            cest_fin[k] = np.sum(np.real(state_curr[:,rowLen] + state_curr[:,2*rowLen])*v1_probs)

        return cest_fin
