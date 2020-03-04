#!/usr/bin/python3

from red_run_fit import fitGroup
import sys, os

def getStarted(infile):
    with open(infile,'r') as inny:
        dRateMat = {}
        dGlobParams = {}
        dExpParams = {}
        for line in inny:
            if not line.startswith('#') and len(line.split())>0:
                line = line.split()
                if str.lower(line[0])=='infol':
                    try:
                        infol = line[1]
                    except:
                        infol = './'
                        print('assuming input files are in current directory')
                if str.lower(line[0])=='outfol':
                    try:
                        outfol = line[1]
                    except:
                        outfol = './'
                        print('outputs will be placed in current directory')
                if str.lower(line[0])=='resilist':
                    resiList = list(line[1:])
                if str.lower(line[0])=='fitparams':
                    fitParams = list(line[1:])
                if str.lower(line[0])=='fitparams_rel':
                    fitParams_rel = list(line[1:])
                if str.lower(line[0])=='r_cz':
                    dRateMat['r_Cz'] = float(line[1])
                if str.lower(line[0])=='r_nz':
                    dRateMat['r_Nz'] = float(line[1])
                if str.lower(line[0])=='r_cz':
                    dRateMat['r_Nxy'] = float(line[1])
                if str.lower(line[0])=='j':
                    dGlobParams['J'] = float(line[1])
                if str.lower(line[0])=='j':
                    dGlobParams['J'] = float(line[1])
                if str.lower(line[0])=='kex':
                    dGlobParams['kex'] = float(line[1])
                if str.lower(line[0])=='pb':
                    dGlobParams['pb'] = float(line[1])
                if str.lower(line[0])=='deltao':
                    dGlobParams['deltaO'] = float(line[1])
                if str.lower(line[0])=='cest_time':
                    dExpParams['cest_time'] = float(line[1])
                if str.lower(line[0])=='inhom_num':
                    dExpParams['inhom_num'] = float(line[1])
                if str.lower(line[0])=='phase':
                    dExpParams['phase'] = float(line[1])

    return infol,outfol,resiList,fitParams,fitParams_rel,dRateMat,dExpParams,dGlobParams

def doFit():
    args = sys.argv
    print('\nMQ-CEST data fitting\n')
    print('Gogulan Karunanithy, UCL, 2020\n')
    if len(args)>2:
        print('Only expecting one argument additional arguments will be ignored...')
    if len(args)==1:
        print('Require an input file to start fits. See example...')
        sys.exit()

    infile = args[1]
    infol,outfol,resiList,fitParams,fitParams_rel,dRateMat,dExpParams, dGlobParams = getStarted(infile)
    if os.path.exists(outfol)==0:
        os.mkdir(outfol)
    for resi in resiList:
        fitGroup(infol+'/'+resi+'.in',outfol,fitParams=fitParams,
            fitParams_rel = fitParams_rel, dExpParams = dExpParams,
            dGlobParams = dGlobParams, dRateMat = dRateMat)

doFit()
