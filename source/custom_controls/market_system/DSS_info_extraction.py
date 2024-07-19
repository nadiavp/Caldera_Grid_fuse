#2023.07.21 Modified by Mingzhi to support the three phase load extraction (original code will report error when there are three phase load in the circuit)
# 2020.09.09 Add the codes for load info extraction, be careful with the three phase connected loads, delata and wye loads are single phase load, for the delta load, a-b connected, the phase is set as 0, b-c connected, phase is set as 1.
# The main file of interest is Network_3ph_pf.py - this is used to create a Network_3ph object. The Network_3ph_pf.setup_network_ieee13() sets up a network object as the ieee13 bus. As discussed, the main idea would be to have an alternative function which could replace this method, to load in bus & line data from an external file.

# Two main panda dataframes need to be set:
# bus_df, with columns ['name','number','load_type','connect','Pa','Pb','Pc','Qa','Qb','Qc’], and 
# line_df with columns ['busA','busB','Zaa','Zbb','Zcc','Zab','Zac','Zbc','Baa','Bbb','Bcc','Bab','Bac','Bbc’]

# The file zbus_3ph_pf_test.py can be used for testing the Network_3ph class.


# ADVISED ACTIONS prior to conversion =======================
# 1. Transformers at the head of the feeder to 1:1 at the secondary voltage
# 2. Transformer source impedance based on xfmr ratio 
# - Have this as a 'source' line impedance as in ieee_tn (from previous ptech18 folder)
# If there are downstream transformers:
# 3. Convert to 1:1 on primary
# 4. Convert downstream impedances and voltages
# ============================================================
# WARNING: be very careful about making sure scripts set the default base frequency if there are capacitance values in the circuit that are not trivial; particularly when changing from script to script (OpenDSS uses a global base frequency which changes when you run, for example, the EU LV network, and it is not trivial to keep track of what it is).

import sys, os
# path = os.path.dirname(os.path.dirname(sys.argv[0]))

import pandas as pd
#import win32com.client
import numpy as np
import opendssdirect as dss
from scipy.sparse import csc_matrix
#import matplotlib.pyplot as plt
#import matplotlib
#from Network_3ph_pf import Network_3ph
from .miscDssFncs import get_ckt, tp2mat, tp_2_ar, get_sYsD

def dssConvert(feeder,dssName, save_model=True):
    path = os.path.dirname(os.path.dirname(sys.argv[0]))
    if os.path.exists(os.path.join(dssName)):
        fn_ckt = os.path.join(dssName)
        dir0 = os.path.dirname(dssName)
        sn0 = os.path.dirname(dssName)
    elif feeder[0]!='n':
        ntwxs = os.path.join(path,'Data','networks')
        fn_ckt = os.path.join(ntwxs,'opendssNetworks',feeder,dssName)
        dir0 = os.path.join(ntwxs,feeder+'_dss')
        sn0 = os.path.join(dir0,feeder)
    else:
        ntwk_no = feeder.split('_')[0][1:]
        fdr_no = feeder.split('_')[1]
        ntwxs = os.path.join(path,'Data','networks')
        fn_ckt = os.path.join(ntwxs,'opendssNetworks','network_'+ntwk_no,'Feeder_'+fdr_no,dssName)
        dir0 = os.path.join(ntwxs,feeder+'_dss')
        sn0 = os.path.join(dir0,feeder)
    if os.path.exists(os.path.dirname(fn_ckt)):
        print('Start extract info from feeder: ', feeder)
        #saveModel = True
        #saveModel = False

        dss.Command(f'Redirect {fn_ckt}') #DSSObj = win32com.client.Dispatch("OpenDSSEngine.DSS")
        dss.Command('Set Controlmode=off') # turn off regs/cap banks
        dss.Solution.Solve()#dss.Command('Solve')
        print('opendss engine started')
        #DSSText= DSSObj.Text
        #DSSCircuit = DSSObj.ActiveCircuit
        #DSSSolution=DSSCircuit.Solution

        #DSSText.Command = 'Compile ('+fn_ckt+ ')'
        #DSSText.Command='Compile ('+fn_ckt+'.dss)'
        
        #DSSText.Command='Set Controlmode=off' # turn all regs/caps off
        #DSSText.Command='solve'
        

        LDS = dss.Loads#DSSCircuit.Loads
        CAP = dss.Capacitors#DSSCircuit.Capacitors
        LNS = dss.Lines#DSSCircuit.Lines
        TRN = dss.Transformers#DSSCircuit.Transformers
        ACE = dss.CktElement#DSSCircuit.ActiveElement
        SetACE = dss.Circuit.SetActiveElement#DSSCircuit.SetActiveElement
        #SetABE = dss.Circuit.SetActiveBus#DSSCircuit.SetActiveBus
        #ADE = dss.Circuit.ActiveDSSElement#DSSCircuit.ActiveDSSElement
        #ABE = dss.ActiveBus#DSSCircuit.ActiveBus
        SRC = dss.Vsources#DSSCircuit.Vsources

        # create solution DF
        slnColumns = ['bus','vLN','vPu','sInjkW']     
        bus_voltage_magni_pu = np.array(dss.Circuit.AllBusMagPu())
        YNodeV = tp_2_ar(dss.Circuit.YNodeVArray())  # complex tuples to array. Node voltage in sequence of Y node 
        sY,sD,iY,iD,yzD,iTot,H = get_sYsD(dss.Circuit) # sY sD are Y/Delta connected power injections, iY,iD are Y/Delta connected current injections
        sInj = 1e-3*YNodeV*(iTot.conj())
        sList = list(sInj)  # Node power injection in Y node order
        yList = list(dss.Circuit.YNodeOrder()) # Node sequence in Y node order
        vList = tp_2_ar(dss.Circuit.YNodeVArray()).tolist()
        yvList = list(map(list, zip(*[yList,vList,bus_voltage_magni_pu,sList]))) # list of list: [node name, node voltage(complex number format), node power injection]
        solution_df = pd.DataFrame(data=yvList,columns=slnColumns)  # node name, node volatge, node power injection (in Y node order)           

        # create source DF
        SRC.First()
        srcColumns = ['kvBaseLL','pu']
        kvbase = SRC.BasekV()
        print(f'voltage source name {SRC.Name()}')
        src_df = pd.DataFrame(data=[[kvbase,SRC.PU()]],columns=srcColumns,index=[SRC.Name()])

        nLds = LDS.Count()
        nCap = CAP.Count()
        nTrn = TRN.Count()
        nLns = LNS.Count()

        lineColumns = ['busA','busB','Zaa','Zbb','Zcc','Zab','Zac','Zbc','Baa','Bbb','Bcc','Bab','Bac','Bbc']
        line_df = pd.DataFrame(data=np.zeros((nLns,len(lineColumns)),dtype=complex), index=LNS.AllNames(), columns=lineColumns)

        # Create Line DF from lines and transformers ====
        i = LNS.First()
        Yprm = {}
        Yprm0 = {}
        YprmErrDiag = []

        while i:
            lineName = LNS.Name()
            line_df.loc[lineName,'busA']=LNS.Bus1().split('.')[0]
            line_df.loc[lineName,'busB']=LNS.Bus2().split('.')[0]
            lineLen = LNS.Length()
            if len(LNS.Geometry())>0:
                zmat0 = (tp2mat(LNS.RMatrix()) + 1j*tp2mat(LNS.XMatrix())) # ohms
                bmat0 = 1j*2*np.pi*60*tp2mat(LNS.CMatrix())*1e-9 # ohms
            else:
                zmat0 = (tp2mat(LNS.RMatrix()) + 1j*tp2mat(LNS.XMatrix()))*lineLen # ohms
                bmat0 = 1j*2*np.pi*60*tp2mat(LNS.CMatrix())*1e-9*lineLen # ohms
            
            SetACE('Line.'+LNS.Name())
            
            nPh = ACE.NumPhases()
            phs = list(map(int,ACE.BusNames()[0].split('.')[1:]))
            phsIdx = np.array(phs)-1
            
            Zmat = np.zeros((3,3),dtype=complex)
            Bmat = np.zeros((3,3),dtype=complex)
            
            Yprm0_Y = tp_2_ar(LNS.Yprim())
            Yprm0[lineName] = Yprm0_Y.reshape((np.sqrt(Yprm0_Y.shape)[0].astype('int32'),np.sqrt(Yprm0_Y.shape)[0].astype('int32')))[0:nPh,0:nPh]
            
            if nPh==1:
                Zmat[phs[0]-1,phs[0]-1] = zmat0[0,0]
                Bmat[phs[0]-1,phs[0]-1] = bmat0[0,0]
                Yprm[lineName] = 1/Zmat[phsIdx,phsIdx]
            if nPh==2:
                Zmat[phs[0]-1,phs[0]-1] = zmat0[0,0]
                Zmat[phs[1]-1,phs[0]-1] = zmat0[0,1]
                Zmat[phs[0]-1,phs[1]-1] = zmat0[0,1]
                Zmat[phs[1]-1,phs[1]-1] = zmat0[1,1]
                Yprm[lineName] = np.linalg.inv(Zmat[phsIdx][:,phsIdx])
                
                Bmat[phs[0]-1,phs[0]-1] = bmat0[0,0]
                Bmat[phs[1]-1,phs[0]-1] = bmat0[0,1]
                Bmat[phs[0]-1,phs[1]-1] = bmat0[0,1]
                Bmat[phs[1]-1,phs[1]-1] = bmat0[1,1]
            if nPh==3:
                Zmat = zmat0
                Bmat = bmat0
                Yprm[lineName] = np.linalg.inv(Zmat)
            
            line_df.loc[lineName,'Zaa']=Zmat[0,0]
            line_df.loc[lineName,'Zbb']=Zmat[1,1]
            line_df.loc[lineName,'Zcc']=Zmat[2,2]
            line_df.loc[lineName,'Zab']=Zmat[0,1]
            line_df.loc[lineName,'Zac']=Zmat[0,2]
            line_df.loc[lineName,'Zbc']=Zmat[1,2]
            
            line_df.loc[lineName,'Baa']=Bmat[0,0]
            line_df.loc[lineName,'Bbb']=Bmat[1,1]
            line_df.loc[lineName,'Bcc']=Bmat[2,2]
            line_df.loc[lineName,'Bab']=Bmat[0,1]
            line_df.loc[lineName,'Bac']=Bmat[0,2]
            line_df.loc[lineName,'Bbc']=Bmat[1,2]
            
            Yprm0_Y = tp_2_ar(LNS.Yprim())
            Yprm0[lineName] = Yprm0_Y.reshape((np.sqrt(Yprm0_Y.shape)[0].astype('int32'),np.sqrt(Yprm0_Y.shape)[0].astype('int32')))[0:nPh,0:nPh]
            
            YprmDiag = np.linalg.inv(zmat0) + 0.5*bmat0
            YprmErrDiag= YprmErrDiag + [np.linalg.norm(Yprm0[lineName] - YprmDiag)/np.linalg.norm(Yprm0[lineName])] # for validation/checking
            i = LNS.Next()

        # print('Primitive impedance errors:',YprmErrDiag) # error checking


        if nTrn > 0:
            trnColumns = ['busA','busB','typeA','typeB','Zseries','Zshunt','kV_A','kV_B']
            trn_df = pd.DataFrame(data=np.zeros((nTrn,len(trnColumns)),dtype=complex), index=TRN.AllNames(), columns=trnColumns)

            # self.transformer_df.append({'busA':'633','busB':'634','typeA':'wye-g','typeB':’wye-g','Zseries':0.381+0.692j,'Zshunt':0},ignore_index=True)  
            # The connection types (‘typeA’ and ‘typeB’) can be 'wye-g', ‘wye' or ‘delta’.

            # Now: go through each of the transformers and use to create the transformer dataframe.

            i = TRN.First() # only Lines implemented at this stage.
            while i:
                trnName = TRN.Name()
                SetACE('Transformer.'+TRN.Name())  # Set the element as active for info extraction
                
                # Winding 1
                TRN.Wdg(1)
                #print(f'buses: {ACE.BusNames()} connected to transformer {trnName}')
                trn_df.loc[trnName,'busA'] = ACE.BusNames()[0]#ADE.Properties('bus').Val
                trn_df.loc[trnName,'kV_A'] = TRN.kV() #ADE.Properties('kV').Val
                if TRN.IsDelta():
                    conn1 = 'Delta'#ADE.Properties('conn').Val
                else:
                    conn1 = 'wye'
                grnded1 = TRN.Rneut()#float(ADE.Properties('Rneut').Val)
                kV1 = TRN.kV()#float(ADE.Properties('kV').Val)
                S1 = float(TRN.kVA())#float(ADE.Properties('kva').Val)
                zbase1 = (kV1**2)/(S1*1e-3)
                nodelosses = TRN.LossesByType() # total, load, no load
                zSh1Pct = nodelosses[4] + 1j*nodelosses[5]#float(ADE.Properties('%Noloadloss').Val) + 1j*float(ADE.Properties('%imag').Val)
                zSr1Pct =TRN.R() + 0.5*1j*TRN.Xhl()#float(ADE.Properties('%r').Val) + 0.5*1j*float(ADE.Properties('xhl').Val)
                
                # find how many nodes there are at the connected bus
                #SetABE(TRN.BusNames())#ADE.Properties('bus').Val)
                numNodes = len(ACE.BusNames())#ABE.Nodes)
                
                if conn1=='wye ' and grnded1==-1 and numNodes==4:
                    trn_df.loc[trnName,'typeA'] = 'wye'
                elif conn1=='wye ':
                    trn_df.loc[trnName,'typeA'] = 'wye-g'
                elif conn1=='Delta ': # Be careful, the connection type returned by OpenDSS is 'Delta '
                    trn_df.loc[trnName,'typeA'] = 'delta'
                
                #ADE.Properties('wdg').Val  # 似乎没有什么用处！
                
                # Winding 2
                TRN.Wdg(2)
                trn_df.loc[trnName,'busB'] = ACE.BusNames()[1] #ADE.Properties('bus').Val
                trn_df.loc[trnName,'kV_B'] = TRN.kV() #ADE.Properties('kV').Val
                is_delta = TRN.IsDelta()
                if is_delta:
                    conn2='Delta'
                else:
                    conn2 = 'wye'#ADE.Properties('conn').Val
                #grnded2 = float(ADE.Properties('Rneut').Val)
                kV2 = TRN.kV()#float(ADE.Properties('kV').Val)
                S2 = TRN.kVA()#float(ADE.Properties('kva').Val)
                zbase2 = (kV2**2)/(S2*1e-3)
                zSh2Pct = TRN.LossesByType()[4] + 1j*TRN.LossesByType()[4]#float(ADE.Properties('%Noloadloss').Val) + 1j*float(ADE.Properties('%imag').Val)
                zSr2Pct = TRN.R() + 0.5*1j*TRN.Xhl()#float(ADE.Properties('%r').Val) + 0.5*1j*float(ADE.Properties('xhl').Val) # only half of leakage
                
                if conn1=='wye ' and grnded1==-1 and numNodes==4:
                    trn_df.loc[trnName,'typeB'] = 'wye'
                elif conn2=='wye ':
                    trn_df.loc[trnName,'typeB'] = 'wye-g'
                elif conn2=='Delta ':
                    trn_df.loc[trnName,'typeB'] = 'delta'
                
                Zsr = 0.01*(zbase1*zSr1Pct + zbase2*zSr2Pct)
                Zsh = 0.01*(zbase1*zSh1Pct + zbase2*zSh2Pct)
                
                trn_df.loc[trnName,'Zseries'] = Zsr
                if Zsh==0:
                    trn_df.loc[trnName,'Zshunt'] = 0
                
                i = TRN.Next()


        # create bus_df from loads and capacitors ======
        nBus = len(dss.Circuit.AllBusNames())
        nLoads = len(LDS.AllNames())

        Load_names= LDS.AllNames()
        
        busColumns = ["name","number","v_base","load_type","connect","Pa","Pb","Pc","Qa","Qb","Qc"]
        bus_df = pd.DataFrame(data=np.zeros((nBus,len(busColumns))), index=dss.Circuit.AllBusNames(), columns=busColumns)
        #bus_df = pd.DataFrame(index=DSSCircuit.AllBusNames, columns=busColumns)
        #
        # extract the voltage base of each bus
        BusNameList = dss.Circuit.AllBusNames()
        bus_voltage_base=[]
        for bus_name in BusNameList:
            #SetABE(bus_name)
            dss.Circuit.SetActiveBus(bus_name)
            bus_voltage_base.append(dss.Bus.kVBase())#ABE.kVBase*1000)
        
        bus_df.loc[:,'v_base'] = bus_voltage_base        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        LoadColumns = ["name","number","bus_name","connect","P","Q","phase_number","phase"]
        #load_df = pd.DataFrame(data=np.zeros((nLoads,len(LoadColumns))), index=Load_names, columns=LoadColumns)
        load_df = pd.DataFrame(index=Load_names, columns=LoadColumns)

        bus_df['name'] = dss.Circuit.AllBusNames()
        bus_df['number'] = np.arange((nBus))

        load_df['name'] = Load_names
        load_df['number'] = np.arange((nLoads))

        # Find the slack bus:
        VSRC = dss.Vsources
        VSRC.First()
        SetACE('Vsource.'+VSRC.Name())
        bus_df.loc[:,'connect'] = 'Y'
        bus_df.loc[ACE.BusNames()[0],'load_type'] = 'S'
        

        i = LDS.First()

        while i:
            SetACE('Load.'+LDS.Name())
            actBus = ACE.BusNames()[0].split('.')[0]
            actLoad = LDS.Name()
            load_df.loc[actLoad,'bus_name'] = actBus
            
            if LDS.Model()==1:
                load_type = 'PQ'
            elif LDS.Model()==2:
                load_type = 'Z'
            elif LDS.Model()==5:
                load_type = 'I'
            else:
                print('Warning! Load: ',LDS.Name(),'Load model not a ZIP load. Setting as PQ.')
                load_type = 'PQ'
            
            if bus_df.loc[actBus,'load_type']==0 or bus_df.loc[actBus,'load_type']==load_type:
                bus_df.loc[actBus,'load_type'] = load_type
            else:
                bus_df.loc[actBus,'load_type'] = 'Mxd'

            nPh = ACE.NumPhases()
            phs = ACE.BusNames()[0].split('.')[1:]
            load_df.at[actLoad,'phase_number'] = nPh            
            # load_df.at[actLoad,'P'] = LDS.kW + load_df.loc[actLoad,'P']
            # load_df.at[actLoad,'Q'] = LDS.kvar + load_df.loc[actLoad,'Q']      
            load_df.at[actLoad,'P'] = LDS.kW() 
            load_df.at[actLoad,'Q'] = LDS.kvar() 
            
            
            if LDS.IsDelta():
                bus_df.loc[actBus,'connect'] = 'D'
                load_df.loc[actLoad,'connect'] = 'D'
                if nPh==1:
                    if '1' in phs and '2' in phs:
                        bus_df.loc[actBus,'Pa'] = LDS.kW() + bus_df.loc[actBus,'Pa']
                        bus_df.loc[actBus,'Qa'] = LDS.kvar() + bus_df.loc[actBus,'Qa']
                        load_df.loc[actLoad,'phase'] = [0]
                    if '2' in phs and '3' in phs:
                        bus_df.loc[actBus,'Pb'] = LDS.kW() + bus_df.loc[actBus,'Pb']
                        bus_df.loc[actBus,'Qb'] = LDS.kvar() + bus_df.loc[actBus,'Qb']
                        load_df.loc[actLoad,'phase'] = [1]
                    if '3' in phs and '1' in phs:
                        bus_df.loc[actBus,'Pc'] = LDS.kW() + bus_df.loc[actBus,'Pc']
                        bus_df.loc[actBus,'Qc'] = LDS.kvar() + bus_df.loc[actBus,'Qc']
                        load_df.loc[actLoad,'phase'] = [2]
                if nPh==3:
                    bus_df.loc[actBus,'Pa'] = LDS.kW()/3 + bus_df.loc[actBus,'Pa']
                    bus_df.loc[actBus,'Pb'] = LDS.kW()/3 + bus_df.loc[actBus,'Pb']
                    bus_df.loc[actBus,'Pc'] = LDS.kW()/3 + bus_df.loc[actBus,'Pc']
                    bus_df.loc[actBus,'Qa'] = LDS.kvar()/3 + bus_df.loc[actBus,'Qa']
                    bus_df.loc[actBus,'Qb'] = LDS.kvar()/3 + bus_df.loc[actBus,'Qb']
                    bus_df.loc[actBus,'Qc'] = LDS.kvar()/3 + bus_df.loc[actBus,'Qc']
                    load_df.at[actLoad,'phase'] = [0, 1, 2]	
                    
                if nPh==2:
                    print('Warning! Load: ',LDS.Name(),'2 phase Delta loads not yet implemented.')
            else:
                bus_df.loc[actBus,'connect'] = 'Y'
                load_df.loc[actLoad,'connect'] = 'Y'
                if '1' in phs or phs==[]:
                    bus_df.loc[actBus,'Pa'] = LDS.kW()/nPh + bus_df.loc[actBus,'Pa']
                    bus_df.loc[actBus,'Qa'] = LDS.kvar()/nPh + bus_df.loc[actBus,'Qa']
                    load_df.loc[actLoad,'phase'] = [0]
                if '2' in phs or phs==[]:
                    bus_df.loc[actBus,'Pb'] = LDS.kW()/nPh + bus_df.loc[actBus,'Pb']
                    bus_df.loc[actBus,'Qb'] = LDS.kvar()/nPh + bus_df.loc[actBus,'Qb']
                    load_df.loc[actLoad,'phase'] = [1]
                if '3' in phs or phs==[]:
                    bus_df.loc[actBus,'Pc'] = LDS.kW()/nPh + bus_df.loc[actBus,'Pc']
                    bus_df.loc[actBus,'Qc'] = LDS.kvar()/nPh + bus_df.loc[actBus,'Qc']
                    load_df.loc[actLoad,'phase'] = [2]
                    if nPh==2:
                        print('Warning! Load: ',LDS.Name(),'2 phase Wye loads not yet implemented.')
                    if nPh==3:
                        load_df.at[actLoad,'phase'] = [0,1,2]				
            i = LDS.Next()

        # create bus_df from loads and capacitors ======

        if nCap>0:
            capColumns = ["name","number","bus","kVln","connect","Qa","Qb","Qc"]
            cap_df = pd.DataFrame(data=np.zeros((nCap,len(capColumns))), index=CAP.AllNames(), columns=capColumns)
            cap_df.loc[:,'name'] = CAP.AllNames()
            i = CAP.First()
            while i:
                capName = CAP.Name()
                cap_df.loc[capName,'number'] = i-1 # not too clear rn what 'number' does?
                
                SetACE('Capacitor.'+capName)
                nPh = ACE.NumPhases()
                phs = ACE.BusNames()[0].split('.')[1:]
                
                actBus = ACE.BusNames()[0].split('.')[0]
                cap_df.loc[capName,'bus'] = actBus
                if CAP.IsDelta():
                    cap_df.loc[capName,'connect'] = 'D'
                    cap_df.loc[capName,'kVln'] = CAP.kV()/np.sqrt(3)
                    if nPh==1:
                        if '1' in phs and '2' in phs:
                            cap_df.loc[capName,'Qa'] = CAP.kvar() + cap_df.loc[capName,'Qa']
                        if '2' in phs and '3' in phs:
                            cap_df.loc[capName,'Qb'] = CAP.kvar() + cap_df.loc[capName,'Qb']
                        if '3' in phs and '1' in phs:
                            cap_df.loc[capName,'Qc'] = CAP.kvar() + cap_df.loc[capName,'Qc']
                    if nPh==3:
                        cap_df.loc[capName,'Qa'] = CAP.kvar()/3 + cap_df.loc[capName,'Qa']
                        cap_df.loc[capName,'Qb'] = CAP.kvar()/3 + cap_df.loc[capName,'Qb']
                        cap_df.loc[capName,'Qc'] = CAP.kvar()/3 + cap_df.loc[capName,'Qc']
                    if nPh==2:
                        print('Warning! Cap: ',CAP.Name(),'2 phase Delta loads not yet implemented.')
                else:
                    cap_df.loc[capName,'connect'] = 'Y'
                    if nPh==3:
                        cap_df.loc[capName,'kVln'] = CAP.kV()/np.sqrt(3) # NB: kV depends on the connection type + no. phases.
                    elif nPh==1:
                        cap_df.loc[capName,'kVln'] = CAP.kV()
                    if '1' in phs or phs==[]:
                        cap_df.loc[capName,'Qa'] = CAP.kvar()/nPh + cap_df.loc[capName,'Qa']
                    if '2' in phs or phs==[]:
                        cap_df.loc[capName,'Qb'] = CAP.kvar()/nPh + cap_df.loc[capName,'Qb']
                    if '3' in phs or phs==[]:
                        cap_df.loc[capName,'Qc'] = CAP.kvar()/nPh + cap_df.loc[capName,'Qc']
                i = CAP.Next()
 
        #%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # create Y matrix  directly from OpenDSS ======
        '''
        Calculate the nodal admittance matrix YMatrix (without including the vource and loads)
            
        Returns
        -------
        numpy array
                It returns a complex array of Y admittance matrix      
        -------        
        What's the difference for directly using Ybus   = dss.Circuit.SystemY(), will the vsources and loads have impacts on Y?
        Main difference is weather or not the voltage sources and loads are included into the Y matrix formulation                   
        '''
        #- disconnect vsources and loads
        #DSSText.Command='vsource.source.enabled = no'
        dss.Command('vsource.source.enabled = no')
        #DSSText.Command='batchedit load..* enabled=no'
        dss.Command('BatchEdit Load..* enabled=no')
        #- extract YMatrix
        #DSSText.Command='Solve'
        dss.Solution.Solve()
        Y =csc_matrix(dss.YMatrix.getYsparse())#DSSCircuit.SystemY  # --- slow
        #nNodes = DSSCircuit.NumNodes
        #
        #Yres = np.reshape(Ybus, (nNodes, nNodes*2))  # --- slow
        #
        ##- populate complex polyphase nodal Y
        #Y = np.zeros((nNodes, nNodes), dtype=complex)
        #for i in range(nNodes):  # --- very very slow
        #    Y[:, i] = Yres[:, 2*i] + 1j*Yres[:, 2*i + 1]
        #- reconnect vsources and loads
        #DSSText.Command='vsource.source.enabled = yes'
        dss.Command('vsource.source.enabled = yes')
        #DSSText.Command='batchedit load..* enabled=yes'
        dss.Command('batchedit load.. enabled=yes')
        #- return to the previous solution   
        #DSSText.Command='Solve'
        dss.Solution.Solve()
        
        nNodes2 = dss.Circuit.NumNodes()

        # Create full Y based on the Y extracted from OpenDSS ======
        #~~~~~~~~~~~~~~~~~~Full Y info extraction_(inclusing all phases)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
        FullPhaseNames=[]
        for i in BusNameList:
            FullPhaseNames.append(i +'.1')
            FullPhaseNames.append(i +'.2')
            FullPhaseNames.append(i +'.3')
        
        N_buses = len(BusNameList)
        N_phases= len(FullPhaseNames)
        
        # OpenDSS Y matrix to Full Y matrix conversion
        phase_index_trans= -1.0 * np.ones(len(FullPhaseNames))
        
        for nodename in FullPhaseNames:    
            if nodename.upper() in yList:
                phase_index_trans[FullPhaseNames.index(nodename)] = yList.index(nodename.upper())
        
        phase_index_trans= phase_index_trans.astype(int)
        Ynet = np.zeros([(N_buses)*3,(N_buses)*3],dtype=np.complex_)
        
        for i in range(N_phases):
             for j in range(N_phases):
                 if phase_index_trans[i]== -1 or phase_index_trans[j]== -1:
                     continue
                 else:
                     Ynet[i,j] = Y[ phase_index_trans[i], phase_index_trans[j] ]
            
        
        if save_model:
            if not os.path.exists(dir0):
                os.makedirs(dir0)
            src_df.to_csv(sn0+"_src_df.csv")
            bus_df.to_csv(sn0+"_bus_df.csv")
            load_df.to_csv(sn0+"_load_df.csv")
            line_df.to_csv(sn0+"_line_df.csv")         
            np.savetxt(sn0+"_Y_dss.csv", Y.data, delimiter=",")
            np.savetxt(sn0+"_Y_dss_full.csv", Ynet.data, delimiter=",")            
            solution_df.to_csv(sn0+"_solution_df.csv")
            if nTrn > 0:
                trn_df.to_csv(sn0+"_trn_df.csv")
            if nCap > 0:
                cap_df.to_csv(sn0+"_cap_df.csv")
            print('\nFiles saved in:\n',dir0+'\\','\n')
        return src_df, bus_df, load_df, line_df, solution_df, Y
        #else:
        #    return src_df, bus_df, load_df, line_df, solution_df, Y


# feeder = "13BusOxEmf"
# dssName = "IEEE13Nodeckt_z_oxemf"

# feeder = "eulv"
# dssName = "Master_z_oxemf"

# # Test 123 bus system info extraction
# feeder = "123Bus_original"
# dssName = "IEEE123Master"

# feeder = "4Bus-YY-Bal-xfmr"
# dssName = "4Bus-YY-Bal"

# for iFdr in range(1,6):
    # for jFdr in range(11):
        # feeder = 'n'+str(iFdr)+'_'+str(jFdr)
        # dssName = "master_oxemf"





# Model conversion test for Xcel Feeder
# feeder = "BTER1342B_latest"  # OpenDSS file folder name
# dssName = "Master" #DSS master file name 

# # # Model conversion test for Xcel Feeder
# feeder = "BTER1349B"  # OpenDSS file folder name
# dssName = "Master" #DSS master file name 

# # # Model conversion test for Xcel Feeder
# feeder = "BTER1356B"  # OpenDSS file folder name
# dssName = "Master" #DSS master file name 

# # # Model conversion test for Xcel Feeder
#feeder = "MEAD2104"  # OpenDSS file folder name
#dssName = "Master_mingzhi" #DSS master file name 

# feeder = "MEAD2059"  # OpenDSS file folder name
# dssName = "Master" #DSS master file name 

# feeder = "MURP1210"  # OpenDSS file folder name MEAD2059
# dssName = "Master" #DSS master file name 



# RUN ME HERE:
#dssConvert(feeder,dssName)




