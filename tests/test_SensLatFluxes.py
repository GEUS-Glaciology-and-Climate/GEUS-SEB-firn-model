import unittest
import __init__
import pandas as pd
import json
import numpy as np

from main_SEB_firn import setConstants
from lib_seb_smb_model import SensLatFluxes_bulk_old, SensLatFluxes_bulk_opt, HHsubsurf, SmoothSurf_opt
import lib_io as io

class SensLatFluxesTestCase(unittest.TestCase):
    '''Test the SensLatFluxes_bulk function'''

    def testMultipleTimes(self):
        '''Run syntethic simulation multiple times, compare result'''

        # Initializing parameters with static, synthetic values 
        # From run k = 4, and k =5999
        WS = [4.88, 4.96]
        nu = [1.5471445341359545e-05, 1.5097488894719079e-05]
        q = [0.002983881449656623, 0.0022233120103667815]
        snowthick = [1.0, 1.0]
        Tsurf = [50.00032543036966, 50.18252389141013]
        theta = [269.23442297512435, 267.1981547860696]
        theta_v = [269.72264072580384, 267.5591781294883]
        pres = [857.105524, 866.39]
        rho_atm = [1.109137923478789, 1.1296708808544773]
        z_WS = [2.9995, 2.358]
        z_T = [2.4995, 1.858]
        z_RH = [2.4995, 1.858]
        z_0 = [0.0013, 0.0013]
        c = setConstants()
       
        L = [0,0]
        LHF = [0,0]
        SHF = [0,0]
        theta_2m = [0,0]
        q_2m = [0,0]
        ws_10m = [0,0]
        Re = [0,0]

        L_opt = [0,0]  
        LHF_opt = [0,0]
        SHF_opt = [0,0]
        theta_2m_opt = [0,0] 
        q_2m_opt = [0,0]
        ws_10m_opt = [0,0] 
        Re_opt = [0,0]
         
        for k in range(2):
            for i in range(1000):

                (
                    L[k], 
                    LHF[k],
                    SHF[k], 
                    theta_2m[k], 
                    q_2m[k], 
                    ws_10m[k], 
                    Re[k]
                ) = SensLatFluxes_bulk_old(
                    WS[k],
                    nu[k],
                    q[k],
                    snowthick[k],
                    Tsurf[k],
                    theta[k],
                    theta_v[k],
                    pres[k],
                    rho_atm[k],
                    z_WS[k],
                    z_T[k],
                    z_RH[k],
                    z_0[k],
                    c
                )

                (
                    L_opt[k], 
                    LHF_opt[k],
                    SHF_opt[k], 
                    theta_2m_opt[k], 
                    q_2m_opt[k], 
                    ws_10m_opt[k], 
                    Re_opt[k]
                ) = SensLatFluxes_bulk_opt(
                    WS[k],
                    nu[k],
                    q[k],
                    snowthick[k],
                    Tsurf[k],
                    theta[k],
                    theta_v[k],
                    pres[k],
                    rho_atm[k],
                    z_WS[k],
                    z_T[k],
                    z_RH[k],
                    z_0[k],
                    c
                ) 

            # Assert that the values from original and optimized code are the same
            assert L == L_opt
            assert LHF == LHF_opt
            assert SHF == SHF_opt
            assert theta_2m == theta_2m_opt
            assert q_2m == q_2m_opt
            assert ws_10m == ws_10m_opt
            assert Re == Re_opt

    # def testSmoothSurf():

    #     # Initializing parameters with static, synthetic values from run k=4
    #     WS= 4.88
    #     nu= 1.5471445341359545e-05
    #     q = 0.002983881449656623
    #     snowthick = 1.0
    #     Tsurf = 50.00032543036966
    #     theta = 269.23442297512435
    #     theta_v = 269.72264072580384
    #     pres = 857.105524
    #     rho_atm = 1.109137923478789
    #     z_WS = 2.9995
    #     z_T = 2.4995
    #     z_RH = 2.4995
    #     z_0 = 0.0013
    #     c = setConstants()

        #res = SmoothSurf_opt(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)
        #res2 = SmoothSurf_opt_2(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)


def writeArray():
    arr = [31.15625, 31.5625, 31.53125, 31.328125, 31.609375, 31.5, 31.40625, 31.484375, 31.34375, 31.359375]
    for i in range(len(arr)):
        print(arr[i])

def writeOutput():
        '''Test and compare output from SensLatFluxes with old code output'''
        c = setConstants()

        with open("parameters.json") as parameter_file:
            parameters = json.load(parameter_file)

        weather_data_input_path = str(parameters['weather_data']['weather_input_path'])
        df_aws = io.load_promice(weather_data_input_path)[:6000]
        df_aws = df_aws.set_index("time").resample("H").mean()

        # Call HHsubsurf, which gets L, LHF, SHF, theta_2m, q_2m, ws_10m and Re from SensLatFluxes
        (L,
        LHF,
        SHF,
        theta_2m,
        q_2m,
        ws_10m,
        Re,
        melt_mweq,
        sublimation_mweq,
        snowc,
        snic,
        slwc,
        T_ice,
        zrfrz,
        rhofirn,
        zsupimp,
        dgrain,
        zrogl,
        Tsurf,
        grndc,
        grndd,
        pgrndcapc,
        pgrndhflx,
        dH_comp,
        snowbkt,
        compaction) = HHsubsurf(df_aws, c)
        
        #Read in the above created files, compare with new calculations
        path = r'C:\Users\brink\Documents\Exjobb\Old_version_code\GEUS-SEB-firn-model\tests\ValuesSensLatFlux'
        old_L = pd.read_csv(path + '\L_OldCodeSensLatFlux.csv', index_col=False)
        new_L = pd.DataFrame(L)
        old_LHF = pd.read_csv(path + '\LHF_OldCodeSensLatFlux.csv', index_col=False)
        new_LHF = pd.DataFrame(LHF)
        old_SHF = pd.read_csv(path + '\SHF_OldCodeSensLatFlux.csv', index_col=False)
        new_SHF = pd.DataFrame(SHF)
        old_theta_2m = pd.read_csv(r'C:\Users\brink\Documents\Exjobb\Old_version_code\GEUS-SEB-firn-model\tests\ValuesSensLatFlux\theta_2m_OldCodeSensLatFlux.csv', index_col=False)
        new_theta_2m = pd.DataFrame(theta_2m)
        old_q_2m = pd.read_csv(path + '\q_2m_OldCodeSensLatFlux.csv', index_col=False)
        new_q_2m = pd.DataFrame(q_2m)
        old_ws_10m = pd.read_csv(path + '\ws_10m_OldCodeSensLatFlux.csv', index_col=False)
        new_ws_10m = pd.DataFrame(ws_10m)
        old_Re = pd.read_csv(path + '\Re_OldCodeSensLatFlux.csv', index_col=False)
        new_Re = pd.DataFrame(Re)
        
        # Write results from SensLatFluxes_bulk to one csv file
        # new_data = new_L
        # new_data = new_data.assign(new_LHF = new_LHF)
        # new_data = new_data.assign(new_SHF = new_SHF)
        # new_data = new_data.assign(new_theta_2m = new_theta_2m)
        # new_data = new_data.assign(new_q_2m = new_q_2m)
        # new_data = new_data.assign(new_ws_10m = new_ws_10m)
        # new_data = new_data.assign(new_Re = new_Re)

        # filename = 'res_SensLatFlux.csv'
        # pd.DataFrame(new_data).to_csv(filename)

        # print("End of comparision.")


        #Print differences between old output to newly computed output
        print("Old L:")
        print(old_L)
        print("New L:")
        print(new_L)

        print("Old LHF: ")
        print(old_LHF)
        print("New LHF: ")
        print(new_LHF)

        print("Old SHF: ")
        print(old_SHF)
        print("New SHF: ")
        print(new_SHF)
        
        print("Old theta_2m: ")
        print(old_theta_2m)
        print("New theta_2m: ")
        print(new_theta_2m)

        print("Old q_2m: ")
        print(old_q_2m)
        print("New q_2m: ")
        print(new_q_2m)

        print("Old ws_10m: ")
        print(old_ws_10m)
        print("New ws_10m: ")
        print(new_ws_10m) 

        print("Old Re: ")
        print(old_Re)
        print("New Re: ")
        print(new_Re)


#writeOutput()

#writeArray()

if __name__ == '__main__':
    unittest.main()