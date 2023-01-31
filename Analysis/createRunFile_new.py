import sys
import os
import csv
import string
import datetime

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['data', 'data_control', 'MC'], help="Specify if data or Monte Carlo")
parser.add_argument("anatype", type=str, choices=['tau3mu', 'control'], help="Specify analysis type")
parser.add_argument("--run", type=str, default='', choices=['2022B', '2022C_0', '2022C_1', '2022C_2', '2022C_3', '2022C_4', '2022C_5', '2022C_6', '2022C_7', '2022D-v1_0', '2022D-v1_1', '2022D-v1_2', '2022D-v1_3', '2022D-v1_4', '2022D-v1_5', '2022D-v1_6', '2022D-v1_7', '2022D-v2_0', '2022D-v2_1', '2022D-v2_2', '2022D-v2_3', '2022D-v2_4', '2022D-v2_5', '2022D-v2_6', '2022D-v2_7', '2022E_0', '2022E_1', '2022E_2', '2022E_3', '2022E_4', '2022E_5', '2022E_6', '2022E_7', '2022F_0', '2022F_1', '2022F_2', '2022F_3', '2022F_4', '2022F_5', '2022F_6', '2022F_7', '2022G_0', '2022G_1', '2022G_2', '2022G_3', '2022G_4', '2022G_5', '2022G_6', '2022G_7'], help="run in data")
# Optional Arguments
parser.add_argument("--outName", type=str, default="test", help="Specify name for output files")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
parser.add_argument("--MCprocess", type=str, default='', choices=['Ds', 'B0', 'Bp', 'DsPhiPi', 'DsPhiMuNu'], help="process in Monte Carlo")
args = parser.parse_args()

#prepare output filename  and option string
if args.dataset == 'data':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+'" "'+args.run+'"'
elif args.dataset == 'data_control':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_control","")+'" "'+args.run+'"'
elif args.dataset == 'MC':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+'" "'+args.MCprocess+'"'

#startTime = datetime.datetime.now().strftime("%Y%m%d_%H%M")

# Create target Directory if don't exist
if args.dataset == 'MC':
   output_name = args.MCprocess+"_"+args.anatype+"_"+args.outName
else: 
   output_name = args.run+"_"+args.anatype+"_"+args.outName

if not os.path.exists(output_name):
    os.mkdir(output_name)
    print('Directory '+output_name+' created\n')
else:    
    print('Directory '+output_name+' already exists\n')

if args.anatype == 'tau3mu':
   #### 2022
   if args.dataset == 'data' and args.run == '2022B':
      path = '' 
   if args.dataset == 'data' and args.run == '2022C_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsTau3mu_2022eraC_stream0_Mini_v6/221215_112428' 
   if args.dataset == 'data' and args.run == '2022C_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsTau3mu_2022eraC_stream1_Mini_v6/221215_112500'
   if args.dataset == 'data' and args.run == '2022C_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsTau3mu_2022eraC_stream2_Mini_v6/221215_112531' 
   if args.dataset == 'data' and args.run == '2022C_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsTau3mu_2022eraC_stream3_Mini_v6/221215_112601' 
   if args.dataset == 'data' and args.run == '2022C_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsTau3mu_2022eraC_stream4_Mini_v6/221215_112633' 
   if args.dataset == 'data' and args.run == '2022C_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsTau3mu_2022eraC_stream5_Mini_v6/221215_112702' 
   if args.dataset == 'data' and args.run == '2022C_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsTau3mu_2022eraC_stream6_Mini_v6/221215_112731' 
   if args.dataset == 'data' and args.run == '2022C_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsTau3mu_2022eraC_stream7_Mini_v6/221215_112800' 
   if args.dataset == 'data' and args.run == '2022D-v1_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraD-v1_stream0_Mini_v6/221215_125019'
   if args.dataset == 'data' and args.run == '2022D-v1_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraD-v1_stream1_Mini_v6/221215_151532'
   if args.dataset == 'data' and args.run == '2022D-v1_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraD-v1_stream2_Mini_v6/221215_151601'
   if args.dataset == 'data' and args.run == '2022D-v1_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraD-v1_stream3_Mini_v6/221215_151631'
   if args.dataset == 'data' and args.run == '2022D-v1_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraD-v1_stream4_Mini_v6/221215_151701'
   if args.dataset == 'data' and args.run == '2022D-v1_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraD-v1_stream5_Mini_v6/221215_151729'
   if args.dataset == 'data' and args.run == '2022D-v1_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraD-v1_stream6_Mini_v6/221215_151758'
   if args.dataset == 'data' and args.run == '2022D-v1_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraD-v1_stream7_Mini_v6/221215_151827'
   if args.dataset == 'data' and args.run == '2022D-v2_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraD-v2_stream0_Mini_v6/221215_145359'
   if args.dataset == 'data' and args.run == '2022D-v2_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraD-v2_stream1_Mini_v6/221215_170419'
   if args.dataset == 'data' and args.run == '2022D-v2_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraD-v2_stream2_Mini_v6/221215_170453'
   if args.dataset == 'data' and args.run == '2022D-v2_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraD-v2_stream3_Mini_v6/221215_170528'
   if args.dataset == 'data' and args.run == '2022D-v2_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraD-v2_stream4_Mini_v6/221215_170603'
   if args.dataset == 'data' and args.run == '2022D-v2_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraD-v2_stream5_Mini_v6/221215_170636'
   if args.dataset == 'data' and args.run == '2022D-v2_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraD-v2_stream6_Mini_v6/221215_170709'
   if args.dataset == 'data' and args.run == '2022D-v2_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraD-v2_stream7_Mini_v6/221215_170744'
   if args.dataset == 'data' and args.run == '2022E_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraE_stream0_Mini_v6/221215_173406'
   if args.dataset == 'data' and args.run == '2022E_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraE_stream1_Mini_v6/221215_180913'
   if args.dataset == 'data' and args.run == '2022E_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraE_stream2_Mini_v6/221215_180946'
   if args.dataset == 'data' and args.run == '2022E_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraE_stream3_Mini_v6/221215_181017'
   if args.dataset == 'data' and args.run == '2022E_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraE_stream4_Mini_v6/221215_181048'
   if args.dataset == 'data' and args.run == '2022E_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraE_stream5_Mini_v6/221215_181121'
   if args.dataset == 'data' and args.run == '2022E_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraE_stream6_Mini_v6/221215_181153'
   if args.dataset == 'data' and args.run == '2022E_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraE_stream7_Mini_v6/221215_181227'
   if args.dataset == 'data' and args.run == '2022F_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraF_stream0_Mini_v6/221216_081751'
   if args.dataset == 'data' and args.run == '2022F_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraF_stream1_Mini_v6/221216_083140'
   if args.dataset == 'data' and args.run == '2022F_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraF_stream2_Mini_v6/221216_083225'
   if args.dataset == 'data' and args.run == '2022F_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraF_stream3_Mini_v6/221216_083254'
   if args.dataset == 'data' and args.run == '2022F_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraF_stream4_Mini_v6/221216_083323'
   if args.dataset == 'data' and args.run == '2022F_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraF_stream5_Mini_v6/221216_083351'
   if args.dataset == 'data' and args.run == '2022F_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraF_stream6_Mini_v6/221216_083419'
   if args.dataset == 'data' and args.run == '2022F_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraF_stream7_Mini_v6/221216_083447'
   if args.dataset == 'data' and args.run == '2022G_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraG_stream0_Mini_v6/221216_082759'
   if args.dataset == 'data' and args.run == '2022G_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraG_stream1_Mini_v6/221216_082828'
   if args.dataset == 'data' and args.run == '2022G_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraG_stream2_Mini_v6/221216_082857'
   if args.dataset == 'data' and args.run == '2022G_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraG_stream3_Mini_v6/221219_090816'
   if args.dataset == 'data' and args.run == '2022G_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraG_stream4_Mini_v6/221216_082954'
   if args.dataset == 'data' and args.run == '2022G_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraG_stream5_Mini_v6/221219_091945'
   if args.dataset == 'data' and args.run == '2022G_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraG_stream6_Mini_v6/221219_092557'
   if args.dataset == 'data' and args.run == '2022G_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraG_stream7_Mini_v6/221219_093100'


if args.anatype == 'control':
   if args.dataset == 'data_control' and args.run == '2022B':
      path = ''
   if args.dataset == 'data_control' and args.run == '2022C_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsPhiPi_2022eraC_stream0_Mini_v6/221216_164245'
   if args.dataset == 'data_control' and args.run == '2022C_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsPhiPi_2022eraC_stream1_Mini_v6/221217_080207'
   if args.dataset == 'data_control' and args.run == '2022C_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsPhiPi_2022eraC_stream2_Mini_v6/221217_080239'
   if args.dataset == 'data_control' and args.run == '2022C_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsPhiPi_2022eraC_stream3_Mini_v6/221217_080311'
   if args.dataset == 'data_control' and args.run == '2022C_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsPhiPi_2022eraC_stream4_Mini_v6/221217_080343'
   if args.dataset == 'data_control' and args.run == '2022C_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsPhiPi_2022eraC_stream5_Mini_v6/221217_080416'
   if args.dataset == 'data_control' and args.run == '2022C_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsPhiPi_2022eraC_stream6_Mini_v6/221217_080448'
   if args.dataset == 'data_control' and args.run == '2022C_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsPhiPi_2022eraC_stream7_Mini_v6/221217_144720'
   if args.dataset == 'data_control' and args.run == '2022D-v1_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsPhiPi_2022eraD-v1_stream0_Mini_v6/221216_170014'
   if args.dataset == 'data_control' and args.run == '2022D-v1_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsPhiPi_2022eraD-v1_stream1_Mini_v6/221217_081546'
   if args.dataset == 'data_control' and args.run == '2022D-v1_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsPhiPi_2022eraD-v1_stream2_Mini_v6/221217_081616'
   if args.dataset == 'data_control' and args.run == '2022D-v1_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsPhiPi_2022eraD-v1_stream3_Mini_v6/221217_081648'
   if args.dataset == 'data_control' and args.run == '2022D-v1_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsPhiPi_2022eraD-v1_stream4_Mini_v6/221217_081719'
   if args.dataset == 'data_control' and args.run == '2022D-v1_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsPhiPi_2022eraD-v1_stream5_Mini_v6/221217_081750'
   if args.dataset == 'data_control' and args.run == '2022D-v1_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsPhiPi_2022eraD-v1_stream6_Mini_v6/221217_081822'
   if args.dataset == 'data_control' and args.run == '2022D-v1_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsPhiPi_2022eraD-v1_stream7_Mini_v6/221217_081854'
   if args.dataset == 'data_control' and args.run == '2022D-v2_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsPhiPi_2022eraD-v2_stream0_Mini_v6/221216_170108'
   if args.dataset == 'data_control' and args.run == '2022D-v2_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsPhiPi_2022eraD-v2_stream1_Mini_v6/221217_083112'
   if args.dataset == 'data_control' and args.run == '2022D-v2_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsPhiPi_2022eraD-v2_stream2_Mini_v6/221217_083145'
   if args.dataset == 'data_control' and args.run == '2022D-v2_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsPhiPi_2022eraD-v2_stream3_Mini_v6/221217_083217'
   if args.dataset == 'data_control' and args.run == '2022D-v2_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsPhiPi_2022eraD-v2_stream4_Mini_v6/221217_083249'
   if args.dataset == 'data_control' and args.run == '2022D-v2_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsPhiPi_2022eraD-v2_stream5_Mini_v6/221217_083321'
   if args.dataset == 'data_control' and args.run == '2022D-v2_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsPhiPi_2022eraD-v2_stream6_Mini_v6/221217_083352'
   if args.dataset == 'data_control' and args.run == '2022D-v2_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsPhiPi_2022eraD-v2_stream7_Mini_v6/221217_083425'
   if args.dataset == 'data_control' and args.run == '2022E_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsPhiPi_2022eraE_stream0_Mini_v6/221216_170414'
   if args.dataset == 'data_control' and args.run == '2022E_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsPhiPi_2022eraE_stream1_Mini_v6/221217_084324'
   if args.dataset == 'data_control' and args.run == '2022E_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsPhiPi_2022eraE_stream2_Mini_v6/221217_084356'
   if args.dataset == 'data_control' and args.run == '2022E_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsPhiPi_2022eraE_stream3_Mini_v6/221217_084428'
   if args.dataset == 'data_control' and args.run == '2022E_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsPhiPi_2022eraE_stream4_Mini_v6/221217_084459'
   if args.dataset == 'data_control' and args.run == '2022E_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsPhiPi_2022eraE_stream5_Mini_v6/221217_084529'
   if args.dataset == 'data_control' and args.run == '2022E_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsPhiPi_2022eraE_stream6_Mini_v6/221217_084558'
   if args.dataset == 'data_control' and args.run == '2022E_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsPhiPi_2022eraE_stream7_Mini_v6/221217_084628'
   if args.dataset == 'data_control' and args.run == '2022F_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsPhiPi_2022eraF_stream0_Mini_v6/221217_145819'
   if args.dataset == 'data_control' and args.run == '2022F_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsPhiPi_2022eraF_stream1_Mini_v6/221217_085156'
   if args.dataset == 'data_control' and args.run == '2022F_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsPhiPi_2022eraF_stream2_Mini_v6/221217_145850'
   if args.dataset == 'data_control' and args.run == '2022F_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsPhiPi_2022eraF_stream3_Mini_v6/221217_145919'
   if args.dataset == 'data_control' and args.run == '2022F_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsPhiPi_2022eraF_stream4_Mini_v6/221217_145948'
   if args.dataset == 'data_control' and args.run == '2022F_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsPhiPi_2022eraF_stream5_Mini_v6/221217_150018'
   if args.dataset == 'data_control' and args.run == '2022F_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsPhiPi_2022eraF_stream6_Mini_v6/221217_150047'
   if args.dataset == 'data_control' and args.run == '2022F_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsPhiPi_2022eraF_stream7_Mini_v6/221217_150115'
   if args.dataset == 'data_control' and args.run == '2022G_0':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass0/SkimDsPhiPi_2022eraG_stream0_Mini_v6/221216_170729'
   if args.dataset == 'data_control' and args.run == '2022G_1':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass1/SkimDsPhiPi_2022eraG_stream1_Mini_v6/221217_085652'
   if args.dataset == 'data_control' and args.run == '2022G_2':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass2/SkimDsPhiPi_2022eraG_stream2_Mini_v6/221217_085724'
   if args.dataset == 'data_control' and args.run == '2022G_3':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass3/SkimDsPhiPi_2022eraG_stream3_Mini_v6/221217_085756'
   if args.dataset == 'data_control' and args.run == '2022G_4':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass4/SkimDsPhiPi_2022eraG_stream4_Mini_v6/221217_085826'
   if args.dataset == 'data_control' and args.run == '2022G_5':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass5/SkimDsPhiPi_2022eraG_stream5_Mini_v6/221217_085855'
   if args.dataset == 'data_control' and args.run == '2022G_6':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass6/SkimDsPhiPi_2022eraG_stream6_Mini_v6/221217_085928'
   if args.dataset == 'data_control' and args.run == '2022G_7':
      path = '/lustre/cms/store/user/caruta/ParkingDoubleMuonLowMass7/SkimDsPhiPi_2022eraG_stream7_Mini_v6/221217_085959'

if args.dataset == 'MC' and args.MCprocess == 'Ds':
   path = '/lustre/cms/store/user/caruta/Pythia8_DsTau3mu_Run3_2022/SkimTau3mu_MCRun3_Ds_Mini_v4/221213_163730'
if args.dataset == 'MC' and args.MCprocess == 'Bp':
   path = '/lustre/cms/store/user/caruta/Pythia8_BuTau3mu_Run3_2022/SkimTau3mu_MCRun3_Bu_Mini_v4/221214_073101'
if args.dataset == 'MC' and args.MCprocess == 'B0':
   path = '/lustre/cms/store/user/caruta/Pythia8_BdTau3mu_Run3_2022/SkimTau3mu_MCRun3_Bd_Mini_v4/221214_073125'
if args.dataset == 'MC' and args.MCprocess == 'DsPhiMuNu':
   path = '/lustre/cms/store/user/caruta/Pythia8_DsPhiMuNu_Run3_2022/SkimTau3mu_MCRun3_DsPhiMuNu_Mini_v4/230114_172745'
if args.dataset == 'MC' and args.MCprocess == 'DsPhiPi':
   path = '/lustre/cms/store/user/caruta/Pythia8_DsPhiPi_MuMuPi_Run3_2022/SkimPhiPi_MCRun3_Mini_v4/221213_155107'


#generating the list of all .root files in given directory and subdirectories
fileList = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            fileList.append(os.path.join(r, file))
            print file

#prepare final script
#final_script = open("submit_analysis_"+startTime+".sh", "w")
final_script = open("submit_analysis_"+output_name+".sh", "w")
final_script.write("#!/bin/bash\n")
final_script.write("chmod 777 -R *\n")
final_script.write("cd "+output_name+"\n")

#loop to generate one .cpp+executable+batch system conf file for each group of "n" files
n_chunk = len(fileList)//args.n
print('Number of files is {0:2d}'.format(len(fileList)))
print('Number of jobs is {0:2d}'.format(n_chunk+1))
for file_index in range(n_chunk+1):
      chunk = '' 
      for idx, l in enumerate(fileList):
         if idx < args.n*(file_index+1) and idx >= args.n*file_index:
             l = l.rstrip()
             l = '        chain->AddFile("{}");\n'.format(l)
             chunk = chunk + l

      #analysis.cpp template
      with open("templates/Analysis_template.cpp", "r") as in_file:
          buf = in_file.readlines()

      cpp_filename = "Analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_chunk"+str(file_index)+".cpp"
      with open(cpp_filename, "w") as out_file:
          for lb in buf:
              if lb == '        //AddFile_'+args.dataset+args.MCprocess+'_'+args.anatype+'\n':
                  #write group of files
                  out_file.write(chunk)
              elif lb == '        //OutFile_'+args.dataset+args.MCprocess+'_'+args.anatype+'\n':
                  #write output file name
                  out_file.write('        fileout = "'+out_filename+str(file_index)+'.root";\n')
              else: out_file.write(lb)

              #elif lb == '            TString fileout = "AddOutput_'+args.dataset+args.MCprocess+'_'+args.anatype+'.root";\n':
                  #write output file name
               #   out_file.write('        TString fileout = "'+out_filename+str(file_index)+'.root";\n')
              #else: out_file.write(lb)

      #executable template
      with open("templates/launch_analysis_template.job", "r") as launch_infile:
          buf2 = launch_infile.readlines()

      launch_filename = "launch_analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_name+"/"+launch_filename, "w") as launch_outfile:
          for lb2 in buf2:
              if lb2 == "#compile\n":
                  launch_outfile.write("cd "+output_name+"\n")
                  launch_outfile.write("g++ -I $ROOTSYS/include ../"+cpp_filename+" `root-config --glibs` `root-config --libs` `root-config --cflags` -lTMVA -L $ROOTSYS/lib -o executable"+str(file_index)+"\n")
              elif lb2 == "#execute\n":
                  launch_outfile.write('./executable'+str(file_index)+option_string+'\n')
              else: launch_outfile.write(lb2)

      #myCondor template
      with open("templates/my_HTCondor_template.job", "r") as myCondor_infile:
          buf3 = myCondor_infile.readlines()

      condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_name+"/"+condor_filename, "w") as myCondor_outfile:
          for lb3 in buf3:
              if lb3 == "Executable = launch_analysis_template.job\n":
                  myCondor_outfile.write("Executable = "+launch_filename+"\n")
              else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo condor_submit "+condor_filename+" -name ettore\n")
      final_script.write("condor_submit "+condor_filename+" -name ettore\n")

final_script.close()
#submitName = "submit_analysis_"+startTime+".sh"
#source submitName
