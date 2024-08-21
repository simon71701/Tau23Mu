import sys
import os
import csv
import string
import datetime

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['data', 'data_control','data_phimunu', 'MC'], help="Specify if data or Monte Carlo")
parser.add_argument("anatype", type=str, choices=['tau3mu', 'control','phimunu'], help="Specify analysis type")
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
elif args.dataset == 'data_phimunu':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_phimunu","")+'" "'+args.run+'"'
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
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsTau3mu_2022eraC_stream0_Mini_v3/240725_160353' 
   if args.dataset == 'data' and args.run == '2022C_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsTau3mu_2022eraC_stream1_Mini_v3/240725_160422'
   if args.dataset == 'data' and args.run == '2022C_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsTau3mu_2022eraC_stream2_Mini_v3/240725_160452' 
   if args.dataset == 'data' and args.run == '2022C_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsTau3mu_2022eraC_stream3_Mini_v3/240725_160521' 
   if args.dataset == 'data' and args.run == '2022C_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsTau3mu_2022eraC_stream4_Mini_v3/240725_160550' 
   if args.dataset == 'data' and args.run == '2022C_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsTau3mu_2022eraC_stream5_Mini_v3/240725_160620' 
   if args.dataset == 'data' and args.run == '2022C_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsTau3mu_2022eraC_stream6_Mini_v3/240725_160649' 
   if args.dataset == 'data' and args.run == '2022C_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsTau3mu_2022eraC_stream7_Mini_v3/240725_160718' 
   if args.dataset == 'data' and args.run == '2022D-v1_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraD-v1_stream0_Mini_v3/240725_175346'
   if args.dataset == 'data' and args.run == '2022D-v1_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraD-v1_stream1_Mini_v3/240725_175416'
   if args.dataset == 'data' and args.run == '2022D-v1_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraD-v1_stream2_Mini_v3/240725_175446'
   if args.dataset == 'data' and args.run == '2022D-v1_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraD-v1_stream3_Mini_v3/240725_175515'
   if args.dataset == 'data' and args.run == '2022D-v1_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraD-v1_stream4_Mini_v3/240725_175545'
   if args.dataset == 'data' and args.run == '2022D-v1_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraD-v1_stream5_Mini_v3/240725_175614'
   if args.dataset == 'data' and args.run == '2022D-v1_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraD-v1_stream6_Mini_v3/240725_175644'
   if args.dataset == 'data' and args.run == '2022D-v1_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraD-v1_stream7_Mini_v3/240725_175714'
   if args.dataset == 'data' and args.run == '2022D-v2_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraD-v2_stream0_Mini_v3/240725_174955'
   if args.dataset == 'data' and args.run == '2022D-v2_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraD-v2_stream1_Mini_v3/240725_175025'
   if args.dataset == 'data' and args.run == '2022D-v2_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraD-v2_stream2_Mini_v3/240725_175054'
   if args.dataset == 'data' and args.run == '2022D-v2_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraD-v2_stream3_Mini_v3/240725_175123'
   if args.dataset == 'data' and args.run == '2022D-v2_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraD-v2_stream4_Mini_v3/240725_175156'
   if args.dataset == 'data' and args.run == '2022D-v2_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraD-v2_stream5_Mini_v3/240725_175226'
   if args.dataset == 'data' and args.run == '2022D-v2_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraD-v2_stream6_Mini_v3/240725_175255'
   if args.dataset == 'data' and args.run == '2022D-v2_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraD-v2_stream7_Mini_v3/240725_175324'
   if args.dataset == 'data' and args.run == '2022E_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraE_stream0_Mini_v3/240725_174548'
   if args.dataset == 'data' and args.run == '2022E_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraE_stream1_Mini_v3/240725_174618'
   if args.dataset == 'data' and args.run == '2022E_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraE_stream2_Mini_v3/240725_174647'
   if args.dataset == 'data' and args.run == '2022E_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraE_stream3_Mini_v3/240725_174717'
   if args.dataset == 'data' and args.run == '2022E_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraE_stream4_Mini_v3/240725_174746'
   if args.dataset == 'data' and args.run == '2022E_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraE_stream5_Mini_v3/240725_174816'
   if args.dataset == 'data' and args.run == '2022E_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraE_stream6_Mini_v3/240725_174847'
   if args.dataset == 'data' and args.run == '2022E_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraE_stream7_Mini_v3/240725_174918'
   if args.dataset == 'data' and args.run == '2022F_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraF_stream0_Mini_v3/240725_174201'
   if args.dataset == 'data' and args.run == '2022F_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraF_stream1_Mini_v3/240725_174231'
   if args.dataset == 'data' and args.run == '2022F_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraF_stream2_Mini_v3/240725_174301'
   if args.dataset == 'data' and args.run == '2022F_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraF_stream3_Mini_v3/240725_174331'
   if args.dataset == 'data' and args.run == '2022F_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraF_stream4_Mini_v3/240725_174400'
   if args.dataset == 'data' and args.run == '2022F_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraF_stream5_Mini_v3/240725_174429'
   if args.dataset == 'data' and args.run == '2022F_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraF_stream6_Mini_v3/240725_174458'
   if args.dataset == 'data' and args.run == '2022F_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraF_stream7_Mini_v3/240725_174528'
   if args.dataset == 'data' and args.run == '2022G_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimTau3mu_2022eraG_stream0_Mini_v3/240725_173710'
   if args.dataset == 'data' and args.run == '2022G_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimTau3mu_2022eraG_stream1_Mini_v3/240725_173739'
   if args.dataset == 'data' and args.run == '2022G_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimTau3mu_2022eraG_stream2_Mini_v3/240725_173809'
   if args.dataset == 'data' and args.run == '2022G_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimTau3mu_2022eraG_stream3_Mini_v3/240725_173838'
   if args.dataset == 'data' and args.run == '2022G_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimTau3mu_2022eraG_stream4_Mini_v3/240725_173907'
   if args.dataset == 'data' and args.run == '2022G_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimTau3mu_2022eraG_stream5_Mini_v3/240725_173937'
   if args.dataset == 'data' and args.run == '2022G_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimTau3mu_2022eraG_stream6_Mini_v3/240725_174006'
   if args.dataset == 'data' and args.run == '2022G_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimTau3mu_2022eraG_stream7_Mini_v3/240725_174036'


if args.anatype == 'control':
   if args.dataset == 'data_control' and args.run == '2022B':
      path = ''
   if args.dataset == 'data_control' and args.run == '2022C_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraC_stream0_Mini_v3/221216_164245'
   if args.dataset == 'data_control' and args.run == '2022C_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraC_stream1_Mini_v3/221217_080207'
   if args.dataset == 'data_control' and args.run == '2022C_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraC_stream2_Mini_v3/221217_080239'
   if args.dataset == 'data_control' and args.run == '2022C_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraC_stream3_Mini_v3/221217_080311'
   if args.dataset == 'data_control' and args.run == '2022C_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraC_stream4_Mini_v3/221217_080343'
   if args.dataset == 'data_control' and args.run == '2022C_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraC_stream5_Mini_v3/221217_080416'
   if args.dataset == 'data_control' and args.run == '2022C_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraC_stream6_Mini_v3/221217_080448'
   if args.dataset == 'data_control' and args.run == '2022C_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraC_stream7_Mini_v3/221217_144720'
   if args.dataset == 'data_control' and args.run == '2022D-v1_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraD-v1_stream0_Mini_v3/221216_170014'
   if args.dataset == 'data_control' and args.run == '2022D-v1_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraD-v1_stream1_Mini_v3/221217_081546'
   if args.dataset == 'data_control' and args.run == '2022D-v1_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraD-v1_stream2_Mini_v3/221217_081616'
   if args.dataset == 'data_control' and args.run == '2022D-v1_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraD-v1_stream3_Mini_v3/221217_081648'
   if args.dataset == 'data_control' and args.run == '2022D-v1_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraD-v1_stream4_Mini_v3/221217_081719'
   if args.dataset == 'data_control' and args.run == '2022D-v1_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraD-v1_stream5_Mini_v3/221217_081750'
   if args.dataset == 'data_control' and args.run == '2022D-v1_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraD-v1_stream6_Mini_v3/221217_081822'
   if args.dataset == 'data_control' and args.run == '2022D-v1_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraD-v1_stream7_Mini_v3/221217_081854'
   if args.dataset == 'data_control' and args.run == '2022D-v2_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraD-v2_stream0_Mini_v3/221216_170108'
   if args.dataset == 'data_control' and args.run == '2022D-v2_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraD-v2_stream1_Mini_v3/221217_083112'
   if args.dataset == 'data_control' and args.run == '2022D-v2_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraD-v2_stream2_Mini_v3/221217_083145'
   if args.dataset == 'data_control' and args.run == '2022D-v2_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraD-v2_stream3_Mini_v3/221217_083217'
   if args.dataset == 'data_control' and args.run == '2022D-v2_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraD-v2_stream4_Mini_v3/221217_083249'
   if args.dataset == 'data_control' and args.run == '2022D-v2_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraD-v2_stream5_Mini_v3/221217_083321'
   if args.dataset == 'data_control' and args.run == '2022D-v2_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraD-v2_stream6_Mini_v3/221217_083352'
   if args.dataset == 'data_control' and args.run == '2022D-v2_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraD-v2_stream7_Mini_v3/221217_083425'
   if args.dataset == 'data_control' and args.run == '2022E_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraE_stream0_Mini_v3/221216_170414'
   if args.dataset == 'data_control' and args.run == '2022E_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraE_stream1_Mini_v3/221217_084324'
   if args.dataset == 'data_control' and args.run == '2022E_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraE_stream2_Mini_v3/221217_084356'
   if args.dataset == 'data_control' and args.run == '2022E_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraE_stream3_Mini_v3/221217_084428'
   if args.dataset == 'data_control' and args.run == '2022E_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraE_stream4_Mini_v3/221217_084459'
   if args.dataset == 'data_control' and args.run == '2022E_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraE_stream5_Mini_v3/221217_084529'
   if args.dataset == 'data_control' and args.run == '2022E_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraE_stream6_Mini_v3/221217_084558'
   if args.dataset == 'data_control' and args.run == '2022E_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraE_stream7_Mini_v3/221217_084628'
   if args.dataset == 'data_control' and args.run == '2022F_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraF_stream0_Mini_v3/221217_145819'
   if args.dataset == 'data_control' and args.run == '2022F_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraF_stream1_Mini_v3/221217_085156'
   if args.dataset == 'data_control' and args.run == '2022F_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraF_stream2_Mini_v3/221217_145850'
   if args.dataset == 'data_control' and args.run == '2022F_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraF_stream3_Mini_v3/221217_145919'
   if args.dataset == 'data_control' and args.run == '2022F_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraF_stream4_Mini_v3/221217_145948'
   if args.dataset == 'data_control' and args.run == '2022F_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraF_stream5_Mini_v3/221217_150018'
   if args.dataset == 'data_control' and args.run == '2022F_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraF_stream6_Mini_v3/221217_150047'
   if args.dataset == 'data_control' and args.run == '2022F_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraF_stream7_Mini_v3/221217_150115'
   if args.dataset == 'data_control' and args.run == '2022G_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraG_stream0_Mini_v3/221216_170729'
   if args.dataset == 'data_control' and args.run == '2022G_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraG_stream1_Mini_v3/221217_085652'
   if args.dataset == 'data_control' and args.run == '2022G_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraG_stream2_Mini_v3/221217_085724'
   if args.dataset == 'data_control' and args.run == '2022G_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraG_stream3_Mini_v3/221217_085756'
   if args.dataset == 'data_control' and args.run == '2022G_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraG_stream4_Mini_v3/221217_085826'
   if args.dataset == 'data_control' and args.run == '2022G_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraG_stream5_Mini_v3/221217_085855'
   if args.dataset == 'data_control' and args.run == '2022G_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraG_stream6_Mini_v3/221217_085928'
   if args.dataset == 'data_control' and args.run == '2022G_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraG_stream7_Mini_v3/221217_085959'

if args.anatype == 'phimunu':
   if args.dataset == 'data_phimunu' and args.run == '2022B':
      path = ''
   if args.dataset == 'data_phimunu' and args.run == '2022C_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraC_stream0_Mini_v3/240725_200149/'
   if args.dataset == 'data_phimunu' and args.run == '2022C_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraC_stream1_Mini_v3/240725_200219'
   if args.dataset == 'data_phimunu' and args.run == '2022C_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraC_stream2_Mini_v3/240725_200253'
   if args.dataset == 'data_phimunu' and args.run == '2022C_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraC_stream3_Mini_v3/240725_200324'
   if args.dataset == 'data_phimunu' and args.run == '2022C_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraC_stream4_Mini_v3/240725_200354'
   if args.dataset == 'data_phimunu' and args.run == '2022C_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraC_stream5_Mini_v3/240725_200424'
   if args.dataset == 'data_phimunu' and args.run == '2022C_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraC_stream6_Mini_v3/240725_200454'
   if args.dataset == 'data_phimunu' and args.run == '2022C_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraC_stream7_Mini_v3/240725_200524'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraD-v1_stream0_Mini_v3/240725_200903'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraD-v1_stream1_Mini_v3/240725_200933'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraD-v1_stream2_Mini_v3/240725_201006'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraD-v1_stream3_Mini_v3/240725_201036'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraD-v1_stream4_Mini_v3/240725_201108'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraD-v1_stream5_Mini_v3/240725_201137'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraD-v1_stream6_Mini_v3/240725_201208'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v1_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraD-v1_stream7_Mini_v3/240725_201237'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraD-v2_stream0_Mini_v3/240725_201407'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraD-v2_stream1_Mini_v3/240725_201437'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraD-v2_stream2_Mini_v3/240725_201508'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraD-v2_stream3_Mini_v3/240725_201538'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraD-v2_stream4_Mini_v3/240725_201609'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraD-v2_stream5_Mini_v3/240725_201640'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraD-v2_stream6_Mini_v3/240725_201710'
   if args.dataset == 'data_phimunu' and args.run == '2022D-v2_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraD-v2_stream7_Mini_v3/240725_201739'
   if args.dataset == 'data_phimunu' and args.run == '2022E_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraE_stream0_Mini_v3/240725_202012'
   if args.dataset == 'data_phimunu' and args.run == '2022E_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraE_stream1_Mini_v3/240725_202043'
   if args.dataset == 'data_phimunu' and args.run == '2022E_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraE_stream2_Mini_v3/240725_202115'
   if args.dataset == 'data_phimunu' and args.run == '2022E_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraE_stream3_Mini_v3/240725_202146'
   if args.dataset == 'data_phimunu' and args.run == '2022E_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraE_stream4_Mini_v3/240725_202216'
   if args.dataset == 'data_phimunu' and args.run == '2022E_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraE_stream5_Mini_v3/240725_202245'
   if args.dataset == 'data_phimunu' and args.run == '2022E_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraE_stream6_Mini_v3/240725_202316'
   if args.dataset == 'data_phimunu' and args.run == '2022E_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraE_stream7_Mini_v3/240725_202345'
   if args.dataset == 'data_phimunu' and args.run == '2022F_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraF_stream0_Mini_v3/240725_202614'
   if args.dataset == 'data_phimunu' and args.run == '2022F_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraF_stream1_Mini_v3/240725_202645'
   if args.dataset == 'data_phimunu' and args.run == '2022F_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraF_stream2_Mini_v3/240725_202716'
   if args.dataset == 'data_phimunu' and args.run == '2022F_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraF_stream3_Mini_v3/240725_202746'
   if args.dataset == 'data_phimunu' and args.run == '2022F_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraF_stream4_Mini_v3/240725_202817'
   if args.dataset == 'data_phimunu' and args.run == '2022F_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraF_stream5_Mini_v3/240725_202846'
   if args.dataset == 'data_phimunu' and args.run == '2022F_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraF_stream6_Mini_v3/240725_202915'
   if args.dataset == 'data_phimunu' and args.run == '2022F_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraF_stream7_Mini_v3/240725_202945'
   if args.dataset == 'data_phimunu' and args.run == '2022G_0':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass0/SkimDsPhiMuNu_2022eraG_stream0_Mini_v3/240725_203125'
   if args.dataset == 'data_phimunu' and args.run == '2022G_1':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass1/SkimDsPhiMuNu_2022eraG_stream1_Mini_v3/240725_203155'
   if args.dataset == 'data_phimunu' and args.run == '2022G_2':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass2/SkimDsPhiMuNu_2022eraG_stream2_Mini_v3/240725_203225'
   if args.dataset == 'data_phimunu' and args.run == '2022G_3':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass3/SkimDsPhiMuNu_2022eraG_stream3_Mini_v3/240725_203256'
   if args.dataset == 'data_phimunu' and args.run == '2022G_4':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass4/SkimDsPhiMuNu_2022eraG_stream4_Mini_v3/240725_203326'
   if args.dataset == 'data_phimunu' and args.run == '2022G_5':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass5/SkimDsPhiMuNu_2022eraG_stream5_Mini_v3/240725_203358'
   if args.dataset == 'data_phimunu' and args.run == '2022G_6':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass6/SkimDsPhiMuNu_2022eraG_stream6_Mini_v3/240725_203429'
   if args.dataset == 'data_phimunu' and args.run == '2022G_7':
      path = '/store/user/jschulte/ParkingDoubleMuonLowMass7/SkimDsPhiMuNu_2022eraG_stream7_Mini_v3/240725_203459'




if args.dataset == 'MC' and args.MCprocess == 'Ds':
   path = '/store/user/jschulte/Pythia8_DsTau3mu_Run3_2022/SkimTau3mu_MCRun3_Ds_Mini_v4/221213_163730'
if args.dataset == 'MC' and args.MCprocess == 'Bp':
   path = '/store/user/jschulte/Pythia8_BuTau3mu_Run3_2022/SkimTau3mu_MCRun3_Bu_Mini_v4/221214_073101'
if args.dataset == 'MC' and args.MCprocess == 'B0':
   path = '/store/user/jschulte/Pythia8_BdTau3mu_Run3_2022/SkimTau3mu_MCRun3_Bd_Mini_v4/221214_073125'
if args.dataset == 'MC' and args.MCprocess == 'DsPhiMuNu':
    if args.anatype == "phimunu":
        path = '/store/user/jschulte/DstoPhiMuNu_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimDsPhiMuNu_MCRun3_DsPhiMuNu_Mini_v3/240725_195013/'
    else:
        path = '/store/user/jschulte/DstoPhiMuNu_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/SkimTau3mu_MCRun3_DsPhiMuNu_Mini_v3/240725_180330/'
        
if args.dataset == 'MC' and args.MCprocess == 'DsPhiPi':
   path = '/store/user/jschulte/Pythia8_DsPhiPi_MuMuPi_Run3_2022/SkimPhiPi_MCRun3_Mini_v4/221213_155107'


#generating the list of all .root files in given directory and subdirectories
fileList = []
for r, d, f in os.walk('/eos/purdue'+path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            fileList.append(os.path.join(r.split('/eos/purdue')[-1], file))
            print (file)

#prepare final script
#final_script = open("submit_analysis_"+startTime+".sh", "w")
final_script = open("submit_analysis_"+output_name+".sh", "w")
final_script.write("#!/bin/bash\n")
final_script.write("cd "+output_name+"\n")
final_script.write("export VOMS_PATH=$(echo $(voms-proxy-info | grep path) | sed 's/path.*: //')\n")
final_script.write("export VOMS_USERID=$(echo $(voms-proxy-info | grep path) | sed 's/.*p_u//')\n")
final_script.write("export VOMS_TRG=/home/$USER/x509up_u$VOMS_USERID\n")
final_script.write("cp $VOMS_PATH $VOMS_TRG\n")
final_script.write("export X509_USER_PROXY=$VOMS_TRG\n")




#loop to generate one .cpp+executable+batch system conf file for each group of "n" files
n_chunk = len(fileList)//args.n
print('Number of files is {0:2d}'.format(len(fileList)))
print('Number of jobs is {0:2d}'.format(n_chunk+1))
for file_index in range(n_chunk+1):
      chunk = '' 
      for idx, l in enumerate(fileList):
         if idx < args.n*(file_index+1) and idx >= args.n*file_index:
             l = l.rstrip()
             l = '        chain->AddFile("root://af-a00.cms.rcac.purdue.edu/{}");\n'.format(l)
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
      #with open("templates/my_HTCondor_template.job", "r") as myCondor_infile:
      #    buf3 = myCondor_infile.readlines()

      #condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      #with open(output_name+"/"+condor_filename, "w") as myCondor_outfile:
      #    for lb3 in buf3:
      #        if lb3 == "Executable = launch_analysis_template.job\n":
      #            myCondor_outfile.write("Executable = "+launch_filename+"\n")
      #        else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo sbatch --mem=4G -N1 -n1 --time=12:00:00 --account=cms-express "+launch_filename+" \n")
      final_script.write("sbatch --mem=4G -N1 -n1 --time=12:00:00 --account=cms-express "+launch_filename+" \n")

final_script.close()
#submitName = "submit_analysis_"+startTime+".sh"
#source submitName
