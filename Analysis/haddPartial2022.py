import subprocess, argparse, sys, glob

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--anaType')
parser.add_argument('-o', '--outName')

args = parser.parse_args()

ANALYSISTYPE=args.anaType
OUT_NAME=args.outName


# Control channel
if ANALYSISTYPE == "control":
    print ("DsPhiPi analysis")
    DATATYPE='data_control'
elif ANALYSISTYPE == "tau3mu":
    print ("Tau3mu analysis")
    DATATYPE='data'
elif ANALYSISTYPE == "phimunu":
    print ("Tau3mu phimunu")
    DATATYPE='data_phimunu'
else:
    print ("unsupported analysis type")
    sys.exit()

files = []

#for i in ['2022C_0', '2022C_1', '2022C_2', '2022C_3', '2022C_4', '2022C_5', '2022C_6', '2022C_7', '2022D-v1_0', '2022D-v1_1', '2022D-v1_2', '2022D-v1_3', '2022D-v1_4', '2022D-v1_5', '2022D-v1_6', '2022D-v1_7', '2022D-v2_0', '2022D-v2_1', '2022D-v2_2', '2022D-v2_3', '2022D-v2_4', '2022D-v2_5', '2022D-v2_6', '2022D-v2_7', '2022E_0', '2022E_1', '2022E_2', '2022E_3', '2022E_4', '2022E_5', '2022E_6', '2022E_7', '2022F_0', '2022F_1', '2022F_2', '2022F_3', '2022F_4', '2022F_5', '2022F_6', '2022F_7', '2022G_0', '2022G_1', '2022G_2', '2022G_3', '2022G_4', '2022G_5', '2022G_6', '2022G_7']:
for i in ['2022G_0', '2022G_1', '2022G_2', '2022G_3', '2022G_4', '2022G_5', '2022G_6', '2022G_7']:

    tmpFiles = glob.glob("%s_%s_%s/*.root"%(i,ANALYSISTYPE,OUT_NAME))
    for f in tmpFiles:
        if "merged" in f: continue
        #fileString+=f
        #fileString+=" "
        files.append(f)

call = ["hadd", '-f', "AnalysedTree_%s_2022G_%s_merged_%s.root"%(DATATYPE, ANALYSISTYPE, OUT_NAME)] + files
print (call)
subprocess.call(call)
#hadd AnalysedTree_${DATATYPE}_${i}_${ANALYSISTYPE}_merged_${OUT_NAME}.root AnalysedTree_${DATATYPE}_${i}_${ANALYSISTYPE}*.root

