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


eras = ['2022F_0', '2022F_1', '2022F_2', '2022F_3', '2022F_4', '2022F_5', '2022F_6', '2022F_7', '2022G_0', '2022G_1', '2022G_2', '2022G_3', '2022G_4', '2022G_5', '2022G_6', '2022G_7']

if ANALYSISTYPE == "control":
    eras = ['2022F_0', '2022F_1', '2022F_2', '2022F_3', '2022F_4', '2022F_5', '2022F_6', '2022F_7', '2022G_0', '2022G_1', '2022G_2', '2022G_3', '2022G_4', '2022G_5', '2022G_6', '2022G_7']

for i in eras:

    tmpFiles = glob.glob("%s_%s_%s/*.root"%(i,ANALYSISTYPE,OUT_NAME))
    for f in tmpFiles:
        if "merged" in f: continue
        #fileString+=f
        #fileString+=" "
        files.append(f)

call = ["hadd", '-f', "AnalysedTree_%s_2022postEE_%s_merged_%s.root"%(DATATYPE, ANALYSISTYPE, OUT_NAME)] + files
print (call)
subprocess.call(call)
#hadd AnalysedTree_${DATATYPE}_${i}_${ANALYSISTYPE}_merged_${OUT_NAME}.root AnalysedTree_${DATATYPE}_${i}_${ANALYSISTYPE}*.root

