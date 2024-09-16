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

for i in ['2023C-v1_0', '2023C-v1_1', '2023C-v1_2', '2023C-v1_3', '2023C-v1_4', '2023C-v1_5', '2023C-v1_6', '2023C-v1_7', '2023C-v2_0', '2023C-v2_1', '2023C-v2_2', '2023C-v2_3', '2023C-v2_4', '2023C-v2_5', '2023C-v2_6', '2023C-v2_7', '2023C-v3_0', '2023C-v3_1', '2023C-v3_2', '2023C-v3_3', '2023C-v3_4', '2023C-v3_5', '2023C-v3_6', '2023C-v3_7', '2023C-v4_0', '2023C-v4_1', '2023C-v4_2', '2023C-v4_3', '2023C-v4_4', '2023C-v4_5', '2023C-v4_6', '2023C-v4_7', '2023D-v1_0', '2023D-v1_1', '2023D-v1_2', '2023D-v1_3', '2023D-v1_4', '2023D-v1_5', '2023D-v1_6', '2023D-v1_7', '2023D-v2_0', '2023D-v2_1', '2023D-v2_2', '2023D-v2_3', '2023D-v2_4', '2023D-v2_5', '2023D-v2_6', '2023D-v2_7']:

    tmpFiles = glob.glob("%s_%s_%s/*.root"%(i,ANALYSISTYPE,OUT_NAME))
    for f in tmpFiles:
        if "merged" in f: continue
        #fileString+=f
        #fileString+=" "
        files.append(f)

call = ["hadd", '-f', "AnalysedTree_%s_2023_%s_merged_%s.root"%(DATATYPE, ANALYSISTYPE, OUT_NAME)] + files
print (call)
subprocess.call(call)
#hadd AnalysedTree_${DATATYPE}_${i}_${ANALYSISTYPE}_merged_${OUT_NAME}.root AnalysedTree_${DATATYPE}_${i}_${ANALYSISTYPE}*.root

