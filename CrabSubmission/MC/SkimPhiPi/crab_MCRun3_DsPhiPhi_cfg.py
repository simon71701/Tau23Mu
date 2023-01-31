from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'SkimPhiPi_MCRun3_Mini_v3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'


config.JobType.psetName = '/lustrehome/aruta/Tau3mu-Run3_Analysis_12_4_11_patch3/CMSSW_12_4_11_patch3/src/SkimTools/SkimPhiPi/test/run_MC2022_DsPhiPiSkimAndTree_cfg.py'
config.Data.inputDataset = '/Pythia8_DsPhiPi_MuMuPi_Run3_2022/caruta-124X_mcRun3_2022_realistic_v12_MINIAODSIM-40e7f3028b76c221fdb27dc79a76e5ce/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 2500
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimPhiPi_MCRun3_Mini_v3'
config.JobType.allowUndistributedCMSSW = True 
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist  = True

