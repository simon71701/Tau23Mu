from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'SkimTau3Mu_MC2022_Ds_Mini_v2bis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'


config.JobType.psetName = '/lustrehome/aruta/Tau3mu-Run3_Analysis/CMSSW_12_4_7/src/SkimTools/SkimTau3Mu/test/run_MC2022_PatAndTree_cfg.py'

config.Data.inputDataset = '/Pythia8_DsToTauTo3Mu_Run3_2022_GEN-SIM/caruta-MCRun3_ToTauTo3Mu_MINIAODSIM-680117774a8e2ab3f00b0e4cd4298b53/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 9999
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimTau3Mu_MC2022_Ds_Mini_v2bis'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist  = True
