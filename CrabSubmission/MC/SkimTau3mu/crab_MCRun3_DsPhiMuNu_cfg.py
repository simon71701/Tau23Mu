from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'SkimTau3mu_MCRun3_DsPhiMuNu_Miniv4_2022_NewSamples_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'


config.JobType.psetName = '/home/schul105/depot/Tau3Mu/analysis/el8/CMSSW_13_0_21/src/SkimTools/SkimTau3Mu/test/run_MC2022_PatAndTree_cfg.py'
config.Data.inputDataset = '/DstoPhiMuNu_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-FilterFix_130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 2500
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimTau3mu_MCRun3_DsPhiMuNu_Miniv4_2022_NewSamples_v1'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_US_Purdue'
config.Site.ignoreGlobalBlacklist  = True

