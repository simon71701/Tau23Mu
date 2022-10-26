from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'SkimDsPhiPi_2022eraC_stream2_Mini_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'


config.JobType.psetName = '/lustrehome/aruta/Tau3mu-Run3_Analysis/CMSSW_12_4_7/src/SkimTools/SkimPhiPi/test/run_Data2022_DsPhiPiSkimAndTree_cfg.py'

config.Data.inputDataset = '/ParkingDoubleMuonLowMass2/Run2022C-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 150
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_357900_Golden.json'
#config.Data.runRange = '193093-193999' # '193093-194075'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimDsPhiPi_2022eraC_stream2_Mini_v2'
config.JobType.allowUndistributedCMSSW = True 
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist  = True

