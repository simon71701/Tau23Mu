import FWCore.ParameterSet.Config as cms

process = cms.Process('Tau3MuSkim')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load("SkimTools.SkimTau3Mu.Tau3MuSkim_miniAOD_cff")
process.load("SkimTools.SkimTau3Mu.Tau3MuSkim_miniAOD_noHLT_TriggerStudy-Run3_cff")

#process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
#process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v16' #data2018D
process.GlobalTag.globaltag = '124X_dataRun3_Prompt_v4' #data2022C

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Begin processing the 25271st record. Run 320012, Event 56448719, LumiSection 36 on stream 0 at 20-Apr-2020 18:53:30.862 CEST
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#Run2022E/EphemeralZeroBias0
        #'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/542/00000/58159c29-c7f6-4efb-85b1-555ccad24889.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/567/00000/d74d517e-b2ac-41a5-a1e0-352998f02182.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/661/00000/073cec82-8f9b-460a-b693-84b98cfdb358.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/661/00000/a4c60078-eae2-45d7-9d3a-5e3b8254d213.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/661/00000/bbd4effa-e1d4-437e-ac37-90bfadf75e61.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/762/00000/75ab2492-6fb6-4b96-b828-88aa28fb5697.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/871/00000/76c8ec4f-d87e-4d02-b544-7e525e165e11.root', 'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/359/998/00000/d21e3b8d-a35f-4f59-88ab-660b1478856b.root',

	#Run2022C/EphemeralZeroBias0
	#'root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/355/869/00000/7f15c759-5c1d-43d5-bb2b-a8601b641ee4.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/355/870/00000/169caf4b-4a20-4af2-a668-2cfa0a6b864c.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/355/871/00000/fea34784-0a04-4826-ad48-02c31936b4e7.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/355/872/00000/a5019361-9a58-45f8-bfd9-b9af8f380b17.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/355/933/00000/02956603-ab26-4e46-8b3e-6108ee966fb1.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralZeroBias0/MINIAOD/PromptReco-v1/000/355/933/00000/2a803da4-7b84-4f57-9c6b-8105e1c47b9e.root',

	#Run2022C/EphemeralHLTPhysics
	#'root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralHLTPhysics0/MINIAOD/PromptReco-v1/000/355/807/00000/bb1f9410-8fbc-4ca1-9792-dc5a72028768.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralHLTPhysics0/MINIAOD/PromptReco-v1/000/355/862/00000/3c0cc909-1614-49d4-b8f9-7323f2a03889.root','root://xrootd-cms.infn.it///store/data/Run2022C/EphemeralHLTPhysics0/MINIAOD/PromptReco-v1/000/355/863/00000/579808a0-5f9b-4de7-82b7-9c43035e3f45.root',

	#Run2022E/EphemeralHLTPhysics
	'root://xrootd-cms.infn.it///store/data/Run2022E/EphemeralHLTPhysics0/MINIAOD/PromptReco-v1/000/359/661/00000/01cf7db0-14b6-414c-9724-f9f5c813fd3e.root',


	#'root://xrootd-cms.infn.it///store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/355/863/00000/389f9ca1-f590-4691-b7f2-41e0146a8a79.root',
    ),
            #eventsToProcess = cms.untracked.VEventRange('320012:56448719')
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TreeData_triggerStudy.root"))


process.unpackedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
    patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ),
    triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
    unpackFilterLabels = cms.bool(True)
)

process.TreeMakerBkg = cms.EDAnalyzer("MiniAnaTau3Mu",
                                      isMcLabel = cms.untracked.bool(False),
                                      isAnaLabel = cms.untracked.bool(True),
                                      is2016Label = cms.untracked.bool(False),
                                      is2017Label = cms.untracked.bool(False),
                                      is2018Label = cms.untracked.bool(True),
                                      isBParkingLabel = cms.untracked.bool(False),
                                      muonLabel=cms.InputTag("looseMuons"),
                                      photonLabel=cms.InputTag("slimmedPhotons"),
                                      VertexLabel=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                      genParticleLabel=cms.InputTag("prunedGenParticles"),
                                      pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
                                      Cand3MuLabel=cms.InputTag("ThreeMuonsVtxKalmanFit"),
                                      triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                      objects = cms.InputTag("unpackedPatTrigger"),
                                      AlgInputTag = cms.InputTag( "gtStage2Digis" ),
                                      algInputTag = cms.InputTag( "gtStage2Digis" ),
                                      extInputTag = cms.InputTag( "gtStage2Digis" )
                                      
)




process.Tau3MuSkim = cms.Path(#process.ThreeMuonSelSeq
			     process.ThreeMuonSelSeq*
                             process.unpackedPatTrigger*
                             process.TreeMakerBkg
                     )





