from CRABClient.UserUtilities import config
config = config()

date = "260124"
requestNames = ["DoubleMu2016G", "DoubleMu2016H",
                "SingleMu2016G", "SingleMu2016H",
                "DoubleEl2016G", "DoubleEl2016H",
                "SingleEl2016G", "SingleEl2016H",
                "MuEG2016G", "MuEG2016H"]

jsonInputs = ["doublemu_2016G.json", "doublemu_2016H.json",
              "singlemu_2016G.json", "singlemu_2016H.json",
              "doubleel_2016G.json", "doubleel_2016H.json",
              "singleel_2016G.json", "singleel_2016H.json",
              "mueg_2016G.json", "mueg_2016H.json"]

datasets = ["/DoubleMuon/Run2016G-21Feb2020_UL2016-v1/AOD",
            "/DoubleMuon/Run2016H-21Feb2020_UL2016-v1/AOD",
            "/SingleMuon/Run2016G-21Feb2020_UL2016-v1/AOD",
            "/SingleMuon/Run2016H-21Feb2020_UL2016-v1/AOD",
            "/DoubleEG/Run2016G-21Feb2020_UL2016-v1/AOD",
            "/DoubleEG/Run2016H-21Feb2020_UL2016-v1/AOD",
            "/SingleElectron/Run2016G-21Feb2020_UL2016-v1/AOD",
            "/SingleElectron/Run2016H-21Feb2020_UL2016-v2/AOD",
            "/MuonEG/Run2016G-21Feb2020_UL2016-v1/AOD",
            "/MuonEG/Run2016H-21Feb2020_UL2016-v1/AOD"]

sampleid = 1
config.General.requestName = requestNames[sampleid] + "_crab" + date + "_skimmer"
config.General.workArea = "Outreach_CrabSkimmer"
config.General.transferLogs = False
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "cmssw_edm_event_skim.py"
config.JobType.inputFiles = [jsonInputs[sampleid]]
config.JobType.pyCfgParams = ["eventsJSON=" + jsonInputs[sampleid]]

config.Data.inputDataset = datasets[sampleid]
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
# config.Data.unitsPerJob = 5
# config.Data.totalUnits = 100
config.Data.publication = False
config.Data.outputDatasetTag = "Outreach_" + requestNames[sampleid] + "_" + date

config.Site.storageSite = "T2_UK_SGrid_RALPP"