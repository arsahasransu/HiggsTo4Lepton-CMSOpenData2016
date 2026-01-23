import json

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


def load_eventid_strings(json_path):
    with open(json_path, 'r') as fh:
        data = json.load(fh)

    out = []
    if isinstance(data, list) and data and all(isinstance(x, int) for x in data):
        for ev in data:
            out.append("0:0:{0}".format(int(ev)))
        return out

    if isinstance(data, list) and all(isinstance(x, dict) for x in data):
        for d in data:
            evt = int(d.get('event') or d.get('Event') or 0)
            run = int(d.get('run') or d.get('Run') or 0)
            lumi = int(d.get('lumi') or d.get('luminosityBlock') or d.get('LuminosityBlock') or 0)
            out.append("{0}:{1}:{2}".format(run, lumi, evt))
        return out

    raise ValueError('Unsupported JSON format for events. Use list of ints or list of dicts.')


# parse command-line options via VarParsing so cmsRun can override
options = VarParsing('analysis')
options.register('eventsJSON', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'JSON file with events')

options.parseArguments()

if not options.eventsJSON:
    raise RuntimeError('ERROR: please provide eventsJSON=path/to/events.json')

eventid_list = load_eventid_strings(options.eventsJSON)

process = cms.Process('SKIM')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Configure MessageLogger to report progress every 10000 events
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring("root://xrootd-cms.infn.it///store/data/Run2016H/DoubleMuon/AOD/21Feb2020_UL2016-v1/250000/3F1F25E0-A088-0049-9A19-A64FF39CAAD0.root"),
    eventsToProcess = cms.untracked.VEventRange(*eventid_list)
)

process.p = cms.Path()  # no processing modules required for pure event selection

process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string("skimmed_AOD.root"),
    outputCommands = cms.untracked.vstring('keep *')
)

process.end = cms.EndPath(process.out)

