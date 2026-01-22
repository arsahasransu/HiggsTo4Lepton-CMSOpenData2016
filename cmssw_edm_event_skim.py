#!/usr/bin/env python
"""
Generated with Copilot GPT
cmsRun config to skim MINIAOD/miniAOD ROOT files by event numbers specified in a JSON file.

This configuration reads a JSON file (list of ints or list of dicts with 'run','lumi','event'),
converts entries to CMS EventIDs and assigns them to `process.source.eventsToProcess` so that
cmsRun only processes those events and the `PoolOutputModule` writes them to a new EDM file.

Usage example:
  cmsRun miniaod_event_skim.py inputFiles="file1.root file2.root" eventsJSON=events.json outputFile=selected.root

JSON formats supported:
  - [12345, 23456, 34567]
  - [{"event":12345}, {"event":23456}]
  - [{"run":1,"lumi":2,"event":12345}, ...]

Notes:
  - When JSON contains only event numbers, run and lumi are set to 0 (useful when you only care about event id).
  - The output file will contain all branches (outputCommands = 'keep *').
"""

import json
import sys
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms


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
options.register('inputList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Text file containing input ROOT file paths, one per line')

options.parseArguments()

if not options.eventsJSON:
    print('ERROR: please provide eventsJSON=path/to/events.json')
    sys.exit(1)

eventid_list = load_eventid_strings(options.eventsJSON)

process = cms.Process('SKIM')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Configure MessageLogger to report progress every 10000 events
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

infiles = []
# If an input list file is provided, read file paths from it (one per line), ignore blank lines and comments
if getattr(options, 'inputList', ''):
    try:
        with open(options.inputList, 'r') as lf:
            for line in lf:
                s = line.strip()
                if not s:
                    continue
                if s.startswith('#'):
                    continue
                infiles.append(s)
    except Exception as e:
        print('ERROR: failed to read inputList file {0}: {1}'.format(options.inputList, e))
        sys.exit(1)
elif options.inputFiles:
    infiles = options.inputFiles
else:
    print('ERROR: provide inputFiles="file1.root file2.root" or inputList=/path/to/list.txt')
    sys.exit(1)

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(*infiles),
    eventsToProcess = cms.untracked.VEventRange(*eventid_list)
)

process.p = cms.Path()  # no processing modules required for pure event selection

process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring('keep *')
)

process.end = cms.EndPath(process.out)

