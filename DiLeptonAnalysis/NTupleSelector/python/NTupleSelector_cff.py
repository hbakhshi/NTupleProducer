import FWCore.ParameterSet.Config as cms

NTupleSelector = cms.EDFilter( 
    "NTupleSelector",
    min_nelectrons = cms.untracked.int32( -100 ),
    min_nmuons = cms.untracked.int32( -100 ),
    min_ntaus = cms.untracked.int32( -100 ),
    min_nleptons = cms.untracked.int32( -100 ),
    min_njets = cms.untracked.int32( -100 ),
    min_HT = cms.untracked.double( 0.0 ),
    HLTPaths = cms.vstring( ['',''] )
)
