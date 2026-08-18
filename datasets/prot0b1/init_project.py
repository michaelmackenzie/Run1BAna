#!/usr/bin/python

from local_classes import *
import os

class Project(ProjectBase):
    def init_datasets(self):
#------------------------------------------------------------------------------
# init datasets - for stage 4, the input is a dataset produced at stage 3
#------------------------------------------------------------------------------
        self.add_dataset(Dataset('sim.mu2e.prot0b1s21r0000.Run1BAna.art','prot0b1s21r0000','local'));
        self.add_dataset(Dataset('dig.mu2e.prot0b1s41r0003.Run1BAna.art','prot0b1s41r0003','local'));


    def __init__(self,idsid=None):
        familyID  = 'prot0b1'
        user      = os.getenv('USER')

        ProjectBase.__init__(self,project='Run1BAna',family_id='prot0b1',idsid=idsid);
        self.init_datasets();

#------------------------------------------------------------------------------
# s3:sim
#------------------------------------------------------------------------------
        s                            = self.new_stage('s3');
        job                          = s.new_job('sim', 'prot0b1s21r0000');

        job.fNInputFiles             = -1                       # number of the job segments

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       = 100000                   # 
        job.fMaxSegments             = 100
        job.fResample                = 'yes'                    # yes/no
        job.fResamplingModuleLabel   = 'TargetStopResampler'
        job.fRunNumber               = 1470
        job.fMaxMemory               = '4000MB'
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'ifdh'                   # ifdh/xrootd
        job.fOutputPath              = [ 'EndPath' ]

        reco_version                 = 'r0003'
        job.fOutputStream            = [ 'PrimaryOutput'                    ]
        job.fOutputDsID              = [ familyID+s.name()+'1'+reco_version ]
        job.fOutputFnPattern         = [ 'dts.'+user+'.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                              ]
#------------------------------------------------------------------------------
# s5:reco_trig_nt
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_trig_nt','prot0b1s41r0003'); #idsid);

        job.fNInputFiles             = 200                      # number of the job segments

        job.fMaxInputFilesPerSegment =  10                      # MC generator
        job.fMaxSegments             = 1000
        # job.fNEventsPerSegment       =  -1                    # defined by the input dataset
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        reco_version                 = 'r0003'
        job.fOutputStream            = [ 'defaultOutput'                    ]
        job.fOutputDsID              = [ familyID+s.name()+'1'+reco_version ]
        job.fOutputFnPattern         = [ 'nts.'+user+'.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'root'                             ]
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------
