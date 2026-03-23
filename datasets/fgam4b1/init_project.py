#!/usr/bin/python

from local_classes import *
import os

class Project(ProjectBase):
    def init_datasets(self):
#------------------------------------------------------------------------------
# init datasets - for stage 4, the input is a dataset produced at stage 3
#------------------------------------------------------------------------------
        self.add_dataset(Dataset('dig.mu2e.fgam4b1s41r0000.Run1BAna.art','fgam4b1s41r0000','local'));


    def __init__(self,idsid=None):
        familyID  = 'fgam4b1'
        user      = os.getenv('USER')

        ProjectBase.__init__(self,project='Run1BAna',family_id='fgam4b1',idsid=idsid);
        self.init_datasets();

#------------------------------------------------------------------------------
# s5:reco_trig_nt
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_trig_nt','fgam4b1s41r0000'); #idsid);

        job.fNInputFiles             = 1949                     # number of the job segments

        job.fMaxInputFilesPerSegment =  20                      # MC generator
        job.fMaxSegments             = 1000
        # job.fNEventsPerSegment       =  -1                    # defined by the input dataset
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        reco_version                 = 'r0002'
        job.fOutputStream            = [ 'defaultOutput'                    ]
        job.fOutputDsID              = [ familyID+s.name()+'1'+reco_version ]
        job.fOutputFnPattern         = [ 'nts.'+user+'.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'root'                             ]
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------
