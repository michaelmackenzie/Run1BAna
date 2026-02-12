// Process the input datasets and write out histogram files
#include "Mu2eEvtAna/scripts/datasets.C"
#include "Mu2eEvtAna/scripts/functions.C"

// Count histogramming processes currently running
int count_processes() {
  TString res = gSystem->GetFromPipe("ps -elf | grep ${USER} | grep make_histograms.C | grep -v grep | wc -l");
  return max(0, res.Atoi() - 1); // subtract one for the master process
}

/**
   processes: Number of parallel processes to submit
   dataset  : Specific dataset to process, if empty it will process all enabled datasets
   mode     : Histogramming mode
   function : histogramming processing function, defined in ana/scripts/
 **/
int make_histograms(int processes = 1, TString dataset = "", const int mode = 1,
                    const char* function = "mu2e_ana") {

  if(processes > 3) {
    cout << "Requested " << processes << " parallel processes, but this exceeds the interactive maximum of about 2-3!\n";
    return 1;
  }

  // All defined datasets
  auto datasets = DATA::datasets();

  vector<TString> logs;
  for(auto config : datasets) {
    if(dataset == "" && !config.process_) continue;
    if(dataset != "" && config.name_ != dataset) continue;
    if(processes > 1) {
      while(count_processes() >= processes) sleep(10); // sleep until a job finishes
      gSystem->Exec("[ ! -d log ] && mkdir log");
      TString command = Form("root.exe -q -b \"${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/make_histograms.C(0, \\\"%s\\\", %i, \\\"%s\\\")\" >| log/out_%s.log 2>&1 &",
                             config.name_.Data(), mode, function, config.name_.Data());
      printf(" Submitting %-20s histogramming...\n", config.name_.Data());
      logs.push_back(Form("log/out_%s.log", config.name_.Data()));
      gSystem->Exec(command.Data());
    } else {
      // Process this dataset
      gInterpreter->ProcessLine(Form("%s(\"%s\", %i);",
                                     function, config.name_.Data(), mode));
    }
  }

  // wait for the jobs to finish
  if(processes > 1) {
    while(count_processes() > 0) sleep(10); // sleep until all jobs finish
  }
  printf("Finished histogramming!\n");

  // for(auto log : logs) {
  //   gSystem->Exec(Form("echo %s; tail -n 1 %s", log.Data(), log.Data()));
  // }

  return 0;
}
