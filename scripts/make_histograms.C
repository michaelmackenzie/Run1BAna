// Process the input datasets and write out histogram files

struct config_t {
  TString dataset;
  TString period;
  bool    process;
  config_t(TString name, TString period = "", bool proc = true) : dataset(name), period(period), process(proc) {}
};

int count_processes() {
  TString res = gSystem->GetFromPipe("ps -elf | grep ${USER} | grep make_histograms.C | grep -v grep | wc -l");
  return max(0, res.Atoi() - 1); // subtract one for the master process
}

/**
   processes: Number of parallel processes to submit
   dataset  : Specific dataset to process, if empty it will process all enabled datasets
   mode     : Histogramming mode
   function : TStnAna processing function, defined in ana/scripts/
 **/
int make_histograms(int processes = 1, TString dataset = "", const int mode = 1,
                    const char* function = "trg_ana") {

  if(processes > 3) {
    cout << "Requested " << processes << " parallel processes, but this exceeds the interactive maximum of about 2-3!\n";
    return 1;
  }
  vector<config_t> datasets = {
    config_t("cele0b1s5r0100", "", false), // MDC2025 samples
    config_t("cele1b1s5r0100", "", false),
    config_t("cpos0b1s5r0100", "", false),
    config_t("cpos1b1s5r0100", "", false),
    config_t("cry4ab1s5r0100", "", false),
    config_t("fele0b1s5r0100", "", false),
    config_t("fpos0b1s5r0100", "", false),
    config_t("pbar1b2s5r0100", "", false),
    config_t("mnbs0b1s5r0100", "", false),
    config_t("mnbs0b2s5r0100", "", false),
    config_t("cele0b1s5r0004", "", false), // MDC2020 samples
    config_t("cele1b1s5r0004", "",  true),
    config_t("cpos0b1s5r0004", "", false),
    config_t("cpos1b1s5r0004", "",  true),
    config_t("cry4ab1s5r0001", "", false),
    config_t("fele0b1s5r0001", "", false),
    config_t("fpos0b1s5r0000", "", false),
    config_t("pbar1b2s5r0000", "", false),
    config_t("mnbs0b1s5r0004", "", false),
    config_t("mnbs0b2s5r0004", "", false)
  };

  vector<TString> logs;
  for(auto config : datasets) {
    if(dataset == "" && !config.process) continue;
    if(dataset != "" && config.dataset != dataset) continue;
    if(processes > 1) {
      while(count_processes() >= processes) sleep(10); // sleep until a job finishes
      gSystem->Exec("[ ! -d log ] && mkdir log");
      TString command = Form("root.exe -q -b \"${MUSE_WORK_DIR}/ConvAna/scripts/make_trig_histograms.C(0, \\\"%s\\\", %i, \\\"%s\\\")\" >| log/trig_out_%s.log 2>&1 &",
                             config.dataset.Data(), mode, function, config.dataset.Data());
      printf(" Submitting %-20s histogramming...\n", config.dataset.Data());
      logs.push_back(Form("log/trig_out_%s.log", config.dataset.Data()));
      gSystem->Exec(command.Data());
    } else {
      stnana("ConvAna" ,config.dataset.Data(),"","",Form("ConvAna_%s()/save=ConvAna.%s.%s.hist",
                                                         function, function, config.dataset.Data()));
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
