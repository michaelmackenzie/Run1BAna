//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script
//               compiler in Mu2e environment. The name of this macro file
//               is defined by the .rootrc file
//
// Jul 08 2014 P.Murat
//------------------------------------------------------------------------------
{
                                // the line below tells rootcling where to look for
                                // the include files

  gInterpreter->AddIncludePath("./build/include");
  gInterpreter->AddIncludePath(gSystem->Getenv("CLHEP_INC"));
  gInterpreter->AddIncludePath(Form("%s/include",gSystem->Getenv("ROOTSYS")));

//-----------------------------------------------------------------------------
//  check batch mode
//-----------------------------------------------------------------------------
  const char* opt ;
  int batch_mode = 0;

  int nargs = gApplication->Argc();

  for (int i=1; i<nargs; i++) {
    opt  = gApplication->Argv(i);
    if (strcmp(opt,"-b") == 0) {
      batch_mode = 1;
      break;
    }
  }

  printf("   batch_mode = %i\n",batch_mode);
//-----------------------------------------------------------------------------
// always need libStntuple_loop, but the other 2 libs should be loaded in
// only if we're running bare root
//-----------------------------------------------------------------------------
  const char* exec_name = gApplication->Argv(0);

  printf(" nargs = %2i exec_name = %s\n",nargs, exec_name);

  if (exec_name) {
    if (TString(exec_name) == "root.exe") {
//-----------------------------------------------------------------------------
// assume STNTUPLE  analysis job
//-----------------------------------------------------------------------------
      if (batch_mode == 1) gSystem->Load("$ROOTSYS/lib/root/libGui.so");
//-----------------------------------------------------------------------------
// Mu2e Offline libraries
//-----------------------------------------------------------------------------
//     //     gSystem->Load("lib/libmu2e_Mu2eInterfaces.so");
//     //     gSystem->Load("lib/libmu2e_CalorimeterGeom.so");
//
      cout << "Loading Stntuple libraries\n";
      gSystem->Load("$MUSE_BUILD_DIR/Stntuple/lib/libStntuple_base.so");
      gSystem->Load("$MUSE_BUILD_DIR/Stntuple/lib/libStntuple_obj.so");
      gSystem->Load("$MUSE_BUILD_DIR/Stntuple/lib/libStntuple_loop.so");
      gSystem->Load("$MUSE_BUILD_DIR/Stntuple/lib/libStntuple_alg.so");
      gSystem->Load("$MUSE_BUILD_DIR/Stntuple/lib/libStntuple_ana.so");
      gSystem->Load("$MUSE_BUILD_DIR/Stntuple/lib/libStntuple_val.so");

      cout << "Loading analysis libraries\n";
                                        //insert_user_libs_here
                                        // print overflows/underflows in the stat box
      gStyle->SetOptStat(11111111);
                                        // print fit results in the stat box
      gStyle->SetOptFit(1110);
      TArrow::SetDefaultArrowSize(0.015);
    }
//-----------------------------------------------------------------------------
//  databases
//-----------------------------------------------------------------------------
//   gSystem->Load("libStntuple_oracle.so");

    if (gSystem->Exec("ls $HOME/root_macros/set_style.C &> /dev/null") == 0) {
      gInterpreter->ExecuteMacro("$HOME/root_macros/set_style.C");
    }
  }
//-----------------------------------------------------------------------------
// report the process ID which simplifies debugging
//-----------------------------------------------------------------------------
  printf(" process ID: %i\n",gSystem->GetPid());
  TAuthenticate::SetGlobalUser(gSystem->Getenv("USER"));
  gInterpreter->ProcessLine(".! ps | grep root");
}
