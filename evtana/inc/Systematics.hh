//
// Run1BEvtAna systematics definitions
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_SYSTEMATICS_HH
#define RUN1BEVTANA_SYSTEMATICS_HH

//c++ includes
#include <map>
#include <iostream>

//ROOT includes
#include "TString.h"

//local includes
#include "Run1BAna/evtana/inc/GlobalConstants.h"

namespace Run1BEvtAna {
  class Systematics {
  public:

    //-----------------------------------------------------------------------------
    //  data members
    //-----------------------------------------------------------------------------
  public:

    // Define a systematic input
    struct Data_t {
      int     num_;
      TString name_;
      bool    up_;
      Data_t(int num = -1, TString name = "", bool up = true) : num_(num), name_(name), up_(up) {}
    };

    std::map<int, Data_t> idToData_;
    std::map<TString, Data_t> nameToData_;
    //-----------------------------------------------------------------------------
    //  functions
    //-----------------------------------------------------------------------------
  public:
    Systematics() {
      //initialize the defined systematics information
      for(int isys = 0; isys < kMaxSystematics; ++isys) {
        Data_t data = getData(isys);
        if(data.name_ == "") continue; //only store defined systematics
        idToData_[isys] = data;
        if(data.up_) nameToData_[data.name_] = data; //map to the up value
      }
    }
    ~Systematics() {}

    TString GetName(const int isys) {
      if(idToData_.find(isys) != idToData_.end()) return idToData_[isys].name_;
      return TString("");
    }

    int GetNum(TString name) {
      if(nameToData_.find(name) != nameToData_.end()) return nameToData_[name].num_;
      std::cout << "Systematics::" << __func__ << ": Undefined systematic " << name.Data() << std::endl;
      return -1;
    }

    bool IsUp(const int isys) {
      if(idToData_.find(isys) != idToData_.end()) return idToData_[isys].up_;
      return false;
    }
  private:
    static Data_t getData(const int isys) {
      switch(isys) {
      case   0: return Data_t(isys, "Nominal"       , true );
      case   1: return Data_t(isys, "Scale"         , true );
      case   2: return Data_t(isys, "Scale"         , false);
      default: break;
      }

      //return the default result if no systematic is defined
      return Data_t();
    }
  };
}
#endif
