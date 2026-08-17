// Build a model for the cluster energy fit

struct Component_t {
  TString name;
  TH1* h;
  double norm;
  bool is_signal;
  int color;
};

vector<Component_t> getComponents(const vector<Process_t>& processes, const char* name, const int set) {
  vector<Component_t> components;
  for(const auto& process : processes) {
    TH1* h = (TH1*) process.f->Get(Form("hist_%i/%s", set + process.set_offset, name));
    if(!h) {
      Error(__func__, "Could not retrieve histogram %s from process %s", name, process.name.Data());
      continue;
    }
    h = (TH1*) h->Clone(Form("h_%s_%i_%s", name, set, process.name.Data()));
    h->Scale(process.norm);
    components.push_back({process.name, h, process.norm, process.is_signal, process.color});
  }
  return components;
}

vector<Component_t> buildModel(const vector<Process_t>& processes, const char* name, const int set,
                               const double xmin = 60., const double xmax = 100.) {
  vector<Component_t> components = getComponents(processes, name, set);
  if(components.empty()) {
    Error(__func__, "No components found for model!");
  }

  // For each component, create a smooth model based on the expected shape
  for(auto& comp : components) {
    if(!comp.h) continue;
    // Smooth the histogram to create a model
    TH1* h_smooth = (TH1*) comp.h->Clone(Form("%s_smooth", comp.name.Data()));
    if(comp.name.Contains("Cosmic")) {
        // Assume cosmics are roughly flat
        const double rate = h_smooth->Integral(h_smooth->FindBin(xmin), h_smooth->FindBin(xmax));
        h_smooth->Reset();
        for(int i = 1; i <= h_smooth->GetNbinsX(); ++i) {
            h_smooth->SetBinContent(i, rate / (xmax - xmin)*h_smooth->GetBinWidth(i));
        }
    } else if(comp.is_signal) {
        // For signal, no need to smooth (stats are high)
    } else {
        // For all other backgrounds, assume an exponential falloff and smooth the histogram
        // Find the tail
        int bins_found = 0; int start_bin = 1;
        for(int bin = h_smooth->GetNbinsX(); bin > 0; --bin) {
            if(h_smooth->GetBinContent(bin) > 0) {
                ++bins_found;
                start_bin = bin;
                if(bins_found > 8) break; // enough bins for a git
            }
        }
        if(h_smooth->GetBinCenter(start_bin) < xmax) {
            TF1 f("f", "expo(0)", h_smooth->GetBinCenter(start_bin), h_smooth->GetBinCenter(h_smooth->GetNbinsX()));
            h_smooth->Fit(&f, "RX0");
            for(int bin = start_bin; bin <= h_smooth->GetNbinsX(); ++bin) {
                h_smooth->SetBinContent(bin, f.Eval(h_smooth->GetBinCenter(bin)));
            }
        }
    }
    comp.h = h_smooth;
  }

  return components;
}

void plotModel(const vector<Process_t>& processes, const char* name, const int set) {
  const double xmin = 60., xmax = 100.;
  vector<Component_t> components = buildModel(processes, name, set, xmin, xmax);
  if(components.empty()) return;

  // Create a canvas to plot the model
  TCanvas c("c", "Model", 800, 600);
  THStack h_stack("h_stack", "Model Components");
  TH1* h_signal = nullptr;
  TLegend legend(0.15, 0.70, 0.85, 0.89);
  legend.SetNColumns(3);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  const int rebin = 2;
  for(const auto& comp : components) {
    if(!comp.h) continue;
    comp.h->SetLineColor(comp.color);
    comp.h->Rebin(rebin);
    if(comp.is_signal) {
        comp.h->SetLineWidth(3);
        comp.h->SetFillStyle(3004);
        comp.h->SetFillColor(comp.color);
        if(!h_signal) {
            h_signal = (TH1*) comp.h->Clone("h_signal");
            legend.AddEntry(h_signal, comp.name, "F");
        } else {
            h_signal->Add(comp.h);
        }
    } else {
        comp.h->SetFillColor(comp.color);
        comp.h->SetFillStyle(kSolid);
        comp.h->SetLineColor(kBlack);
        comp.h->SetLineWidth(1);
        h_stack.Add(comp.h);
        legend.AddEntry(comp.h, comp.name, "F");
    }
  }
  if(!h_signal) {
    Error(__func__, "No signal component found for model!");
    return;
  }
  h_signal->Draw("hist");
  h_stack.Draw("hist same noclear");
  h_signal->Draw("hist same");
  legend.Draw();
  double max_val = std::max(h_signal->GetMaximum(), h_stack.GetMaximum());
  h_signal->GetYaxis()->SetRangeUser(0., 1.3*max_val);
  h_signal->GetXaxis()->SetRangeUser(65., xmax);
  c.SaveAs(Form("%s/model_%s_%i.png", dir_.Data(), name, set));
  h_signal->GetYaxis()->SetRangeUser(1.e-6*max_val, 200*max_val);
  c.SetLogy();
  c.SaveAs(Form("%s/model_%s_%i_log.png", dir_.Data(), name, set));
}