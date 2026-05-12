// Optimize the hole shape for the Run 1B target + Run 1A running

TString fig_dir_ = "figures";
bool smooth_inputs_ = true; // Smooth the 2D distributions

//--------------------------------------------------------------------------------------------------------------------
// Evaluate Run 1A metric
double evaluate_run1a_metric(const double run1a_stops) {
  return run1a_stops; // assume background is constant, so just maximize the signal rate
}

//--------------------------------------------------------------------------------------------------------------------
// Evaluate Run 1B metric
double evaluate_run1b_metric(const double run1b_muons, const double calo_stops) {
  constexpr double dio_eff = 1.e-3; // Approximate efficiency for Run 1B target -> calo disk 0
  return (run1b_muons <= 0.) ? 0. : dio_eff*run1b_muons / sqrt(std::pow(run1b_muons*dio_eff, 2) + 2.*dio_eff*run1b_muons*calo_stops + std::pow(calo_stops, 2));
}

//--------------------------------------------------------------------------------------------------------------------
// Evaluate the metrics for a given set of histograms
std::pair<double,double> evaluate_metrics(TH2* h_run1a_stops, TH2* h_run1b_muons, TH2* h_run1b_calo_stops) {
  const double calo_stops  = (h_run1b_calo_stops ? h_run1b_calo_stops->Integral() : 0.);
  const double run1b_muons = h_run1b_muons     ->Integral();
  const double run1a_stops = h_run1a_stops     ->Integral();

  const double run1b_metric = evaluate_run1b_metric(run1b_muons, calo_stops);
  const double run1a_metric = evaluate_run1a_metric(run1a_stops);

  return std::make_pair(run1a_metric, run1b_metric);
}

//--------------------------------------------------------------------------------------------------------------------
// Metric-ordered bin scan
TGraph* metric_ordered_bin_scan(TH2* h_run1a_stops, TH2* h_run1b_muons, TH2* h_run1b_calo_stops,
                                const double run1a_metric, const double run1b_metric, TGraph* g_circle_hole_metric) {
  TH2* h_run1a_metrics      = (TH2*) h_run1a_stops     ->Clone("h_run1a_metrics");
  TH2* h_run1b_metrics      = (TH2*) h_run1b_muons     ->Clone("h_run1b_metrics");
  TH2* total_metrics        = (TH2*) h_run1a_stops     ->Clone("total_metrics");
  for(int ix = 1; ix <= h_run1a_stops->GetNbinsX(); ++ix) {
    for(int iy = 1; iy <= h_run1a_stops->GetNbinsY(); ++iy) {
        const double run1a_stops = h_run1a_stops->GetBinContent(ix, iy);
        const double run1b_muons = h_run1b_muons->GetBinContent(ix, iy);
        const double calo_stops  = h_run1b_calo_stops->GetBinContent(ix, iy);
        const double run1a_metric_bin = evaluate_run1a_metric(run1a_stops) / run1a_metric; // Relative to optimal
        const double run1b_metric_bin = evaluate_run1b_metric(run1b_muons, calo_stops);
        const double total_metric_bin = (calo_stops > 0.) ? run1a_stops / calo_stops : (run1a_stops > 0.) ? 10. : 0.;
        h_run1a_metrics->SetBinContent(ix, iy, run1a_metric_bin);
        h_run1b_metrics->SetBinContent(ix, iy, run1b_metric_bin);
        total_metrics->SetBinContent(ix, iy, total_metric_bin);
    }
  }

  //---------------------------------------------------------------------
  // Plot the histograms
  //---------------------------------------------------------------------

  TCanvas c_metrics("c_metrics", "Metrics", 1000, 1000);
  h_run1a_stops->Draw("COLZ");
  c_metrics.SaveAs(fig_dir_ + "/run1a_stops.png");
  h_run1a_metrics->Draw("COLZ");
  c_metrics.SaveAs(fig_dir_ + "/run1a_metric_bins.png");
  h_run1b_muons->Draw("COLZ");
  c_metrics.SaveAs(fig_dir_ + "/run1b_muons.png");
  h_run1b_metrics->Draw("COLZ");
  c_metrics.SaveAs(fig_dir_ + "/run1b_metric_bins.png");
  total_metrics->Draw("COLZ");
  total_metrics->SetTitle("Run 1A Stops / Run 1B Calo Stops by bin;X (mm);Y (mm)");
  c_metrics.SaveAs(fig_dir_ + "/total_metric_bins.png");
  h_run1b_calo_stops->Draw("COLZ");
  c_metrics.SaveAs(fig_dir_ + "/run1b_calo_stops.png");

  //---------------------------------------------------------------------
  // Add hole bins ordered by total metric
  //---------------------------------------------------------------------

  TGraph* g_run1a_metric_vs_run1b_metric = new TGraph();
  TH2* h_run1a_stops_hole      = (TH2*) h_run1a_stops     ->Clone("h_run1a_stops_hole");
  TH2* h_run1b_muons_hole      = (TH2*) h_run1b_muons     ->Clone("h_run1b_muons_hole");
  TH2* h_run1b_calo_stops_hole = (TH2*) h_run1b_calo_stops->Clone("h_run1b_calo_stops_hole");
  h_run1a_stops_hole->Reset();
  h_run1b_calo_stops_hole->Reset();
  TH2* h_run1a_best = nullptr;
  TH2* h_run1b_best = nullptr;
  TH2* h_calo_stops_best = nullptr;
  const int max_bins = std::min(1000, h_run1a_stops->GetNbinsX() * h_run1a_stops->GetNbinsY()); // Limit to top 1000 bins for speed
  for(int bin = 1; bin <= max_bins; ++bin) {
    // Get maximum bin
    const int max_bin = total_metrics->GetMaximumBin();
    // Add the bin as a hole
    h_run1a_stops_hole->SetBinContent(max_bin, h_run1a_stops->GetBinContent(max_bin));
    h_run1b_calo_stops_hole->SetBinContent(max_bin, h_run1b_calo_stops->GetBinContent(max_bin));
    h_run1b_muons_hole->SetBinContent(max_bin, 0.); // not stopped in the plate
    // Re-evaluate the metrics
    auto [run1a_metric_step, run1b_metric_step] = evaluate_metrics(h_run1a_stops_hole, h_run1b_muons_hole, h_run1b_calo_stops_hole);
    g_run1a_metric_vs_run1b_metric->AddPoint(run1a_metric_step / run1a_metric, run1b_metric_step / run1b_metric);
    // Set the bin to 0 in the total metric to avoid picking it again
    total_metrics->SetBinContent(max_bin, 0.);
    if(!h_run1a_best && run1a_metric_step > 0.5*run1a_metric) {
      h_run1a_best = (TH2*) h_run1a_stops_hole->Clone("h_run1a_best");
      h_run1b_best = (TH2*) h_run1b_muons_hole->Clone("h_run1b_best");
      h_calo_stops_best = (TH2*) h_run1b_calo_stops_hole->Clone("h_calo_stops_best");
    }
  }

  //---------------------------------------------------------------------
  // Plot the Run 1A and Run 1B metrics
  //----------------------------------------------------------------------
  TCanvas c_metric_correlation("c_metric_correlation", "Run 1A vs Run 1B Metric", 800, 600);
  g_run1a_metric_vs_run1b_metric->Draw("AP");
  g_run1a_metric_vs_run1b_metric->SetMarkerStyle(20);
  g_run1a_metric_vs_run1b_metric->SetMarkerSize(0.8);
  g_run1a_metric_vs_run1b_metric->SetMarkerColor(kBlack);
  g_run1a_metric_vs_run1b_metric->SetTitle("Run 1A vs Run 1B Metric;Run 1A Metric (relative);Run 1B Metric (relative)");
  g_circle_hole_metric->SetLineColor(kRed);
  g_circle_hole_metric->SetLineWidth(2);
  g_circle_hole_metric->Draw("L SAME");
  TLegend leg(0.6, 0.75, 0.89, 0.89);
  leg.SetBorderSize(0); leg.SetFillStyle(0);
  leg.AddEntry(g_run1a_metric_vs_run1b_metric, "Metric-ordered Scan", "p");
  leg.AddEntry(g_circle_hole_metric, "Hole Radius Scan", "l");
  leg.Draw();
  c_metric_correlation.SaveAs(fig_dir_ + "/metric_ordered_correlation.png");

  if(h_run1a_best && h_run1b_best && h_calo_stops_best) {
    // Calculate the effective hole area/center
    double x(0.), y(0.), area(0.); int n_bins_hole(0);
    for(int ix = 1; ix <= h_run1a_stops->GetNbinsX(); ++ix) {
      for(int iy = 1; iy <= h_run1a_stops->GetNbinsY(); ++iy) {
        if(h_run1a_best->GetBinContent(ix, iy) > 0.) {
          double bin_x = h_run1a_best->GetXaxis()->GetBinCenter(ix);
          double bin_y = h_run1a_best->GetYaxis()->GetBinCenter(iy);
          x += bin_x + 3904.; // Shift to solenoid center
          y += bin_y;
          const double bin_width = h_run1a_best->GetXaxis()->GetBinWidth(ix);
          const double bin_height = h_run1a_best->GetYaxis()->GetBinWidth(iy);
          area += bin_width * bin_height;
          n_bins_hole++;
        }
      }
    }
    if(n_bins_hole > 0) {
      x /= n_bins_hole;
      y /= n_bins_hole;
    }
    const double radius = std::sqrt(area / M_PI);

    TEllipse* hole_ellipse = new TEllipse(x-3904., y, radius);
    hole_ellipse->SetLineColor(kRed);
    hole_ellipse->SetLineWidth(1);
    hole_ellipse->SetFillStyle(0);

    TCanvas c_best("c_best", "Best Hole", 1200, 400);
    c_best.Divide(3, 1);
    auto pad = c_best.cd(1); h_run1a_best->Draw("COLZ"); c_best.Update(); pad->SetGrid();
    hole_ellipse->Draw("SAME");
    c_best.cd(2); h_run1b_best->Draw("COLZ"); c_best.Update();
    hole_ellipse->Draw("SAME");
    c_best.cd(3); h_calo_stops_best->Draw("COLZ"); c_best.Update();
    hole_ellipse->Draw("SAME");
    c_best.SaveAs(fig_dir_ + "/best_hole.png");

    // Print the best hole rates/information
    const double run1a_stops      = h_run1a_stops     ->Integral();
    const double run1b_muons      = h_run1b_muons     ->Integral();
    const double run1b_calo_stops = h_run1b_calo_stops->Integral();
    const double best_calo_stops  = h_calo_stops_best->Integral();
    const double best_run1b_muons = h_run1b_best->Integral();
    const double best_run1a_stops = h_run1a_best->Integral();
    const double best_run1b_metric = evaluate_run1b_metric(best_run1b_muons, best_calo_stops);
    const double best_run1a_metric = evaluate_run1a_metric(best_run1a_stops);
    printf("-----------------------------------------------------------------\n");
    printf("Best hole metric information:\n");
    printf("Run 1A Stops     : %.3e (%.1f%% of optimal)\n", best_run1a_stops, 100.*best_run1a_stops/run1a_stops);
    printf("Run 1B Muons     : %.3e (%.1f%% of optimal)\n", best_run1b_muons, 100.*best_run1b_muons/run1b_muons);
    printf("Run 1B Calo Stops: %.3e (%.1f%% removed)\n", best_calo_stops, 100.*(1. - best_calo_stops/run1b_calo_stops));
    printf("Run 1A metric    : %.3e (%.1f%% of optimal)\n", best_run1a_metric, 100.*best_run1a_metric/run1a_metric);
    printf("Run 1B metric    : %.3e (%.1f%% of optimal)\n", best_run1b_metric, 100.*best_run1b_metric/run1b_metric);
    printf("Effective hole area: %.1f mm^2 --> radius = %.2f\n", area, radius);
    printf("Effective hole center: (%.1f, %.1f) mm\n", x, y);
    printf("-----------------------------------------------------------------\n");
  }
  return nullptr;
}

//--------------------------------------------------------------------------------------------------------------------
// Hole radius optimization for the Run 1B target
TGraph* optimize_hole_radius(const TH2* h_run1a_stops, const TH2* h_run1b_muons, const TH2* h_run1b_calo_stops,
                             const double run1a_metric, const double run1b_metric) {
  const double r_min =   0.; // Minimum hole radius (mm)
  const double r_max = 100.; // Maximum hole radius (mm)
  const int n_bins   = 100 ; // Number of steps in radius to take
  TH1* h_run1a_metric_vs_r = new TH1D("h_run1a_metric_vs_r", "Run 1A Metric vs Hole Radius;Hole Radius (mm);Run 1A Metric", n_bins, r_min, r_max);
  TH1* h_run1b_metric_vs_r = new TH1D("h_run1b_metric_vs_r", "Run 1B Metric vs Hole Radius;Hole Radius (mm);Run 1B Metric", n_bins, r_min, r_max);
  TGraph* g_run1a_metric_vs_run1b_metric = new TGraph();

  for(int istep = 0; istep < n_bins; ++istep) {
    double r = r_min + (r_max - r_min) * (istep + 0.5) / n_bins;
    TH2* h_run1a_stops_hole      = (TH2*) h_run1a_stops     ->Clone(Form("h_run1a_stops_hole_%.1fmm", r));
    TH2* h_run1b_muons_hole      = (TH2*) h_run1b_muons     ->Clone(Form("h_run1b_muons_hole_%.1fmm", r));
    TH2* h_run1b_calo_stops_hole = (TH2*) h_run1b_calo_stops->Clone(Form("h_run1b_calo_stops_hole_%.1fmm", r));

    for(int ix = 1; ix <= h_run1a_stops->GetNbinsX(); ++ix) {
      for(int iy = 1; iy <= h_run1a_stops->GetNbinsY(); ++iy) {
        double x = h_run1a_stops->GetXaxis()->GetBinCenter(ix) + 3904.; // Shift to solenoid center
        double y = h_run1a_stops->GetYaxis()->GetBinCenter(iy);
        double r_event = sqrt(x*x + y*y);
        if(r_event < r) { // Inside the hole
          h_run1b_muons_hole     ->SetBinContent(ix, iy, 0.); // not stopped in the plate
        } else { // outside the hole
          h_run1a_stops_hole     ->SetBinContent(ix, iy, 0.); // stopped in the plate
          h_run1b_calo_stops_hole->SetBinContent(ix, iy, 0.);
        }
      }
    }
    auto [run1a_metric_step, run1b_metric_step] = evaluate_metrics(h_run1a_stops_hole, h_run1b_muons_hole, h_run1b_calo_stops_hole);
    printf("Hole radius: %.1f mm, Run 1A Metric: %.3e, Run 1B Metric: %.3e\n", r, run1a_metric_step / run1a_metric, run1b_metric_step / run1b_metric);
    h_run1a_metric_vs_r->SetBinContent(istep + 1, run1a_metric_step / run1a_metric); // Relative to no hole case
    h_run1b_metric_vs_r->SetBinContent(istep + 1, run1b_metric_step / run1b_metric); // Relative to no hole case
    g_run1a_metric_vs_run1b_metric->AddPoint(run1a_metric_step / run1a_metric, run1b_metric_step / run1b_metric);
    delete h_run1a_stops_hole;
    delete h_run1b_muons_hole;
    delete h_run1b_calo_stops_hole;
  }

  //---------------------------------------------------------------------
  // Plot the metrics vs. hole radius
  //----------------------------------------------------------------------

  gStyle->SetOptStat(0);
  TCanvas c_metrics("c_metrics", "Metrics vs Hole Radius", 800, 600);
  h_run1a_metric_vs_r->SetLineColor(kBlue);
  h_run1b_metric_vs_r->SetLineColor(kRed);
  h_run1a_metric_vs_r->SetLineWidth(2);
  h_run1b_metric_vs_r->SetLineWidth(2);
  h_run1a_metric_vs_r->Draw("HIST");
  h_run1b_metric_vs_r->Draw("HIST SAME");
  h_run1a_metric_vs_r->SetTitle("Relative metrics vs. hole radius; Hole radius (mm);Relative metric");

  TLine line_50pc = TLine(r_min, 0.5, r_max, 0.5);
  line_50pc.SetLineStyle(kDashed);
  line_50pc.SetLineColor(kGray+1);
  line_50pc.SetLineWidth(2);
  line_50pc.Draw();

  const double max_val = max(h_run1a_metric_vs_r->GetMaximum(), h_run1b_metric_vs_r->GetMaximum());
  TLegend legend(0.6, 0.75, 0.89, 0.89);
  legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.AddEntry(h_run1a_metric_vs_r, "Run 1A Metric", "l");
  legend.AddEntry(h_run1b_metric_vs_r, "Run 1B Metric", "l");
  legend.Draw();

  h_run1a_metric_vs_r->GetYaxis()->SetRangeUser(0., 1.2*max_val);
  c_metrics.SaveAs(fig_dir_ + "/metrics_vs_hole_radius.png");
  h_run1a_metric_vs_r->GetYaxis()->SetRangeUser(max_val/1.e5, 20.*max_val);
  c_metrics.SetLogy();
  c_metrics.SaveAs(fig_dir_ + "/metrics_vs_hole_radius_log.png");

  // Plot the Run 1A metric vs. Run 1B metric
  TCanvas c_metric_correlation("c_metric_correlation", "Run 1A vs Run 1B Metric", 800, 600);
  g_run1a_metric_vs_run1b_metric->SetMarkerStyle(20);
  g_run1a_metric_vs_run1b_metric->SetMarkerSize(0.8);
  g_run1a_metric_vs_run1b_metric->SetMarkerColor(kBlack);
  g_run1a_metric_vs_run1b_metric->SetTitle("Run 1A vs Run 1B Metric;Run 1A Metric (relative);Run 1B Metric (relative)");
  g_run1a_metric_vs_run1b_metric->Draw("AP");
  c_metric_correlation.SaveAs(fig_dir_ + "/run1a_metric_vs_run1b_metric.png");

  return g_run1a_metric_vs_run1b_metric;
}

void test_random_geometries(TH2* h_run1a_stops, TH2* h_run1b_muons, TH2* h_run1b_calo_stops,
                           const double run1a_metric, const double run1b_metric,
                           TGraph* g_run1a_metric_vs_run1b_metric) {

  const int n_random = 1000;
  TGraph* g_random_metrics = new TGraph();
  TRandom3 rand(90);
  for(int i = 0; i < n_random; ++i) {
    const double p_kill = (i == 0) ? 0. : i == 1 ? 1. : rand.Uniform(); // Test no hole, full hole, and random intermediate cases
    TH2* h_run1a_stops_hole      = (TH2*) h_run1a_stops     ->Clone("h_run1a_stops_random");
    TH2* h_run1b_muons_hole      = (TH2*) h_run1b_muons     ->Clone("h_run1b_muons_random");
    TH2* h_run1b_calo_stops_hole = (TH2*) h_run1b_calo_stops->Clone("h_run1b_calo_stops_random");
    h_run1b_muons_hole->Reset(); // Start with no muons passing through the hole, then randomly add some back in

    for(int ix = 1; ix <= h_run1a_stops->GetNbinsX(); ++ix) {
      for(int iy = 1; iy <= h_run1a_stops->GetNbinsY(); ++iy) {
        const double p = rand.Uniform();
        if(p > p_kill) continue; // Keep this event
        h_run1a_stops_hole     ->SetBinContent(ix, iy, 0.); // stopped in the plate
        h_run1b_calo_stops_hole->SetBinContent(ix, iy, 0.);
        h_run1b_muons_hole     ->SetBinContent(ix, iy, h_run1b_muons->GetBinContent(ix, iy)); // not stopped in the plate
      }
    }

    auto [run1a_metric_step, run1b_metric_step] = evaluate_metrics(h_run1a_stops_hole, h_run1b_muons_hole, h_run1b_calo_stops_hole);
    g_random_metrics->AddPoint(run1a_metric_step / run1a_metric, run1b_metric_step / run1b_metric);
    delete h_run1a_stops_hole;
    delete h_run1b_muons_hole;
    delete h_run1b_calo_stops_hole;
  }

  //----------------------------------------------------------------------
  // Plot the random geometries
  //----------------------------------------------------------------------

  TCanvas c_random_metrics("c_random_metrics", "Random Geometries", 800, 600);
  g_random_metrics->SetMarkerStyle(20);
  g_random_metrics->SetMarkerSize(0.8);
  g_random_metrics->SetMarkerColor(kBlue);
  g_random_metrics->SetTitle("Random Geometries;Run 1A Metric (relative);Run 1B Metric (relative)");
  g_run1a_metric_vs_run1b_metric->Draw("AL");
  g_random_metrics->Draw("P SAME");
  g_run1a_metric_vs_run1b_metric->SetLineColor(kRed);
  g_run1a_metric_vs_run1b_metric->SetLineWidth(2);
  g_run1a_metric_vs_run1b_metric->Draw("L SAME");

  TLegend leg(0.6, 0.75, 0.89, 0.89);
  leg.SetBorderSize(0); leg.SetFillStyle(0);
  leg.AddEntry(g_random_metrics, "Random Geometries", "p");
  leg.AddEntry(g_run1a_metric_vs_run1b_metric, "Hole Radius Scan", "l");
  leg.Draw();
  c_random_metrics.SaveAs(fig_dir_ + "/random_geometries.png");
}

//--------------------------------------------------------------------------------------------------------------------
// Test random holes
void test_random_holes(TH2* h_run1a_stops, TH2* h_run1b_muons, TH2* h_run1b_calo_stops,
                       const double run1a_metric, const double run1b_metric, TGraph* g_run1a_metric_vs_run1b_metric) {
  TRandom rand(90);
  const int n_random = 50000;
  TGraph* g_random_metrics = new TGraph();
  const double r_min  =  10.;
  const double r_max  = 150.;
  const double dx_max = 100.;
  const double dy_max = 100.;
  for(int i = 0; i < n_random; ++i) {
    if(i % 1000 == 0) printf("Testing random hole %d / %d\n", i, n_random);
    const double r_in = rand.Uniform(r_min, r_max);
    const double r_out = rand.Uniform(r_in, r_max);
    const double dx = rand.Uniform(-dx_max, dx_max);
    const double dy = rand.Uniform(-dy_max, dy_max);
    const double phi_0 = rand.Uniform(0., 2.*M_PI);
    const double dphi  = rand.Uniform(0., 2.*M_PI);
    TH2* h_run1a_stops_hole      = (TH2*) h_run1a_stops     ->Clone(Form("h_run1a_stops_random_%d", i));
    TH2* h_run1b_muons_hole      = (TH2*) h_run1b_muons     ->Clone(Form("h_run1b_muons_random_%d", i));
    TH2* h_run1b_calo_stops_hole = (TH2*) h_run1b_calo_stops->Clone(Form("h_run1b_calo_stops_random_%d", i));

    for(int ix = 1; ix <= h_run1a_stops->GetNbinsX(); ++ix) {
      for(int iy = 1; iy <= h_run1a_stops->GetNbinsY(); ++iy) {
        const double x = h_run1a_stops->GetXaxis()->GetBinCenter(ix) + 3904. - dx; // Shift to solenoid center and apply random offset
        const double y = h_run1a_stops->GetYaxis()->GetBinCenter(iy) - dy;
        double phi = atan2(y, x) + phi_0; // Apply random rotation
        if(phi > 2.*M_PI) phi -= 2.*M_PI;
        if(phi < 0.) phi += 2.*M_PI;
        const double r_event = sqrt(x*x + y*y);
        if(r_event > r_in && r_event < r_out && phi < dphi) { // Inside the hole
          h_run1b_muons_hole     ->SetBinContent(ix, iy, 0.); // not stopped in the plate
        } else { // outside the hole
          h_run1a_stops_hole     ->SetBinContent(ix, iy, 0.); // stopped in the plate
          h_run1b_calo_stops_hole->SetBinContent(ix, iy, 0.);
        }
      }
    }
    auto [run1a_metric_step, run1b_metric_step] = evaluate_metrics(h_run1a_stops_hole, h_run1b_muons_hole, h_run1b_calo_stops_hole);
    g_random_metrics->AddPoint(run1a_metric_step / run1a_metric, run1b_metric_step / run1b_metric);
    delete h_run1a_stops_hole;
    delete h_run1b_muons_hole;
    delete h_run1b_calo_stops_hole;
  }

  //----------------------------------------------------------------------
  // Plot the random geometries
  //----------------------------------------------------------------------

  TCanvas c_random_metrics("c_random_metrics", "Random Geometries", 800, 600);
  g_random_metrics->SetMarkerStyle(20);
  g_random_metrics->SetMarkerSize(0.8);
  g_random_metrics->SetMarkerColor(kBlue);
  g_random_metrics->SetTitle("Random Geometries;Run 1A Metric (relative);Run 1B Metric (relative)");
  g_run1a_metric_vs_run1b_metric->Draw("AL");
  g_random_metrics->Draw("P SAME");
  g_run1a_metric_vs_run1b_metric->SetLineColor(kRed);
  g_run1a_metric_vs_run1b_metric->SetLineWidth(2);
  g_run1a_metric_vs_run1b_metric->Draw("L SAME");

  TLegend leg(0.6, 0.75, 0.89, 0.89);
  leg.SetBorderSize(0); leg.SetFillStyle(0);
  leg.AddEntry(g_random_metrics, "Random Holes", "p");
  leg.AddEntry(g_run1a_metric_vs_run1b_metric, "Hole Radius Scan", "l");
  leg.Draw();
  c_random_metrics.SaveAs(fig_dir_ + "/random_hole_geometries.png");
}

//--------------------------------------------------------------------------------------------------------------------
// Main function to optimize the hole shape
int optimize_hole() {

  //---------------------------------------------------------------------
  // Get the input data
  //----------------------------------------------------------------------

  TString base_dir = "/exp/mu2e/app/users/jmott/analysis/NoDSPlans/PlanB/run/CaloStopsTS5/Plots/";
  TFile* f_run1a_stops      = TFile::Open(base_dir + "FieldOn/outputHists_Run1A_TargetStops.root", "READ");
  TFile* f_run1b_muons      = TFile::Open(base_dir + "FieldOff/outputHists_Run1B_AllMuons.root"  , "READ");
  TFile* f_run1b_calo_stops = TFile::Open(base_dir + "FieldOff/outputHists_Run1B_CaloStops.root" , "READ");

  if(!f_run1a_stops || !f_run1b_muons || !f_run1b_calo_stops) {
    std::cerr << "Error opening files!" << std::endl;
    return 1;
  }

  TH2* h_run1a_stops      = (TH2*) f_run1a_stops     ->Get("hAftTS5Coll");
  TH2* h_run1b_muons      = (TH2*) f_run1b_muons     ->Get("hAftTS5Coll");
  TH2* h_run1b_calo_stops = (TH2*) f_run1b_calo_stops->Get("hAftTS5Coll");

  if(!h_run1a_stops || !h_run1b_muons || !h_run1b_calo_stops) {
    std::cerr << "Error retrieving histograms!" << std::endl;
    return 1;
  }

  h_run1a_stops      = (TH2*) h_run1a_stops     ->Clone("h_run1a_stops");
  h_run1b_muons      = (TH2*) h_run1b_muons     ->Clone("h_run1b_muons");
  h_run1b_calo_stops = (TH2*) h_run1b_calo_stops->Clone("h_run1b_calo_stops");
  h_run1a_stops      ->SetDirectory(0);
  h_run1b_muons      ->SetDirectory(0);
  h_run1b_calo_stops ->SetDirectory(0);
  f_run1a_stops     ->Close();
  f_run1b_muons     ->Close();
  f_run1b_calo_stops->Close();

  if(smooth_inputs_) {
    h_run1a_stops     ->Smooth(1);
    h_run1b_muons     ->Smooth(1);
    h_run1b_calo_stops->Smooth(1);
  }

  //----------------------------------------------------------------------
  // Make the output directory for the plots
  //----------------------------------------------------------------------

  fig_dir_ = "figures";
  gSystem->Exec(Form("mkdir -p %s", fig_dir_.Data()));

  //---------------------------------------------------------------------
  // Normalize the histograms given the sim number of protons on target
  //----------------------------------------------------------------------

  const double n_protons_run1a = 15.3e6*(5e7 / 639084.); // Simulated 15.3M events
  const double n_protons_run1b = 15.3e6*(5e7 / 639084.);

  h_run1a_stops      ->Scale(1.0 / n_protons_run1a);
  h_run1b_muons      ->Scale(1.0 / n_protons_run1b);
  h_run1b_calo_stops ->Scale(1.0 / n_protons_run1b);

  //---------------------------------------------------------------------
  // Effective total rates
  //----------------------------------------------------------------------

  const double rate_run1a_stops      = h_run1a_stops     ->Integral();
  const double rate_run1b_muons      = h_run1b_muons     ->Integral();
  const double rate_run1b_calo_stops = h_run1b_calo_stops->Integral();
  // optimal case with 100% Run 1B/Run 1A acceptance, no calo muon stops
  auto [run1a_metric, run1b_metric] = evaluate_metrics(h_run1a_stops, h_run1b_muons, nullptr);

  printf("-----------------------------------------------------------------\n");
  printf("Total rates per proton on target:\n");
  printf("Run 1A Target Stops:   %.3e\n", rate_run1a_stops);
  printf("Run 1B Muons:          %.3e\n", rate_run1b_muons);
  printf("Run 1B Calo Stops:     %.3e\n", rate_run1b_calo_stops);
  printf("Optimal Run 1A Metric: %.3e\n", run1a_metric);
  printf("Optimal Run 1B Metric: %.3e\n", run1b_metric);
  printf("-----------------------------------------------------------------\n");

  //---------------------------------------------------------------------
  // Evaluate the metrics vs. hole radius
  //----------------------------------------------------------------------

  TGraph* g_run1a_metric_vs_run1b_metric = optimize_hole_radius(h_run1a_stops, h_run1b_muons, h_run1b_calo_stops, run1a_metric, run1b_metric);

  //----------------------------------------------------------------------
  // Evaluate the metrics for each bin and plot
  //----------------------------------------------------------------------

  TGraph* g_metric_ordered_scan = metric_ordered_bin_scan(h_run1a_stops, h_run1b_muons, h_run1b_calo_stops, run1a_metric, run1b_metric, g_run1a_metric_vs_run1b_metric);

  //---------------------------------------------------------------------
  // Test random geometries
  //----------------------------------------------------------------------

  // test_random_geometries(h_run1a_stops, h_run1b_muons, h_run1b_calo_stops, run1a_metric, run1b_metric, g_run1a_metric_vs_run1b_metric);
  // test_random_holes     (h_run1a_stops, h_run1b_muons, h_run1b_calo_stops, run1a_metric, run1b_metric, g_run1a_metric_vs_run1b_metric);

  return 0;
}
