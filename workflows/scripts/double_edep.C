// Evaluate the double edep distribution

//----------------------------------------------------------------------------------------------------------
// Double a histogram by sampling from the original distribution
TH1* mix_histograms(const TH1* h_1, TH1* h_2, TRandom3& rnd) {
    if(!h_1 || !h_2) {
        std::cerr << "Error: Null histogram pointer passed to mix_histograms" << std::endl;
        return nullptr;
    }

    // Create a new histogram with the same binning as the input histograms
    TH1* h_mixed = dynamic_cast<TH1*>(h_1->Clone("mixed_histogram"));
    h_mixed->Reset();

    // Sample each histogram and fill the mixed histogram
    const double bin_width = h_1->GetBinWidth(1); // Assuming uniform binning
    const long max_entries = 1e6;
    for(long i = 0; i < max_entries; ++i) {
        const double x_1 = h_1->GetRandom(&rnd);
        const double x_2 = h_2->GetRandom(&rnd);
        const double doubled = (x_1 + x_2);
        h_mixed->Fill(doubled);
    }
    h_mixed->Scale(1. / max_entries); // Normalize the histogram

    return h_mixed;
}

//----------------------------------------------------------------------------------------------------------
// Main function
int double_edep(const char* file_name,
                const double ngen = 100000.,
                const double n_per_event = 1000.) {

    // Open the histogram file
    TFile* file = TFile::Open(file_name, "READ");
    if(!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return 1;
    }

    // Get the total calo energy histogram
    TH1* h_edep = dynamic_cast<TH1*>(file->Get("EDepAna/hist_1/total_calo_energy"));
    if(!h_edep) {
        std::cerr << "Error: Histogram not found in file: " << file_name << std::endl;
        file->Close();
        return 1;
    }

    const double eff = h_edep->GetEntries() / ngen;
    const double n_exp = eff * n_per_event;
    h_edep->Scale(1./h_edep->Integral()); // Normalize to probability distribution

    // FIXME: Guess of p(time overlap) and p(space overlap)
    const double p_time_overlap = 0.05;
    const double p_space_overlap = 50.*50. / (3.14*(700.*700. - 350.*350.)); // area of 50x50 cm^2 over area of annulus between 350 and 700 cm radius
    const double p_overlap = p_time_overlap * p_space_overlap;


    // Efficiencies for double and triple edep assuming independent events
    const double n_exp_double = n_exp * p_overlap;
    const double n_exp_triple = n_exp * p_overlap * n_exp_double;
    std::cout << "Efficiency for single edep: " << eff << std::endl;
    std::cout << "Expected per event for single edep: " << n_exp << std::endl;
    std::cout << "Expected per event for double edep: " << n_exp_double << std::endl;
    std::cout << "Expected per event for triple edep: " << n_exp_triple << std::endl;

    // Get the double and triple edep histograms
    TRandom3 rnd(90);
    TH1* h_double = mix_histograms(h_edep  , h_edep, rnd);
    TH1* h_triple = mix_histograms(h_double, h_edep, rnd);

    // Plot the result
    gSystem->Exec("mkdir -p figures");
    gStyle->SetOptStat(0);
    TCanvas c("c", "Double Edep", 800, 600);
    h_edep->SetLineColor(kBlue);
    h_double->SetLineColor(kRed);
    h_triple->SetLineColor(kGreen);
    h_edep->SetLineWidth(2);
    h_double->SetLineWidth(2);
    h_triple->SetLineWidth(2);
    h_edep->SetTitle("Total Calo Energy;Energy (MeV);Probability");
    h_edep->Draw("HIST");
    h_double->Draw("HIST SAME");
    h_triple->Draw("HIST SAME");
    TLegend legend(0.6, 0.7, 0.9, 0.9);
    legend.AddEntry(h_edep, "Single Edep", "l");
    legend.AddEntry(h_double, "Double Edep", "l");
    legend.AddEntry(h_triple, "Triple Edep", "l");
    legend.Draw();
    c.SaveAs("figures/double_edep.png");
    h_edep->GetYaxis()->SetRangeUser(1e-6, 5);
    c.SetLogy();
    c.SaveAs("figures/double_edep_log.png");

    file->Close();
    return 0;
}