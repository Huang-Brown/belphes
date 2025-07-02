// CompareBTagScore.cpp
/**
 * CompareBTagScore:
 *   Given a Belphes .root file, select all b-jets in a specified (pt, eta) bin,
 *   plot their pseudo b-tag scores, overlay the sampling PDF, compute χ²,
 *   and save the result as a PDF.
 *
 * Usage:
 *   ./PlotBTagScore_Bin <input_belphes_root> <input_pdf_root>
 *                     <abs_eta_bin_index> <pt_bin_index> [output_folder]
 *
 * Example:
 *   ./PlotBTagScore_Bin belphes.root btag_pdf.root 2 5 ./figs
 * 
 * Compile with
 *   g++ -O3 CompareBTagScore.cpp $(root-config --cflags --libs) -o ./CompareBTagScore
 *
 * 29 June 2025
 */

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

using namespace std;
using namespace ROOT;

//--------------------------------------------------------------------------------
/**
 * Select pseudo b-tag scores for b-flavor jets within given (η, pt) intervals.
 */
RVec<Float_t> GetBTagDeepFlavB(
    const RVec<Float_t>& Jet_Eta,
    const RVec<Float_t>& Jet_PT,
    const vector<Float_t>& abs_eta_interval,
    const vector<Float_t>& pt_interval,
    const RVec<Float_t>& Jet_btagDeepFlavB,
    const RVec<UInt_t>&  Jet_Flavor)
{
    RVec<Float_t> selected;
    for (size_t i = 0; i < Jet_Eta.size(); ++i) {
        Float_t aeta = fabs(Jet_Eta[i]);
        Float_t pt   = Jet_PT[i];
        UInt_t  flav = Jet_Flavor[i];
        if (flav == 5 &&
            aeta >= abs_eta_interval[0] && aeta < abs_eta_interval[1] &&
            pt   >= pt_interval[0]      && pt   < pt_interval[1])
        {
            selected.push_back(Jet_btagDeepFlavB[i]);
        }
    }
    return selected;
}

//--------------------------------------------------------------------------------
/**
 * Convert a two-element interval to a string {min, max} with one decimal.
 */
string interval_to_str(const vector<Float_t>& interval) {
    if (interval.size() < 2) return "{}";
    ostringstream oss;
    oss << fixed << setprecision(1) << interval[0];
    string lo = oss.str();
    oss.str(""); oss.clear();
    oss << fixed << setprecision(1) << interval[1];
    string hi = oss.str();
    return "{" + lo + ", " + hi + "}";
}

//--------------------------------------------------------------------------------
/**
 * Core routine: build and draw histograms, overlay input PDF, compute χ²,
 * and save canvas to output_folder.
 */
int CompareBTagScore(
    const string& input_file,
    const string& pdf_file,
    int           abs_eta_bin,
    int           pt_bin,
    const string& output_folder)
{
    // Define bin edges
    const vector<int>   pt_bins       = {20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,400,600,1000};
    const vector<double> abs_eta_bins = {0.0,0.5,1.0,1.5,2.0,2.5};

    // Determine intervals from bin indices
    vector<Float_t> pt_interval      = {Float_t(pt_bins[pt_bin]),   Float_t(pt_bins[pt_bin+1])};
    vector<Float_t> abs_eta_interval = {Float_t(abs_eta_bins[abs_eta_bin]), Float_t(abs_eta_bins[abs_eta_bin+1])};

    // Create output folder if needed
    gSystem->mkdir(output_folder.c_str(), /*recursive=*/true);

    // Load Delphes tree and filter
    RDataFrame df("Delphes", input_file);
    auto df_sel = df.Define(
            "selBTags",
            [&](const RVec<Float_t>& eta,
                const RVec<Float_t>& pt,
                const RVec<Float_t>& b,
                const RVec<UInt_t>&  f)
            { return GetBTagDeepFlavB(eta,pt,abs_eta_interval,pt_interval,b,f); },
            {"Jet.Eta","Jet.PT","Jet.Jet_btagDeepFlavB","Jet.Flavor"}
        )
        .Filter("selBTags.size()>0");

    // Histogram of selected scores
    string title = "b-tag score on pt=" + interval_to_str(pt_interval)
                 + ", abs(eta)=" + interval_to_str(abs_eta_interval)
                 + ";score;count";
    auto h_sel = df_sel.Histo1D({"hSel", title.c_str(), 100, 0.0, 1.0}, "selBTags");
    double nSel = h_sel->Integral();
    TH1* hPtr = h_sel.GetPtr();

    // Draw selected-data histogram
    TCanvas c("c", "bTagScore", 800,600);
    h_sel->SetStats(0);
    h_sel->SetFillColor(kBlue-10);
    h_sel->Draw("HIST");

    // Overlay input PDF histogram
    TFile* f_pdf = TFile::Open(pdf_file.c_str(), "READ");
    if (!f_pdf || f_pdf->IsZombie()) {
        cerr << "Error: cannot open PDF file " << pdf_file << endl;
        return 1;
    }
    string hist_name = Form("hist_eta%d_pt%d", abs_eta_bin, pt_bin);
    TH1F* h_pdf = dynamic_cast<TH1F*>(f_pdf->Get(hist_name.c_str()));
    if (!h_pdf) {
        cerr << "Error: PDF histogram " << hist_name << " not found in " << pdf_file << endl;
        f_pdf->Close();
        return 1;
    }
    // Normalize both to unity
    h_sel->Scale(1.0/nSel);
    h_pdf->Scale(1.0/h_pdf->GetEntries());

    // Compute chi-square
    double chi2 = 0;
    int bins = hPtr->GetNbinsX();
    for (int b = 1; b <= bins; ++b) {
        double O = hPtr->GetBinContent(b);
        double E = h_pdf->GetBinContent(b);
        if (E > 0) {
            double d = O - E;
            chi2 += (d*d) / E;
        }
    }
    cout << "Normalized chi2 = " << chi2 << endl;

    // Restore scale for plotting
    h_sel->Scale(nSel);
    h_pdf->Scale(nSel);

    h_pdf->SetStats(0);
    h_pdf->SetLineColor(kRed);
    h_pdf->SetLineWidth(2);
    h_pdf->Draw("HIST SAME");

    // Legend with chi2
    TLegend leg(0.15,0.70,0.50,0.85);
    leg.SetHeader(("#chi^{2}=" + to_string(chi2)).c_str(),"C");
    leg.SetBorderSize(1);
    leg.SetFillStyle(0);
    leg.AddEntry(hPtr, "Selected b-jets","f");
    leg.AddEntry(h_pdf, "Input PDF","l");
    leg.Draw();

    // Save canvas
    string base = input_file;
    if (auto p = base.find_last_of("/\\"); p!=string::npos) base = base.substr(p+1);
    if (auto p = base.rfind(".root"); p!=string::npos) base = base.substr(0,p);
    string out = output_folder + "/" + base +
                 "_eta" + interval_to_str(abs_eta_interval) +
                 "_pt"  + interval_to_str(pt_interval) + ".pdf";
    c.SaveAs(out.c_str());

    f_pdf->Close();
    return 0;
}

//--------------------------------------------------------------------------------
/**
 * Entry point: parse arguments and invoke CompareBTagScore().
 */
int main(int argc, char** argv) {
    if (argc < 5 || argc > 6) {
        cout << "Usage: " << argv[0]
             << " <input_belphes_root> <input_pdf_root>"
             << " <abs_eta_bin> <pt_bin> [output_folder]" << endl;
        return 1;
    }
    string input_file   = argv[1];
    string pdf_file     = argv[2];
    int    abs_eta_bin  = stoi(argv[3]);
    int    pt_bin       = stoi(argv[4]);
    string output_folder = (argc == 6 ? argv[5] : string("./figs"));

    return CompareBTagScore(input_file, pdf_file, abs_eta_bin, pt_bin, output_folder);
}


// // PlotBTagScore_Bin.cpp
// /*
//     Given a belphes file (.root), plot the pseudo b-tag score of all b jets that 
//     fall intoa certain interval. Overlay the plot with the distribution from which
//     the interval is sampled.

//     int PlotBTagScore_Bin(
//     const std::string&           input_file_dir,
//     const std::vector<Float_t>&  pt_interval,
//     const std::vector<Float_t>&  abs_eta_interval,
//     const std::string&           output_folder="../figs/")

//     29 June 2025
// */
// #include <ROOT/RDataFrame.hxx>
// #include <TCanvas.h>
// #include <TH1F.h>
// #include <TSystem.h>
// #include <TString.h>
// #include <sstream>
// #include <iomanip>
// #include <vector>
// #include <string>
// #include <cstdlib>

// using namespace std;

// ROOT::RVec<Float_t> GetBTagDeepFlavB(
//     const ROOT::RVec<Float_t>&   Jet_Eta,
//     const ROOT::RVec<Float_t>&   Jet_PT,
//     const std::vector<Float_t>&  abs_eta_interval,
//     const std::vector<Float_t>&  pt_interval,
//     const ROOT::RVec<Float_t>&   Jet_btagDeepFlavB,
//     const ROOT::RVec<UInt_t>&    Jet_Flavor)
// {
//     ROOT::RVec<Float_t> ret;
//     for (size_t i = 0; i < Jet_Eta.size(); ++i) {
//         Float_t aeta = std::abs(Jet_Eta[i]);
//         Float_t pt   = Jet_PT[i];
//         UInt_t  flav = Jet_Flavor[i];
//         if (aeta >= abs_eta_interval[0] && aeta < abs_eta_interval[1] &&
//             pt   >= pt_interval[0]      && pt   < pt_interval[1]      &&
//             flav == 5) {
//             ret.emplace_back(Jet_btagDeepFlavB[i]);
//         }
//     }
//     return ret;
// }

// // A to-string method for pt and eta intervals 
// std::string interval_to_str(const std::vector<Float_t>& my_interval) {
//     if (my_interval.size() < 2) {
//         return "{}";
//     }

//     std::ostringstream oss;
//     oss << std::fixed << std::setprecision(1) << my_interval[0];
//     std::string str0 = oss.str();

//     // clear the buffer and any flags
//     oss.str("");
//     oss.clear();

//     oss << std::fixed << std::setprecision(1) << my_interval[1];
//     std::string str1 = oss.str();

//     return "{" + str0 + ", " + str1 + "}";
// }

// int CompareBTagScore(
//     const std::string&           input_file_dir,
//     const std::string&           B_distribution_dir,
//     int                          abs_eta_bin,
//     int                          pt_bin,
//     const std::string&           output_folder="../figs/")
// {
//     // 0) prep work 
//     /* change this proportion if binning changes! */
//     // Bins for pt and eta
//     const vector<int> pt_bins = {20, 30, 40, 50, 60, 70, 80, 90, 100,
//                                          120, 140, 160, 180, 200,
//                                          250, 300, 400, 600, 1000};
//     const vector<double> abs_eta_bins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

//     std::vector<Float_t> pt_interval      = {static_cast<float>(pt_bins[pt_bin]), static_cast<float>(pt_bins[pt_bin+1])};
//     std::vector<Float_t> abs_eta_interval = {(Float_t)abs_eta_bins[abs_eta_bin], (Float_t)abs_eta_bins[abs_eta_bin+1]};

//     // 1) ensure output folder exists
//     gSystem->mkdir(output_folder.c_str(), /*recursive=*/true);

//     // 2) load the TTree “Delphes” from your file
//     ROOT::RDataFrame df("Delphes", input_file_dir);

//     // 3) define a new column “selBTags” = the vector of scores in the bin
//     auto df2 = df.Define("selBTags",
//         [&](const ROOT::RVec<Float_t>& eta,
//             const ROOT::RVec<Float_t>& pt,
//             const ROOT::RVec<Float_t>& btag,
//             const ROOT::RVec<UInt_t>&  flavor)
//         {
//             return GetBTagDeepFlavB(eta, pt, abs_eta_interval, pt_interval, btag, flavor);
//         },
//         {"Jet.Eta", "Jet.PT", "Jet.Jet_btagDeepFlavB", "Jet.Flavor"})
//       .Filter("selBTags.size()>0");  // keep only events with ≥1 jet in bin

//     // 4) book a histogram of all selected scores
//     //    here: 100 bins from 0 to 1 — adjust as needed
//     std::string title_str = "Jet_btagDeepFlavB on PT=" + interval_to_str(pt_interval) + ", abs(Eta)=" + interval_to_str(abs_eta_interval) + ";Jet_btagDeepFlavB;Count";
//     auto h = df2.Histo1D(
//         {"hBtag", title_str.c_str(), 100, 0.0, 1.0},
//         "selBTags");
    
//     // get raw TH1* from the RResultPtr for legend
//     TH1* hPtr = h.GetPtr();
    
//     double h_count = h->Integral();

//     // 5) draw 
//     TCanvas* c = new TCanvas("c", "Histogram", 800, 600);

//     h->SetStats(0);
//     h->SetFillColor(kBlue-10);
//     h->Draw("HIST");


//     // -------------------------------------------------------------------------
//     // 6) overlay with the corresponding histogram
//     // Open the ROOT file
//     TFile* B_file = TFile::Open(B_distribution_dir.c_str(), "READ");
//     if (!B_file || B_file->IsZombie()) {
//         cerr << "Error: Cannot open file " << B_distribution_dir << endl;
//         return -1;
//     }

//     // Construct the histogram name based on iEta and iPt
//     string histName = Form("hist_eta%d_pt%d", abs_eta_bin, pt_bin);

//     // Retrieve the histogram from the file
//     TH1F* B_hist = dynamic_cast<TH1F*>(B_file->Get(histName.c_str()));
//     if (!B_hist) {
//         cerr << "Error: Histogram " << histName << " not found in file " << B_distribution_dir << endl;
//         B_file->Close();
//         return -1;
//     }

//     int totalEntries = B_hist->GetEntries();
//     if(totalEntries != 0) {
//         B_hist->Scale(1.0/totalEntries); // scale according to the scale of h
//     } else {
//         cerr << "The requested distribution has 0 events!" << endl;
//     }

//     // 7) compute chi-sq
//     double chi2 = 0;
//     int    nbins = h->GetNbinsX();
//     h->Scale(1.0/h_count);

//     // make sure both histograms are defined over the same binning!
//     // the normalized Chi-sq is calculated by normalizing the observed and expected
//     // datasets first, then summing up the chi-sq per bin.
//     for (int bin = 1; bin <= nbins; ++bin) {
//         double O = hPtr->GetBinContent(bin);    // observed (your selected-data)
//         double E = B_hist->GetBinContent(bin);  // expected (input PDF, scaled)
//         if (E > 0) {
//             double d = O - E;
//             chi2 += (d*d) / E;
//         }
//     }
//     cout << "(Normalized) chi-sq: " << chi2 << endl;

//     // scale h and B_hist back 
//     h->Scale(1.0*h_count);
//     B_hist->Scale(1.0*h_count);

//     B_hist->SetStats(0);
//     B_hist->SetLineColor(kRed);
//     B_hist->SetLineWidth(3);
//     B_hist->Draw("HIST SAME L");


//     // Legend 
//     TLegend* leg = new TLegend(0.15, 0.70, 0.50, 0.85);
//     leg->SetHeader(("Normalized #chi^2 = " + std::to_string(chi2)).c_str(),"C");
//     leg->SetBorderSize(1);
//     leg->SetFillColor(0);
//     leg->AddEntry(h.GetPtr(),  "Belphes .root file", "f");
//     leg->AddEntry(B_hist, "Belphes input PDF, scaled", "l");
//     leg->Draw();


//     // get base file name 
//     std::string base = input_file_dir;
//     // strip any directory components
//     size_t sep = base.find_last_of("/\\");
//     if (sep != std::string::npos)
//         base = base.substr(sep + 1);
//     // strip the ".root" extension (case-sensitive)
//     size_t ext = base.rfind(".root");
//     if (ext != std::string::npos)
//         base = base.substr(0, ext);

//     std::string file_name = 
//         output_folder + "/" + base + "_Jet_btagDeepFlavB_distrb" 
//         + "_eta_" + interval_to_str(abs_eta_interval) 
//         + "_pt_" + interval_to_str(pt_interval) + ".pdf";
//     c->Print(file_name.c_str());

//     return 0;
// }