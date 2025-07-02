// PlotHistogram.C
/*
    Select one histogram in a given histogram and plot. Each histogram should follow
    the naming convention "hist_eta%d_pt%d". Also need to check the binning of PT 
    and abs(eta) before plotting.

    Compile with
    g++ -O3 PlotHistogram.cpp $(root-config --cflags --libs) -o ./PlotHistogram

    9 Apr 2025
*/

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>      // for std::stoi
#include <stdexcept>    // for std::invalid_argument, std::out_of_range

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>

using namespace std;

/// Removes any path components and a trailing ".root" extension (case-sensitive) from 'base'.
void stripPathAndRootExtension(string& base) {
    // 1) strip any directory components
    size_t sep = base.find_last_of("/\\");
    if (sep != string::npos) {
        base = base.substr(sep + 1);
    }

    // 2) strip the ".root" extension
    const string ext = ".root";
    if (base.size() >= ext.size() &&
        base.compare(base.size() - ext.size(), ext.size(), ext) == 0) {
        base.erase(base.size() - ext.size());
    }
}

void PlotHistogram(
    const string& input_file_dir, 
    int iEta, int iPt,
    const string& output_folder) 
{
    /* change this proportion if binning changes! */
    // Bins for pt and eta
    const vector<int> pt_bins = {20, 30, 40, 50, 60, 70, 80, 90, 100,
                                         120, 140, 160, 180, 200,
                                         250, 300, 400, 600, 1000};
    const vector<double> abs_eta_bins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    // Open the ROOT file
    TFile* file = TFile::Open(input_file_dir.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << input_file_dir << endl;
        return;
    }

    // Construct the histogram name based on iEta and iPt
    string histName = Form("hist_eta%d_pt%d", iEta, iPt);

    // Retrieve the histogram from the file
    TH1F* hist = dynamic_cast<TH1F*>(file->Get(histName.c_str()));
    if (!hist) {
        cerr << "Error: Histogram " << histName << " not found in file " << input_file_dir << endl;
        file->Close();
        return;
    }

    int totalEntries = hist->GetEntries();
    string totalEntries_str = to_string(totalEntries);
    cout << "The given histogram has " << totalEntries_str << " entries." << endl;

    if(totalEntries != 0) {
        hist->Scale(1.0/(totalEntries*1.0));
    }


    // Set the histogram title to include iEta, iPt, and the number of entries
    string title = Form("Prob. Distrib. for abs(Eta) in [%.1f, %.1f], Pt in [%d, %d] GeV;JetBtagDeepFlavB;",
        abs_eta_bins[iEta], abs_eta_bins[iEta+1],
        pt_bins[iPt], pt_bins[iPt+1]);
    hist->SetTitle(title.c_str());

    // Create a canvas and draw the histogram
    TCanvas* c = new TCanvas("c", "Histogram", 800, 600);
    // hist->SetStats(0);
    // hist->Draw("HIST L"); // include HIST so other style choices apply
    
    // also draw the CDF
    TH1F* cdf = dynamic_cast<TH1F*>(hist->Clone("cdf"));
    cdf->Reset();
    double cumulative = 0;
    int   nBins      = hist->GetNbinsX();
    for (int b = 1; b <= nBins; ++b) {
        cumulative += hist->GetBinContent(b);
        cdf->SetBinContent(b, cumulative);
    }

    cdf->SetStats(0);
    cdf->SetFillColor(kBlue-10);
    cdf->SetLineWidth(0);
    cdf->Draw("HIST");

    hist->Scale(3.0);
    hist->SetStats(0);
    hist->SetLineColor(kRed);
    hist->Draw("HIST L SAME"); // include HIST so other style choices apply

    // Legend 
    TLegend* leg = new TLegend(0.15, 0.70, 0.35, 0.85);
    leg->SetHeader((totalEntries_str+" entries").c_str(),"C");
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->AddEntry(cdf,  "CDF", "f");
    leg->AddEntry(hist, "PDF (x3)", "l");
    leg->Draw();


    // save the file in a designated directory 
    // get base file name 
    string base = input_file_dir;
    stripPathAndRootExtension(base);
    string file_name = output_folder + "/" + base + "_" + histName + ".pdf";
    c->Print(file_name.c_str());

    // Clean up
    file->Close();
    delete c;
}

int main(int argc, char *argv[]) {
    
    if (argc < 4 || argc > 5) {
        cerr << "Usage: " << argv[0]
             << " <input_root_file> <iEta> <iPt> [output_folder]\n";
        return 1;
    }

    string input_file = argv[1];
    int iEta, iPt;
    try {
        iEta = stoi(argv[2]);
        iPt  = stoi(argv[3]);
    }
    catch (const invalid_argument&) {
        cerr << "Error: iEta and iPt must be integers.\n";
        return 1;
    }
    catch (const out_of_range&) {
        cerr << "Error: iEta or iPt out of range.\n";
        return 1;
    }

    string output_folder = (argc == 5 ? argv[4] : ".");

    PlotHistogram(input_file, iEta, iPt, output_folder);

    return 0;
}