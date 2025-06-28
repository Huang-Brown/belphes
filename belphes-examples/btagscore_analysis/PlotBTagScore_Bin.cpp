// PlotBTagScore_Bin.cpp
/*
    Given a belphes file, plot the pseudo b-tag score of all jets that fall into
    a certain interval.
*/
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TString.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

ROOT::RVec<Float_t> GetBTagDeepFlavB(
    const ROOT::RVec<Float_t>&   Jet_Eta,
    const ROOT::RVec<Float_t>&   Jet_PT,
    const std::vector<Float_t>&  abs_eta_interval,
    const std::vector<Float_t>&  pt_interval,
    const ROOT::RVec<Float_t>&   Jet_btagDeepFlavB)
{
    ROOT::RVec<Float_t> ret;
    for (size_t i = 0; i < Jet_Eta.size(); ++i) {
        Float_t aeta = std::abs(Jet_Eta[i]);
        Float_t pt   = Jet_PT[i];
        if (aeta >= abs_eta_interval[0] && aeta < abs_eta_interval[1] &&
            pt   >= pt_interval[0]      && pt   < pt_interval[1]) {
            ret.emplace_back(Jet_btagDeepFlavB[i]);
        }
    }
    return ret;
}

// A to-string method for pt and eta intervals 
std::string interval_to_str(const std::vector<Float_t>& my_interval) {
    if (my_interval.size() < 2) {
        return "{}";
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << my_interval[0];
    std::string str0 = oss.str();

    // clear the buffer and any flags
    oss.str("");
    oss.clear();

    oss << std::fixed << std::setprecision(1) << my_interval[1];
    std::string str1 = oss.str();

    return "{" + str0 + ", " + str1 + "}";
}

int PlotBTagScore_Bin(
    const std::string&           input_file_dir,
    const std::vector<Float_t>&  pt_interval,
    const std::vector<Float_t>&  abs_eta_interval,
    const std::string&           output_folder="../figs/")
{
    // 1) ensure output folder exists
    gSystem->mkdir(output_folder.c_str(), /*recursive=*/true);

    // 2) load the TTree “Delphes” from your file
    ROOT::RDataFrame df("Delphes", input_file_dir);

    // 3) define a new column “selBTags” = the vector of scores in the bin
    auto df2 = df.Define("selBTags",
        [&](const ROOT::RVec<Float_t>& eta,
            const ROOT::RVec<Float_t>& pt,
            const ROOT::RVec<Float_t>& btag)
        {
            return GetBTagDeepFlavB(eta, pt, abs_eta_interval, pt_interval, btag);
        },
        {"Jet.Eta", "Jet.PT", "Jet.Jet_btagDeepFlavB"})
      .Filter("selBTags.size()>0");  // keep only events with ≥1 jet in bin

    // 4) book a histogram of all selected scores
    //    here: 100 bins from 0 to 1 — adjust as needed
    std::string title_str = "Jet_btagDeepFlavB on PT=" + interval_to_str(pt_interval) + ", abs(Eta)=" + interval_to_str(abs_eta_interval) + ";Jet_btagDeepFlavB;Count";
    auto h = df2.Histo1D(
        {"hBtag", title_str.c_str(), 100, 0.0, 1.0},
        "selBTags");

    // 5) draw & save
    TCanvas c;
    h->Draw();
    h->SetStats(0);

    // get base file name 
    std::string base = input_file_dir;
    // strip any directory components
    size_t sep = base.find_last_of("/\\");
    if (sep != std::string::npos)
        base = base.substr(sep + 1);
    // strip the ".root" extension (case-sensitive)
    size_t ext = base.rfind(".root");
    if (ext != std::string::npos)
        base = base.substr(0, ext);

    std::string file_name = 
        output_folder + base + "_Jet_btagDeepFlavB_distrb" 
        + "_eta_" + interval_to_str(abs_eta_interval) 
        + "_pt_" + interval_to_str(pt_interval) + ".pdf";
    c.Print(file_name.c_str());

    return 0;
}