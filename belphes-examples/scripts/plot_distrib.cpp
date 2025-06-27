// plot_distrib.cpp 
/*
    Given a .root file containing multiple distributions, choose one of them, and 
    plot it in a .pdf file.
*/

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>

void plot_distrib(
    const std::string& input_file_dir, 
    const std::string& hist_name)
{
    // 1) Open the file
    TFile *file = TFile::Open(input_file_dir.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open file " << input_file_dir << std::endl;
        return;
    }

    // 2) Retrieve the histogram
    TH1 *hist = dynamic_cast<TH1*>(file->Get(hist_name.c_str()));
    if (!hist) {
        std::cerr << "Error: histogram \"" << hist_name 
                  << "\" not found in file " << input_file_dir << std::endl;
        file->Close();
        return;
    }

    // 3) Draw it
    TCanvas canvas("c1", "Distribution", 800, 600);
    hist->SetLineWidth(2);
    hist->Draw();

    // 4) Construct output PDF name
    std::string figure_folder = "../figs/";

    std::string base = input_file_dir;
    // strip any directory components
    size_t sep = base.find_last_of("/\\");
    if (sep != std::string::npos)
        base = base.substr(sep + 1);
    // strip the ".root" extension (case-sensitive)
    size_t ext = base.rfind(".root");
    if (ext != std::string::npos)
        base = base.substr(0, ext);
    
    std::string output = figure_folder + base + "_" + hist_name + ".pdf";
    std::cout << "File saved at " << output << std::endl;

    // 5) Save as PDF
    canvas.SaveAs(output.c_str());

    file->Close();
}

// Simple main() to invoke it from the command line:
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <hist_name>" << std::endl;
        return 1;
    }
    plot_distrib(argv[1], argv[2]);
    return 0;
}
