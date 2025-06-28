// PlotBTagScore.cpp 
/*
    Given a belphes file, plot the pseudo b-tag score of each jet.
*/
int PlotBTagScore(
    const std::string& input_file_dir,
    const std::string& output_folder="../figs/")
{
    // enable multi-threading
    ROOT::EnableImplicitMT();

    // create dataframe
    ROOT::RDataFrame df("Delphes", input_file_dir);

    TCanvas *c = new TCanvas(); // draw

    std::string col_name1 = "Jet.Jet_btagDeepFlavB";
    const char* col_name1_c_str = col_name1.c_str();

    int xbins   = 200;
    double xmin = -1.0;
    double xmax = 1.0;
    std::string title_string = "Jet_btagDeepFlavB Distribution;Jet_btagDeepFlavB score;Events";

    auto hist = df.Histo1D(
        {col_name1_c_str, title_string.c_str(), xbins, xmin, xmax}, 
        col_name1);

    hist->DrawClone();

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

    std::string file_name = output_folder + base + "_Jet_btagDeepFlavB_distrb.pdf";
    c->Print(file_name.c_str());

    return 0;
}