/*
    Plot the number of jets in the given file. Would expect at least 2 per event
    plus initial / final state radiations for "p p > z > b b~".

    To run, do 
    ```
    root -l -q 'PlotJetSize.cpp'
    ```

    25 June 2025
*/
int PlotJetSize() 
{
    // enable multi-threading
    ROOT::EnableImplicitMT();

    // create dataframe
    const std::string input_file_dir = "pp_z_bb.root";
    ROOT::RDataFrame df("Delphes", input_file_dir);


/******************************************************************************/


    TCanvas *c = new TCanvas(); // draw

    std::string col_name1 = "Jet_size";
    const char* col_name1_c_str = col_name1.c_str();

    int xbins   = 6;
    double xmin = 0.0;
    double xmax = 6.0;
    std::string title_string = "Number of Jets per Event in p p > z > b b~;Jets;Events";

    auto hist = df.Histo1D(
        {col_name1_c_str, title_string.c_str(), xbins, xmin, xmax}, 
        col_name1);

    hist->DrawClone();

    std::string file_name = "number_of_jets.pdf";
    c->Print(file_name.c_str());

    return 0;
}