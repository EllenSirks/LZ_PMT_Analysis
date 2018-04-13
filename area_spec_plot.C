 #include <sstream>
#include <string>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>

void area_spec_plot(std::string folder){
    
    std::string OutputDirectory = "/home/ellen/Desktop/MPhys_Project/Output/" + folder + "_output/specs/peak_area.txt";
    std::string PlotDirectory1 = "/home/ellen/Desktop/MPhys_Project/Output/" + folder + "_output/specs/" + folder + "-area_Spec.root";
    std::string PlotDirectory2 = "/home/ellen/Desktop/MPhys_Project/Output/" + folder + "_output/specs/" + folder + "-area_Spec.png";

    std::ifstream infile(OutputDirectory);
    std::vector<double> v;
    std::string l;

    while (std::getline(infile, l)) {
        std::istringstream iss(l);
        double a;
        if (!(iss >> a)) { 
          break; 
        } 
        v.push_back(a);
    }

	if (v.size() == 0) {
		std::cout << "Input file is empty" << std::endl;
		return;
	}

    auto min = *std::min_element(v.begin(), v.end());

    TH1F * h = new TH1F("h","h", 300, min, 250);
    for (const auto &val : v) {
      h->Fill(val);
    }

    TCanvas * c = new TCanvas("c","c");

    // std::string title = folder;
    std::string title;
    // title.append(" peak area spectrum");
    title.append(";Peak Area (mV*ns) ;Entries");
    h->SetTitle(title.c_str());
    h->Draw();

    gStyle->SetOptStat(11);

    c->SetLogy();
    c->SaveAs(PlotDirectory1.c_str());
    c->SaveAs(PlotDirectory2.c_str());

    // delete c;
    // delete h;
}
