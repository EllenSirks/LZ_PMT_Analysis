#include <sstream>
#include <string>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>

void area_spec_fit(std::string folder){
    
    std::string OutputDirectory1 = "/home/ellen/Desktop/MPhys_Project/Output/";
    std::string PlotDirectory1 = "/home/ellen/Desktop/MPhys_Project/Output/";
    std::string PlotDirectory2 = "/home/ellen/Desktop/MPhys_Project/Output/";

    OutputDirectory1.append(folder);
    OutputDirectory1.append("_output/specs/peak_area.txt");

    std::ifstream infile1(OutputDirectory1);

    PlotDirectory1.append(folder);
    PlotDirectory1.append("_output/specs/");
    PlotDirectory1.append(folder);
    PlotDirectory1.append("-area_Spec.root");

    PlotDirectory2.append(folder);
    PlotDirectory2.append("_output/specs/");
    PlotDirectory2.append(folder);
    PlotDirectory2.append("-area_Spec.png");

    // std::ifstream infile1("/home/ellen/Desktop/MPhys_Project/Output/20180131-0002_output/peak_area.txt");

    std::vector<double> v1;

    std::string l1;
    while (std::getline(infile1, l1)) {
        std::istringstream iss(l1);
        double a;
        if (!(iss >> a)) { 
          break; 
        } 
        v1.push_back(a);
    }

    auto max = 0;
    for (const auto &val : v1) {
      if (val > max) {
        max = val;
      }
    }

    auto min = *std::min_element(v1.begin(), v1.end());

    TH1F * h1 = new TH1F("h1","h1", 300, min, 100);
    for (const auto &val : v1) {
      h1->Fill(val);
    }

    // TCanvas * c1 = new TCanvas("c1","c1");

    TH1F *h2 = (TH1F*) h1->Clone();

    Double_t par1[6];
    // Double_t par1[9];

    TF1 *g1    = new TF1("g1","gaus",min,5);
    // TF1 *g2    = new TF1("g2","gaus",5,46);
    TF1 *g2    = new TF1("g2","gaus",5,40);
    // TF1 *g3    = new TF1("g3","gaus",46,100);

    TF1 *total1 = new TF1("total1","gaus(0)+gaus(3)",min, 40);
    // TF1 *total1 = new TF1("total1","gaus(0)+gaus(3)+gaus(6)", min, 100);

    total1->SetLineColor(2);

    h1->Fit(g1,"R");
    h1->Fit(g2,"R+");
    // h1->Fit(g3,"qR+");

    g1->GetParameters(&par1[0]);
    g2->GetParameters(&par1[3]);
    // g3->GetParameters(&par1[6]);

    total1->SetParameters(par1);

    TCanvas * c2 = new TCanvas("c2","c2");

    std::string title = folder;
    title.append(" peak area spectrum");
    title.append(";Peak Area (mV*ns) ;Entries");
    h2->SetTitle(title.c_str());

    h2->Draw();
    h2->Fit(total1,"qR+");

    TF1 * f1 = h2->GetFunction("total1");
    // gStyle->SetOptFit(1);
    // gStyle->SetOptStat(kFALSE); 

    // std::cout << h2->GetMaximum() << std::endl;

    TLine *line2 = new TLine(g2->GetParameter(1), 0, g2->GetParameter(1), 2700);
    // TLine *line2 = new TLine(total1->GetParameter(4), 0, total1->GetParameter(4), 240);
    line2->SetLineColor(2);
    line2->SetLineWidth(2);
    line2->SetLineStyle(2);
    line2->Draw("same");

    std::cout << g2->GetParameter(0) << std::endl;
    std::cout << g2->GetParameter(1) << std::endl;
    std::cout << g2->GetParameter(2) << std::endl;

    c2->SetLogy();

    c2->SaveAs(PlotDirectory1.c_str());
    c2->SaveAs(PlotDirectory2.c_str());

    // delete c1;
    // delete c2;
    // delete h1;
    // delete h2;
}
