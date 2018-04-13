#include <sstream>
#include <string>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>
#include <math.h>

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
    PlotDirectory1.append("-area_Spec_fit.root");

    PlotDirectory2.append(folder);
    PlotDirectory2.append("_output/specs/");
    PlotDirectory2.append(folder);
    PlotDirectory2.append("-area_Spec_fit.png");

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

    TH1F * h1 = new TH1F("h1","h1", 300, min, 250);
    for (const auto &val : v1) {
      h1->Fill(val);
    }

    TCanvas * c1 = new TCanvas("c1","c1");

    TH1F *h2 = (TH1F*) h1->Clone();

    Double_t par1[6];
    // Double_t par1[9];

    TF1 *g1    = new TF1("g1","gaus", min, 10);
    TF1 *g2    = new TF1("g2","gaus", 40, 80);

    h1->Fit(g1,"R");
    h1->Fit(g2,"R+");

    // TF1 *g3    = new TF1("g3","gaus", 20,100);
    // g3->SetParameter(1, g2->GetParameter(1));

    // h1->Fit(g3,"qR+");

    g1->GetParameters(&par1[0]);
    g2->GetParameters(&par1[3]);
    // g3->GetParameters(&par1[6]);

    TF1 *total1 = new TF1("total1","gaus(0)+gaus(3)", -5, 200);
    // TF1 *total1 = new TF1("total1","gaus(0)+gaus(3)+gaus(6)", -5, 100);
    total1->SetParameters(par1);
    total1->SetLineColor(2);

    // c1->SetLogy();

    TCanvas * c2 = new TCanvas("c2","c2");

    // std::string title = folder;
    std::string title;
    // title.append(" peak area spectrum");
    title.append(";Peak Area (mV*ns) ;Entries");
    h2->SetTitle(title.c_str());

    h2->Draw();
    h2->Fit(total1,"qR+");
    TF1 * f1 = h2->GetFunction("total1");
    gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);



	Double_t chi2 = g2->GetChisquare()/g2->GetNDF();
	std::cout << "chi2:" << chi2 << std::endl;
    chi2 = total1->GetChisquare()/total1->GetNDF();
    std::cout << "chi2:" << chi2 << std::endl;

    std::cout << "mean: " << g2->GetParameter(1) << "+/-" << g2->GetParError(1)<< std::endl;
    std::cout << "sigma: " << g2->GetParameter(2) << "+/-" << g2->GetParError(2)<< std::endl;

    auto legend = new TLegend(0.45, 0.75, 0.9, 0.9);

    std::stringstream ss1, ss2, ss3, ss4;
    ss1 << roundf(g2->GetParameter(1) * 100) / 100;
    ss2 << roundf(g2->GetParError(1) * 100) / 100;
    ss3 << roundf(g2->GetParameter(2) * 100) / 100;
    ss4 << roundf(g2->GetParError(2) * 100) / 100;

    // auto m = 6.95433;
    // auto e = 0.252195;

    // ss1 << roundf(m * 100) / 100;
    // ss2 << roundf(e * 100) / 100;

    std::string myString1 =  "Sphe Peak Mean: " + ss1.str() + "(" + ss2.str() + ")";
    std::string myString2 =  "Sphe Peak Width: " + ss3.str() + "(" + ss4.str() + ")";



    TLine *line = new TLine(g2->GetParameter(1), 0, g2->GetParameter(1), 320);
    // TLine *line = new TLine(m, 0, m, 1700);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("same");


    // std::cout << g2->Eval(g2->GetParameter(1) - g2->GetParameter(2)) << std::endl;
    TLine *line2 = new TLine(g2->GetParameter(1) - g2->GetParameter(2), g2->Eval(g2->GetParameter(1) - g2->GetParameter(2)), g2->GetParameter(1) + g2->GetParameter(2), g2->Eval(g2->GetParameter(1) + g2->GetParameter(2)));
    line2->SetLineWidth(3);
    line2->SetLineStyle(2);
    line2->Draw("same");

    TLine *line3 = new TLine(g2->GetParameter(1) - g2->GetParameter(2), 0, g2->GetParameter(1) - g2->GetParameter(2), g2->Eval(g2->GetParameter(1) + g2->GetParameter(2)));
    line3->SetLineWidth(2);
    // line3->SetLineStyle(2);
    // line3->Draw("same");   


    legend->AddEntry(line, myString1.c_str(), "l");
    legend->AddEntry(line2, myString2.c_str(), "l");

    legend->SetTextSize(0.0375);
    legend->Draw();




    c2->SetLogy();
    c2->SaveAs(PlotDirectory1.c_str());
    c2->SaveAs(PlotDirectory2.c_str());

    // delete c1;
    // delete c2;
    // delete h1;
    // delete h2;

    // TCanvas * c3 = new TCanvas("c3","c3");
    // total1->Draw();



}
