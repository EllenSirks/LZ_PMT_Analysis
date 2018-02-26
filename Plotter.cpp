#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <stdio.h>

/*** Method to Plot Single Waveforms ***/

int Plotter(std::string fileno, std::string folder) {

	gROOT->SetBatch(kTRUE);
	gErrorIgnoreLevel = kWarning;
	gPrintViaErrorHandler = kTRUE;

	std::string OutputDirectory = "/home/ellen/Desktop/MPhys_Project/Output/Plots/";
	std::string DataDirectory = "/media/Ellen/Elements/PMT_Data/" + folder + "/";
	std::string DataFileName = folder + "_";
	std::string filetype_data = ".csv", filetype_plot = ".png";

	std::string x_unit;
	std::string filename = DataDirectory + DataFileName + fileno + filetype_data;

	std::ifstream myfile(filename);
	std::vector<double> C1, C2, C3;
	int line_no = 0;

	for (std::string col1, col2; std::getline(myfile, col1, ',') && std::getline(myfile, col2);) {
		if (line_no == 1) {
			if (col1 == "(us)") {
		    	x_unit = "microseconds";
			} else if (col1 == "(ns)") {
				x_unit = "nanoseconds";
			} else if (col1 == "(ms)") {
				x_unit = "milliseconds";
			}
		} else if (line_no > 2) {
			C1.push_back(strtod(col1.c_str(), NULL));
			C2.push_back(strtod(col2.c_str(), NULL));
		}
		line_no++;
	}

	myfile.close();

	double inf = std::numeric_limits<double>::infinity();
	double neg_inf = -std::numeric_limits<double>::infinity();
	auto average = accumulate(C2.begin(), C2.end(), 0.0)/C2.size();

	if (std::find(C2.begin(), C2.end(), inf) != C2.end() || std::find(C2.begin(), C2.end(), neg_inf) != C2.end() || C1.size() == 0) {
		std::cerr << "error: invalid input.\n";
		return EXIT_FAILURE;
	}

	if (average < 1000) {
		for (Int_t i = 0; i < C2.size() ; i++) {
			C2.at(i) = C2.at(i)*-1;
		}
	}

	TH1F *h = new TH1F("h", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		h->SetBinContent(i, C2.at(i));
	}

	TCanvas *c = new TCanvas("c","c", 10, 10, 1000, 900);
	gStyle->SetOptStat(kFALSE);

	h->GetXaxis()->SetTitle(x_unit.c_str());
	h->GetYaxis()->SetTitle("millivolt");
	h->Draw("same");

	std::string plot_name = OutputDirectory + folder + "_" + fileno + filetype_plot;
	c->SaveAs(plot_name.c_str());

	return 0;
}
