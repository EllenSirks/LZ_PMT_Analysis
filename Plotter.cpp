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

	TH1F *h1 = (TH1F*) h->Clone();
	TH1F *h2 = (TH1F*) h->Clone();

	// c->Divide(3,1, small, small, small);

	// gStyle->SetPadBorderMode(0);
 //   	gStyle->SetFrameBorderMode(0);

	// c->cd(1);
	// gPad->SetRightMargin(small);
	gStyle->SetOptStat(kFALSE);

	TF1 *f2 = new TF1("f2", "[0]", C1.front(), C1.at(400));
	h->Fit(f2,"R");

	auto background1 = f2->GetParameter(0);
	// threshold1 = (300./400.)*4*f2->GetParError(0);


	TCanvas *c = new TCanvas("c","c", 10, 10, 1000, 900);



	TH1F *h3 = new TH1F("h3", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		// h3->SetBinContent(i, C2.at(i) - background1);
		h3->SetBinContent(i, C2.at(i));

	}

	// h3->SetFillColor(3);
	h3->Draw();
	// h3->Fit(f2,"R");

	// auto background1 = f2->GetParameter(0);
	// threshold1 = 4*f2->GetParError(0);
	threshold1 = 1;

	double addition = 0.;
	std::vector<int> above_bkg_starts, above_bkg_ends, below_bkg_starts, below_bkg_ends;


	for (Int_t i = 0; i < C2.size() - 1; i++) {
		if (C2.at(i) > background1 + addition ) {
			for (Int_t j = i + 1; j < C2.size(); j++) {
				if (C2.at(j) < background1 + addition) {
					above_bkg_starts.push_back(i);
					above_bkg_ends.push_back(j);
					i = j;
					break;
				}
			}
		} else if (C2.at(i) < background1) {
			for (Int_t j = i + 1; j < C2.size(); j++) {
				if (C2.at(j) > background1) {
					below_bkg_starts.push_back(i);
					below_bkg_ends.push_back(j);
					i = j;
					break;
				}
			}
		}
	}


	// Peaks above/below threshold

	std::vector<int> peak_starts, peak_ends, peak_starts_neg, peak_ends_neg;

	for (Int_t i = 0; i < above_bkg_starts.size(); i++) {
		for (Int_t j = above_bkg_starts.at(i); j < above_bkg_ends.at(i); j++) {
			if (C2.at(j) > background1 + threshold1) {
				peak_starts.push_back(above_bkg_starts.at(i));
				peak_ends.push_back(above_bkg_ends.at(i));
				break;
			}
		}
	}

	for (Int_t i = 0; i < below_bkg_starts.size(); i++) {
		for (Int_t j = below_bkg_starts.at(i); j < below_bkg_ends.at(i); j++) {
			if (C2.at(j) < background1 - threshold1) {
				peak_starts_neg.push_back(below_bkg_starts.at(i));
				peak_ends_neg.push_back(below_bkg_ends.at(i));
				break;
			}
		}
	}

	double binwidth;

	if ( x_unit == "microseconds") {
		binwidth = h3->GetBinWidth(2)*1000.;
	} else if (x_unit == "milliseconds") {
		binwidth = h3->GetBinWidth(2)*1000000.;
	} else {
		binwidth = h3->GetBinWidth(2);
	}

	TH1F *h4 = new TH1F("h4", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		h4->SetBinContent(i, C2.at(i) - background1);
		// h4->SetBinContent(i, C2.at(i));

	}

	double peak_area = 0;

	for (Int_t i = peak_starts.at(0) - 2; i < peak_ends.at(0) ; i++) {
		peak_area += h4->GetBinContent(i)*binwidth;
	}

	
	auto legend = new TLegend(0.5, 0.7, 0.9, 0.9);


	TLine *l = new TLine(C1.front(), background1, C1.back(), background1);
	l->SetLineColor(2);
	l->SetLineWidth(2);


	TLine *l2 = new TLine(C1.at(400), -11.34, C1.at(400), 2.7);
	l2->SetLineColor(2);
	l2->SetLineStyle(2);
	l2->SetLineWidth(2);


	TLine *line2 = new TLine(C1.front(), background1+threshold1, C1.back(), background1+threshold1);
	line2->SetLineStyle(2);
	// line2->SetLineWidth(2);	


	legend->AddEntry(l, "Background", "l");
	legend->AddEntry(line2, "Threshold", "l");


	TH1F *h6 = new TH1F("", "", C1.size(), C1.front(), C1.back());
	for (Int_t j = 0; j < C2.size(); j++) {
		if (j < peak_ends.at(0) && j > peak_starts.at(0) - 2) {
			h6->SetBinContent(j, 0);
		} else { 
			h6->SetBinContent(j, C2.at(j));
		}
	}

	h6->SetFillColor(2);
	// h6->Draw("same"); 

	for (Int_t i = 0; i < peak_starts.size() ; i++) {
	TH1F *h5 = new TH1F("", "", C1.size(), C1.front(), C1.back());
		for (Int_t j = 0; j < C2.size(); j++) {
			if (j < peak_ends.at(i) && j > peak_starts.at(i) - 2) {
				h5->SetBinContent(j, C2.at(j));
			} else { 
				h5->SetBinContent(j, background1);
			}
		}

	std::stringstream ss;
	ss << peak_area;
	std::string myString1 = "signal area: " + ss.str() + " mV*ns";

	legend->AddEntry(h5, myString1.c_str(), "f");

	h5->SetFillColor(3);
	// h5->Draw("same"); 
	}



	TH1F *h11 = new TH1F("h11", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		if (C2.at(i) >= 0) {
			h11->SetBinContent(i, C2.at(i));
		}
	}

	// h11->SetLineColor(3);
	h11->SetFillColor(3);
	h11->Draw("same");

	TH1F *h12 = new TH1F("h12", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		if (i < peak_ends.at(0) && i > peak_starts.at(0) -2) {
			h12->SetBinContent(i, background1);
		}
	}

	// h12->SetLineColor(3);
	h12->SetFillColor(3);
	h12->Draw("same");

	TH1F *h7 = new TH1F("h7", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		if (C2.at(i) <= 0 && i < peak_ends.at(0) && i > peak_starts.at(0) -2) {
			h7->SetBinContent(i, C2.at(i));
		}
	}

	// h7->SetLineColor(10);
	h7->SetFillColor(10);
	h7->Draw("same");

	// h3->SetFillColor(0);
	// h3->SetFillColorAlpha(0, 1);
	// h3->SetFillStyle(0);
	// h3->Draw();
	int start, end;
	for (Int_t i = 0; i < C2.size() ; i++) {
		if (C2.at(i) >= 0) {
			start = i;
			break;
		}
	}

	for (Int_t i = 0; i < C2.size() ; i++) {
		if (C2.at(i) <= 0 && i > start) {
			end = i;
			break;
		}
	}


	TH1F *h88 = new TH1F("h88", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		if (i < peak_starts.at(0) -2) {
			h88->SetBinContent(i, C2.at(i));
		} else {
			h88->SetBinContent(i, 0);
		}
	}

	h88->SetLineColor(10);
	h88->SetFillColor(10);
	// h88->Draw("same");  

	TH1F *h89 = new TH1F("h89", "", C1.size(), C1.front(), C1.back());
	for (Int_t i = 0; i < C2.size() ; i++) {
		if (i > peak_ends.at(0)) {
			h89->SetBinContent(i, C2.at(i));
		} else {
			h89->SetBinContent(i, 0);
		}
	}

	h89->SetLineColor(10);
	h89->SetFillColor(10);
	// h89->Draw("same");



	TLine *l5 = new TLine(C1.at(start + 1), 0, C1.at(end -2), 0);
	l5->SetLineColor(3);
	l5->SetLineWidth(3);
	l5->Draw("same");


	l->Draw("same");

	l2->Draw("same");

	line2->Draw("same");

	gPad->RedrawAxis();

	h3->GetYaxis()->SetTitle("millivolt");
	h3->GetXaxis()->SetTitle(x_unit.c_str());

	legend->Draw("same");

	std::string plot_name = OutputDirectory + folder + "_" + fileno + filetype_plot;
	c->SaveAs(plot_name.c_str());

	delete h;

	return 0;
}
