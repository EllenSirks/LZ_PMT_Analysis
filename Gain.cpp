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

int Gain(int start_file, int num_files, std::string folder) {

	gROOT->SetBatch(kTRUE);
	gErrorIgnoreLevel = kWarning;
	gPrintViaErrorHandler = kTRUE;

	/*** Create Paths ***/

	std::string OutputDirectory = "/home/ellen/Desktop/MPhys_Project/Output/" + folder + "_output/";
	std::string DataDirectory = "/media/Ellen/Elements/PMT_Data/" + folder + "/";
	std::string DataFileName = folder + "_";
	std::string HistogramName = "histograms/" + folder + "_";
	std::string filetype_data = ".csv", filetype_plot = ".png";

	std::vector<int> skipped_files, no_peaks_per_file, pulse_bkg_files;
	std::vector<double> peaks_fit, peaks_max, peaks_area, total_ap_percentages;
	std::string x_unit;

	double time_per_file, no_devs, threshold, mean_time_per_file, total_time;
	int total_no_peaks = 0;

	/*** Create Correct Filenumber ***/

	for (i = start_file; i < start_file + num_files; i++) {

		std::string fileno;

		if (i < 10) {
			std::string zeros = "0000";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else if (i < 100 && i > 9) {
			std::string zeros = "000";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else if (i < 1000 && i > 99) {
			std::string zeros = "00";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else if (i < 10000 && i > 999 ) {
			std::string zeros = "0";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else {
			fileno = std::to_string(i);
		}

		/*** Create File Path ***/

		std::string filename = DataDirectory + DataFileName + fileno + filetype_data;

		/*** Read in Data ***/

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

		/*** Check Percentage of Files Done ***/

		if (num_files > 99) {
			if ((i - start_file + 1) % (num_files/100) == 0) {
				std::cout << ((i - start_file + 1)/(double)num_files)*100 << "%" << std::endl;
			}
		} else if (num_files < 99 && num_files > 9) {
			if ((i - start_file + 1) % (num_files/10) == 0) {
				std::cout << ((i - start_file + 1)/(double)num_files)*100 << "%" << std::endl;
			}
		}

		/*** Check if Files have Valid Data ***/

		double inf = std::numeric_limits<double>::infinity();
		double neg_inf = -std::numeric_limits<double>::infinity();
		auto average = accumulate(C2.begin(), C2.end(), 0.0)/C2.size();

		if (std::find(C2.begin(), C2.end(), inf) != C2.end() || std::find(C2.begin(), C2.end(), neg_inf) != C2.end() || C1.size() == 0) {
			skipped_files.push_back(i);
			continue;

		} else {

			TCanvas *c1 = new TCanvas("c1","c1", 10, 10, 1000, 900);
			c1->Divide(1,2);
			c1->cd(1);

			/*** Flip Signal ***/

			if (average < 1000) {
				for (Int_t i = 0; i < C2.size() ; i++) {
					if (C2.at(i) != 0) {
						C2.at(i) = C2.at(i)*-1;
					}
				}
			}

			TH1F *h1 = new TH1F("h1", "", C1.size(), C1.front(), C1.back());
			for (Int_t i = 0; i < C2.size() ; i++) {
				h1->SetBinContent(i, C2.at(i));
			}

			std::string title = folder;
			title.append("_"), title.append(fileno);

			h1->GetXaxis()->SetTitle(x_unit.c_str());
			h1->GetYaxis()->SetTitle("millivolt");
			h1->SetTitle(title.c_str());
			h1->Draw();

			/*** Determine Background & Deviation ***/

			double background;
			double threshold;
			double bkg_fit_end;

			no_devs = 4; // Change for different threshold

			TH1F *h2 = (TH1F*) h1->Clone();

			auto sum = accumulate(C2.begin(), C2.begin()+ 400, 0.0)/400.;

			if (abs(sum) < 0.02) {
				background = 0;
				threshold = 5.;

			} else {

				TF1 *f1 = new TF1("f1", "[0]", C1.front(), C1.at(400));
				h2->Fit(f1,"R");

				auto chi2 = f1->GetChisquare();

				/*** Save Plot if Large Deviation in Background ***/

				if (chi2 > 1.0) {

					TF1 *f2 = new TF1("f2", "[0]", C1.front(), C1.at(300));
					h2->Fit(f2,"R");

					bkg_fit_end = C1.at(300);

					background = f2->GetParameter(0);
					threshold = (3./4.)*no_devs*f2->GetParError(0);

					TCanvas *c2 = new TCanvas("c2","c2", 10, 10, 1000, 900);

					h1->GetXaxis()->SetTitle(x_unit.c_str());
					h1->GetYaxis()->SetTitle("millivolt");
					h1->Draw();

					TLine *l = new TLine(C1.front(), f1->GetParameter(0), C1.at(300), f1->GetParameter(0));
					l->SetLineColor(2);
					l->SetLineWidth(2);
					l->Draw("same");

					std::string chi2_plot_name = OutputDirectory + "plots/chi2/" + fileno + filetype_plot;
					// c2->SaveAs(chi2_plot_name.c_str());

					delete c2;

				} else {

					bkg_fit_end = C1.at(400);

					background = f1->GetParameter(0);
					threshold = no_devs*f1->GetParError(0);

				}
			}

			// threshold = 1.25; // Comment out if want to use std dev as threshold

			c1->cd(1);

			TLine *line_pos = new TLine(C1.front(), background + threshold, C1.back(), background + threshold);
			line_pos->SetLineStyle(2);
			line_pos->Draw("same");

			TLine *line_neg = new TLine(C1.front(), background - threshold, C1.back(), background - threshold);
			line_neg->SetLineStyle(2);
			line_neg->Draw("same");

			/*** Find Peaks Above Threshold ***/

			// Peaks above/below background

			std::vector<int> above_bkg_starts, above_bkg_ends, below_bkg_starts, below_bkg_ends;

			double addition = 0.;

			for (Int_t i = 0; i < C2.size() - 1; i++) {
				if (C2.at(i) > background + addition ) {
					for (Int_t j = i + 1; j < C2.size(); j++) {
						if (C2.at(j) < background + addition) {
							above_bkg_starts.push_back(i);
							above_bkg_ends.push_back(j);
							i = j;
							break;
						}
					}
				} else if (C2.at(i) < background) {
					for (Int_t j = i + 1; j < C2.size(); j++) {
						if (C2.at(j) > background) {
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

			int no_peaks = 0;
			for (Int_t i = 0; i < above_bkg_starts.size(); i++) {
				for (Int_t j = above_bkg_starts.at(i); j < above_bkg_ends.at(i); j++) {
					if (C2.at(j) > background + threshold) {
						peak_starts.push_back(above_bkg_starts.at(i));
						peak_ends.push_back(above_bkg_ends.at(i));
						no_peaks += 1;
						total_no_peaks += 1;
						break;
					}
				}
			}

			for (Int_t i = 0; i < below_bkg_starts.size(); i++) {
				for (Int_t j = below_bkg_starts.at(i); j < below_bkg_ends.at(i); j++) {
					if (C2.at(j) < background - threshold) {
						peak_starts_neg.push_back(below_bkg_starts.at(i));
						peak_ends_neg.push_back(below_bkg_ends.at(i));
						break;
					}
				}
			}

			no_peaks_per_file.push_back(no_peaks);

			/*** Check which Peak is in Window ***/

			// Set peak window

			double start = 45.;
			double end = 70.;

			if (no_peaks > 0) {

				auto end_loop = peak_starts.size();

				int index = 0;
				for (Int_t i = 0; i < end_loop; i++) {
					if (C1.at(peak_starts.at(i - index)) < start) {
						peak_starts.erase(peak_starts.begin() + i - index);
						peak_ends.erase(peak_ends.begin() + i - index);
						index += 1;
					}
				}

				end_loop = peak_starts.size();

				index = 0;
				for (Int_t i = 0; i < end_loop; i++) {
					if (C1.at(peak_starts.at(i - index)) > end) {
						peak_starts.erase(peak_starts.begin() + i - index);
						peak_ends.erase(peak_ends.begin() + i - index);
						index += 1;
					}
				}
			}

			/*** Determine Pulse Window if No Peak is Found ***/

			if (peak_starts.size() == 0) { 

				if (x_unit == "microseconds") {
					end = end/pow(10, 3);
				} else if (x_unit == "milliseconds") {
					end = end/pow(10, 6);
				}

				auto it1 = 0;
				auto it2 = 0;

				// Comment out if want zero peak area for noise

				for (Int_t i = 0; i < C1.size() ; i++) {
					if (C1.at(i) < start) {
						it1 = i;
					}
				}

				for (Int_t i = 0; i < C1.size() ; i++) {
					if (C1.at(i) > start && C1.at(i) < end) {
						it2 = i;
					}
				}

				peak_starts.push_back(it1);
				peak_ends.push_back(it2);

			}

			/*** Draw Histograms Set at Baseline ***/
				
			c1->cd(2);

			TH1F *h3 = new TH1F("h3", "", C1.size(), C1.front(), C1.back());
			for (Int_t i = 0; i < C2.size() ; i++) {
				h3->SetBinContent(i, C2.at(i) - background);
			}

			h3->GetXaxis()->SetTitle(x_unit.c_str());
			h3->GetYaxis()->SetTitle("millivolt");
			h3->Draw("same");

			auto binwidth = h3->GetBinWidth(2);

			if (x_unit == "microseconds") {
				binwidth = binwidth * 1000.0;
			} else if (x_unit == "milliseconds") {
				binwidth = binwidth * 1000000.0;
			}

			/*** Determine Peak Area ***/

			double index5 = 0;
			double window = 10; // nanoseconds

			if (x_unit == "microseconds") {
				window = window/pow(10, 3);
			} else if (x_unit == "milliseconds") {
				window = window/pow(10, 6);
			}

			if (peak_ends.at(0) != 0) {
				for (Int_t i = 0; i < C1.size(); i++) {
					if (C1.at(i) > C1.at(peak_starts.at(0)) + window) {
						index5 = i;
						break;
					}
				}
			}

			double peak_area = 0;

			for (Int_t i = peak_starts.at(0); i < peak_ends.at(0) ; i++) {
			// for (Int_t i = peak_starts.at(0); i < index5 ; i++) { // Change if want to set integration range by hand
				peak_area += h3->GetBinContent(i)*binwidth;
			}

			peaks_area.push_back(peak_area);

			/*** Draw Peak Area ***/

			auto legend = new TLegend(0.625, 0.67, 0.9, 0.9);

			TH1F *h4 = new TH1F("h4", "", C1.size(), C1.front(), C1.back());
			for (Int_t i = 0; i < C2.size() ; i++) {
				if (i < peak_ends.at(0) && i > peak_starts.at(0)) { 
				// if (i < index5 && i > peak_starts.at(0)) { // Change if set integration range by hand
					h4->SetBinContent(i, C2.at(i) + abs(background));
				} else { 
					h4->SetBinContent(i, 0);
				}
			}

			std::stringstream ss;
			ss << peak_area;
			std::string myString1 = "signal area: " + ss.str() + " psV";

			legend->AddEntry(h4, myString1.c_str(), "f");
			h4->SetFillColor(42);
			h4->Draw("same");

			legend->Draw();

			Float_t ymax = h3->GetMaximum();
			Float_t ymin = h3->GetMinimum();

			TLine *line1 = new TLine(C1.front(), 0, C1.back(), 0);
			line1->SetLineColor(2);
			line1->SetLineWidth(2);
			line1->Draw("same");

			TLine *line2 = new TLine(C1.front(), threshold, C1.back(), threshold);
			line2->SetLineStyle(2);
			line2->Draw("same");

			TLine *line3 = new TLine(C1.front(), -1*threshold, C1.back(), -1*threshold);
			line3->SetLineStyle(2);
			line3->Draw("same");

			h3->SetStats(kFALSE);
			h4->SetStats(kFALSE);

			/*** Save Plots ***/

			std::string plot_name = OutputDirectory + HistogramName + fileno + "_peak_area" + filetype_plot;
			std::string root_name = OutputDirectory + HistogramName + fileno + "_peak_area.root";

			if (i - start_file + 1 < 21) {
				c1->SaveAs(plot_name.c_str());
				// c1->SaveAs(root_name.c_str());
			}

			/*** Delete Canvasses/Histograms ***/

			delete c1;
			delete h1;
			delete h2;
			delete h3;
			delete h4;
		}
	}

	/*** Create Files ***/

	std::string Out = OutputDirectory + "specs/peak_area.txt";

	std::ofstream out(Out);
	for (const auto &area : peaks_area) {
		out << area << "\n";
	}

	std::cout << "Skipped: " << skipped_files.size() << " Files" << std::endl;
	std::string err_Out = OutputDirectory + "specs/skipped_files.txt";

	std::ofstream err_out(err_Out);
	for (const auto &file : skipped_files) {
		err_out << file << "\n";
	}

	return 0;
}