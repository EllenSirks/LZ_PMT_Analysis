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

int Peak_Find(int start_file, int num_files, std::string folder) {

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
			std::string zeros = "00";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else if (i < 100 && i > 9) {
			std::string zeros = "0";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else if (i < 1000 && i > 99) {
			std::string zeros = "";
			fileno = zeros;
			fileno.append(std::to_string(i));
		}
		else if (i < 10000 && i > 999 ) {
			std::string zeros = "";
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
					C2.at(i) = C2.at(i)*-1;
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

			TH1F *h2 = (TH1F*) h1->Clone();
			TF1 *f1 = new TF1("f1", "[0]", C1.front(), C1.at(250));
			h2->Fit(f1,"R");

			auto background = f1->GetParameter(0);

			no_devs = 3;
			auto threshold = no_devs*f1->GetParError(0); // Decide here how many devations from background is threshold
			// auto threshold = 1.25;

			auto chi2 = f1->GetChisquare();

			c1->cd(1);
			TLine *line = new TLine(C1.front(), background + threshold, C1.back(), background + threshold);
			line->SetLineStyle(2);
			line->Draw("same");

			/*** Save Plot if Large Deviation in Background ***/

			if (chi2 > 10000) {

				TCanvas *c2 = new TCanvas("c2","c2", 10, 10, 1000, 900);

				h1->GetXaxis()->SetTitle(x_unit.c_str());
				h1->GetYaxis()->SetTitle("millivolt");
				h1->Draw();

				TLine *l = new TLine(C1.front(), background, C1.at(250), background);
				l->SetLineColor(2);
				l->SetLineWidth(2);
				l->Draw("same");

				std::string chi2_plot_name = OutputDirectory + "plots/chi2/" + fileno + filetype_plot;
				c2->SaveAs(chi2_plot_name.c_str());

				delete c2;
			}

			/*** Find Peaks Above Threshold ***/

			std::vector<int> above_bkg_starts, above_bkg_ends;

			for (Int_t i = 0; i < C2.size() - 1; i++) {
				if (C2.at(i) > background) {
					for (Int_t j = i + 1; j < C2.size(); j++) {
						if (C2.at(j) < background) {
							above_bkg_starts.push_back(i);
							above_bkg_ends.push_back(j);
							i = j;
							break;
						}
					}
				}
			}

			std::vector<int> peak_starts, peak_ends;

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

			no_peaks_per_file.push_back(no_peaks);

			/*** Check which Peak is Actual Pulse ***/

			double max_of_peak = 0;
			int pulse_it;

			for (Int_t i = 0; i < peak_starts.size(); i++) {
				for (Int_t j = peak_starts.at(i); j < peak_ends.at(i); j++) {
					if (C2.at(j) > max_of_peak) {
						max_of_peak = C2.at(j);
						pulse_it = i;
					}
				}
			}

			if (pulse_it != 0) {
				std::cout << "Pulse in Background!" << std::endl;
				pulse_bkg_files.push_back(i);
			}

			int k = 0;
			while (k < pulse_it) {
				peak_starts.erase(peak_starts.begin());
				peak_ends.erase(peak_ends.begin());
				k++;
			}

			/*** Determine Pulse Window if No Peak is Found ***/

			if (peak_starts.size() == 0) {

				double start = 0.0;
				double end = 0.0;

				if (x_unit == "microseconds") {
					end = 15.0/pow(10, 3);
				} else if (x_unit == "milliseconds") {
					end = 15.0/pow(10, 6);
				} else {
					end = 15.0;
				}

				auto it1 = 0;
				auto it2 = 0;

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

			/*** Save Plot if a Peak is Found ***/

			if (no_peaks != 0) {

				TCanvas *c3 = new TCanvas("c3","c3", 10, 10, 1000, 900);

				h1->GetXaxis()->SetTitle(x_unit.c_str());
				h1->GetYaxis()->SetTitle("millivolt");
				h1->Draw();

				TLine *l1 = new TLine(C1.front(), background + threshold, C1.back(), background + threshold);
				l1->SetLineStyle(2);
				l1->SetLineWidth(2);
				l1->Draw("same");

				TLine *l2 = new TLine(C1.front(), background, C1.at(250), background);
				l2->SetLineColor(2);
				l2->SetLineWidth(2);
				l2->Draw("same");

				std::string count_plot_name = OutputDirectory + "plots/count/" + fileno + filetype_plot;
				c3->SaveAs(count_plot_name.c_str());
				c3->SaveAs("test.root");

				delete c3;
			}

			/*** Determine After Pulsing Window (2 us) ***/

			double ap_end;
			if (x_unit == "nanoseconds") {
				ap_end = C1.at(peak_ends.at(0)) + 2*pow(10, 3);
			} else if (x_unit == "microseconds") {
				ap_end = C1.at(peak_ends.at(0)) + 2;
			} else if (x_unit == "milliseconds") {
				ap_end = C1.at(peak_ends.at(0)) + 2/pow(10, 3);
			}

			int ap_end_it = 0;
			for (Int_t i = peak_ends.at(0) + 1; i < C1.size() ; i++) {
				if (C1.at(i) < ap_end) {
					ap_end_it = i;
				}
			}

			/*** Determine if there is any After Pulsing ***/

			std::vector<int> ap_starts, ap_ends;

			for (Int_t i = 1; i < peak_ends.size(); i++) {
				if (peak_starts.at(i) < ap_end_it) {
					ap_starts.push_back(peak_starts.at(i));
					ap_ends.push_back(peak_ends.at(i));
				}
			}

			/*** Determine Peak Height (Max) ***/

			auto max = *std::max_element(C2.begin(), C2.end());
			peaks_max.push_back(max - f1->GetParameter(0));

			/*** Draw Histograms Set at Baseline ***/
				
			c1->cd(2);

			TH1F *h3 = new TH1F("h3", "", C1.size(), C1.front(), C1.back());
			for (Int_t i = 0; i < C2.size() ; i++) {
				h3->SetBinContent(i, C2.at(i) - f1->GetParameter(0));
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

			double peak_area = 0;

			for (Int_t i = peak_starts.at(0); i < peak_ends.at(0) ; i++) {
				peak_area += h3->GetBinContent(i)*binwidth;
			}

			peaks_area.push_back(peak_area);

			/*** Determine After Pulsing Area ***/

			std::vector<double> ap_areas;
			double total_ap_area = 0;

			for (Int_t i = 0; i < ap_starts.size() ; i++) {
				double ap_area = 0;
				for (Int_t j = ap_starts.at(i); j < ap_ends.at(i); j++) {
					ap_area += h3->GetBinContent(j)*binwidth;
				}
				ap_areas.push_back(ap_area);
				total_ap_area += ap_area;
			}

			double total_ap_percentage = (total_ap_area/peak_area)*100;
			total_ap_percentages.push_back(total_ap_percentage);

			/*** Calculate Percentage of Peak Area for After Pulses ***/

			std::vector<double> ap_percentages;

			for (Int_t i = 0; i < ap_areas.size(); i++ ) {
				double percentage = (ap_areas.at(i)/peak_area)*100;
				ap_percentages.push_back(percentage);
			}

			/*** Draw Peak Area ***/

			auto legend = new TLegend(0.525, 0.625, 0.9, 0.9);

			TH1F *h4 = new TH1F("h4", "", C1.size(), C1.front(), C1.back());
			for (Int_t i = 0; i < C2.size() ; i++) {
				if (i < peak_ends.at(0) && i > peak_starts.at(0)) {
					h4->SetBinContent(i, C2.at(i) + abs(f1->GetParameter(0)));
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

			/*** Draw After Pulsing Area ***/

			for (Int_t i = 0; i < ap_starts.size() ; i++) {
				TH1F *h = new TH1F("", "", C1.size(), C1.front(), C1.back());
				for (Int_t j = 0; j < C2.size(); j++) {
					if (j < ap_ends.at(i) && j > ap_starts.at(i)) {
						h->SetBinContent(j, C2.at(j) + abs(f1->GetParameter(0)));
					} else { 
						h->SetBinContent(j, 0);
					}
				}

				stringstream ss;
				ss << ap_percentages.at(i);
				std::string myString = ss.str() + " % of signal area";

				legend->AddEntry(h, myString.c_str(), "f");
				h->SetFillColor(i+2);
				h->Draw("same"); 
			}

			std::stringstream ss1, ss2;
			ss1 << total_ap_area;
			ss2 << total_ap_percentage;
			std::string myString2 = "summed area: " + ss1.str() + " psV - " + ss2.str() + " % of signal area";

   			legend->AddEntry((TObject*)0, myString2.c_str(), "");

			if (no_peaks != 0) {
				if (ap_areas.size() == 0 ) {
					gStyle->SetLegendTextSize(0.05);
				} else {
					gStyle->SetLegendTextSize(0.);
				}
				legend->Draw();
			}

			Float_t ymax = h3->GetMaximum();
			Float_t ymin = h3->GetMinimum();

			TLine *line1 = new TLine(C1.front(), 0, C1.back(), 0);
			line1->SetLineColor(2);
			line1->SetLineWidth(2);
			line1->Draw("same");
			TLine *line2 = new TLine(C1.front(), threshold, C1.back(), threshold);
			line2->SetLineStyle(2);
			line2->Draw("same");

			h3->SetStats(kFALSE);
			h4->SetStats(kFALSE);

			/*** Save Plots ***/

			std::string plot_name = OutputDirectory + HistogramName + fileno + filetype_plot;
			std::string root_name = OutputDirectory + HistogramName + fileno + ".root";

			if (i - start_file + 1 < 101) {
				c1->SaveAs(plot_name.c_str());
				c1->SaveAs(root_name.c_str());
			}

			/*** Delete Canvasses/Histograms ***/

			delete c1;
			delete h1;
			delete h2;
			delete h3;
			delete h4;

			/*** Count Total Time ***/

			total_time += C1.back() - C1.front();
		}
	}

	/*** Create Files ***/

	std::string Out1 = OutputDirectory + "specs/peak_heights_max.txt";

	std::ofstream out1(Out1);
	for (const auto &peak : peaks_max) {
		out1 << peak << "\n";
	}

	std::string Out2 = OutputDirectory + "specs/peak_area.txt";

	std::ofstream out2(Out2);
	for (const auto &area : peaks_area) {
		out2 << area << "\n";
	}

	std::string Out3 = OutputDirectory + "specs/peaks_per_file.txt";

	std::ofstream out3(Out3);
	for (const auto &peaks : no_peaks_per_file) {
		out3 << peaks << "\n";
	}

	std::string Out4 = OutputDirectory + "specs/ap_percentage.txt";

	std::ofstream out4(Out4);
	for (const auto &percentages : total_ap_percentages) {
		out4 << percentages << " % \n";
	}

	std::string Out5 = OutputDirectory + "specs/pulse_bkg_files.txt";

	std::ofstream out5(Out5);
	for (const auto &file : pulse_bkg_files) {
		out5 << file << "\n";
	}

	std::cout << "Skipped: " << skipped_files.size() << " Files" << std::endl;
	std::string err_Out = OutputDirectory + "specs/skipped_files.txt";

	std::ofstream err_out(err_Out);
	for (const auto &file : skipped_files) {
		err_out << file << "\n";
	}
	
	/*** Dark Count Rate ***/

	double total_no_files = num_files - skipped_files.size();
	
	mean_time_per_file = total_time/total_no_files;

	double time_in_seconds;

	if (x_unit == "milliseconds") {
		time_in_seconds = total_time*pow(10, -3);
	} else if (x_unit == "microseconds") {
		time_in_seconds = total_time*pow(10, -6);
	} else if (x_unit == "nanoseconds") {
		time_in_seconds = total_time*pow(10, -9);
	}

	double dark_count_rate = total_no_peaks/time_in_seconds;
	double dark_count_rate_err = sqrt(total_no_peaks)/time_in_seconds;

	std::cout << "Number of Files: " << total_no_files << std::endl;
	std::cout << "Time per File: " << mean_time_per_file << " " << x_unit << std::endl;
	std::cout << "Total Time: " << total_time << " " << x_unit << std::endl;

	std::cout << "Number of Peaks: " << total_no_peaks << std::endl;
	std::cout << "(Dark) Count Rate: " << dark_count_rate << " +/- " << dark_count_rate_err << " Hz" << std::endl;

	std::string dark_count_Out = OutputDirectory + "specs/dark_count.txt";

	std::ofstream dark_count_out(dark_count_Out);
	dark_count_out << dark_count_rate << " +/- " << dark_count_rate_err << " Hz \n";
	dark_count_out << no_devs << " standard devations";

	return 0;
}