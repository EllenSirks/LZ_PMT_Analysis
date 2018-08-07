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

int Dark_Current(int start_file, int num_files, std::string folder) {

	gROOT->SetBatch(kTRUE);
	gErrorIgnoreLevel = kWarning;
	gPrintViaErrorHandler = kTRUE;

	/*** Create Paths ***/

	std::string OutputDirectory = "/home/ellen/Desktop/MPhys_Project/Output/" + folder + "_output/";
	// std::string DataDirectory = "/media/Ellen/Elements/PMT_Data/" + folder + "/";
	std::string DataDirectory = "/home/ellen/Downloads/" + folder + "/";
	std::string DataFileName = folder + "_";
	std::string filetype_data = ".csv", filetype_plot = ".png";

	std::vector<int> skipped_files, no_peaks_per_file;
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

			/*** Determine Background & Deviation ***/

			double background;
			double threshold;
			double bkg_fit_end;

			no_devs = 4; // Change for different threshold

			TH1F *h2 = (TH1F*) h1->Clone();

			TF1 *f1 = new TF1("f1", "[0]", C1.front(), C1.at(400));
			h2->Fit(f1,"R");

			auto chi2 = f1->GetChisquare();

			if (chi2 > 1.0) {

				TF1 *f2 = new TF1("f2", "[0]", C1.front(), C1.at(300));
				h2->Fit(f2,"R");

				bkg_fit_end = C1.at(300);

				background = f2->GetParameter(0);
				threshold = (3./4.)*no_devs*f2->GetParError(0);

			} else {

				bkg_fit_end = C1.at(400);

				background = f1->GetParameter(0);
				threshold = no_devs*f1->GetParError(0);

			}

			// threshold = 1.5; // Set Threshold by Hand

			/*** Find Peaks Above Threshold ***/

			std::vector<int> above_bkg_starts, above_bkg_ends, below_bkg_starts, below_bkg_ends;

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

			if (no_peaks > 0) {

				TCanvas *c2 = new TCanvas("c2","c2", 10, 10, 1000, 900);

				h1->GetXaxis()->SetTitle(x_unit.c_str());
				h1->GetYaxis()->SetTitle("millivolt");
				h1->Draw();

				TLine *l = new TLine(C1.front(), background, C1.at(400), background);
				l->SetLineColor(2);
				l->SetLineWidth(2);
				l->Draw("same");

				TLine *l2 = new TLine(C1.front(), background+threshold, C1.back(), background+threshold);
				l2->SetLineWidth(2);
				l2->SetLineStyle(2);
				l2->Draw("same");

				std::string count_plot_name = OutputDirectory + input +  "/plots/count/" + fileno + "_" + input + filetype_plot;
				c2->SaveAs(count_plot_name.c_str());

				delete c2;

			}

			/*** Delete Canvasses/Histograms ***/

			delete h1;
			delete h2;

			/*** Count Total Time ***/

			total_time += C1.back() - C1.front();
		}
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

	/*** Create Files ***/

	std::string Out = OutputDirectory + input + "/specs/peaks_per_file_" + input + ".txt";

	std::ofstream out(Out);
	for (const auto &peaks : no_peaks_per_file) {
		out << peaks << "\n";
	}

	std::cout << "Skipped: " << skipped_files.size() << " Files" << std::endl;

	std::string err_Out = OutputDirectory + input + "/specs/skipped_files_" + input + ".txt";
	std::ofstream err_out(err_Out);
	for (const auto &file : skipped_files) {
		err_out << file << "\n";
	}

	std::string dark_count_Out = OutputDirectory + input + "/specs/dark_count_" + input + ".txt";
	std::ofstream dark_count_out(dark_count_Out);
	dark_count_out << dark_count_rate << " +/- " << dark_count_rate_err << " Hz \n";
	dark_count_out << no_devs << " standard devations";
	
	return 0;
}