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

int Add_Waveforms(int start_file, int num_files, std::string folder) {

	gROOT->SetBatch(kTRUE);
	gErrorIgnoreLevel = kWarning;
	gPrintViaErrorHandler = kTRUE;

	/*** Choose Channel ***/

	std::string input = "";

	std::cout << "What channel (A or B) do you wish to analyse?" << std::endl;
	getline(cin, input);
 	std::cout << "You entered channel: " << input << std::endl; 

	/*** Create Paths ***/

	std::string OutputDirectory = "/home/ellen/Desktop/MPhys_Project/Output/" + folder + "_output/";
	// std::string DataDirectory = "/media/Ellen/Elements/PMT_Data/" + folder + "/";
	std::string DataDirectory = "/home/ellen/Downloads/" + folder + "/";
	std::string DataFileName = folder + "_";
	std::string HistogramName = input + "/histograms/" + folder + "_";
	std::string filetype_data = ".csv", filetype_plot = ".png";

	std::vector<int> skipped_files;
	std::vector<double> final_arr, time;
	std::vector< std::vector<double> > arrays;
	std::string x_unit, col_name;
	int size;

	/*** Create Correct Filenumber ***/

	times = 0;

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
		// std::cout << filename << std::endl;

		/*** Read in Data ***/

		std::ifstream myfile(filename);
		std::vector<double> C1, C2;
		int line_no = 0;

		for (std::string col1, col2, col3; std::getline(myfile, col1, ',') && std::getline(myfile, col2, ',') && std::getline(myfile, col3);) {
			if (line_no == 1) {
				if (col1 == "(us)") {
		    		x_unit = "microseconds";
		    		col_name = "(us)";
				} else if (col1 == "(ns)") {
					x_unit = "nanoseconds";
					col_name = "(ns)";
				} else if (col1 == "(ms)") {
					x_unit = "milliseconds";
					col_name = "(ms)";
				}
			} else if (line_no > 1) {
				C1.push_back(strtod(col1.c_str(), NULL));
				if ( input == "A") {
 					C2.push_back(strtod(col2.c_str(), NULL));
 				} else {
 					C2.push_back(strtod(col3.c_str(), NULL));
 				}
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

			arrays.push_back(C2);

			if (times == 0) {
				size = C1.size();
				time = C1;
			}
	
		}

		times += 1;
	}

	for (Int_t i = 0; i < size; i++) {
		final_arr.push_back(0);
	} 

	for (Int_t i = 0; i < size; i++) {
		for (Int_t j = 0; j < arrays.size(); j++) {
			final_arr.at(i) += arrays.at(j).at(i);
		}
	}

	std::string Out_sum = OutputDirectory  + input + "/specs/waveforms_summed_" + input + filetype_data;
    std::ofstream out_sum(Out_sum);
	for (Int_t i = 0; i < final_arr.size(); i++) {
		if (i == 0) {
			out_sum << "Time,Channel " << input << "\n";
			out_sum << col_name << "," << "(mV)" << "\n" << "\n";
		}
		out_sum << time.at(i) << "," << final_arr.at(i) << "\n";
	}

	std::string Out_sum_av = OutputDirectory  + input + "/specs/waveforms_summed_average_" + input + filetype_data;
    std::ofstream out_sum_av(Out_sum_av);
	for (Int_t i = 0; i < final_arr.size(); i++) {
		if ( i == 0 ) {
			out_sum_av << "Time,Channel " << input << "\n";
			out_sum_av << col_name << "," << "(mV)" << "\n" << "\n";
		}
		out_sum_av << time.at(i) << "," << final_arr.at(i)/(num_files - skipped_files.size()) << "\n";
	}

	// std::cout << num_files - skipped_files.size() << std::endl;
	// std::cout << arrays.size() << std::endl;

	TCanvas *c1 = new TCanvas("c1","c1", 10, 10, 1000, 900);

	TH1F *h1 = new TH1F("h1", "", time.size(), time.front(), time.back());
	for (Int_t i = 0; i < final_arr.size() ; i++) {
		h1->SetBinContent(i, -1*final_arr.at(i));
	}

	std::string title1 = "summed";

	h1->GetXaxis()->SetTitle(x_unit.c_str());
	h1->GetYaxis()->SetTitle("millivolt");
	h1->SetTitle(title1.c_str());
	h1->Draw();

	h1->GetYaxis()->SetTitleOffset(1.4);

	std::string plot_name1 = OutputDirectory + HistogramName + "summed_" + input + filetype_plot;
	c1->SaveAs(plot_name1.c_str());

	TCanvas *c2 = new TCanvas("c2","c2", 10, 10, 1000, 900);

	TH1F *h2 = new TH1F("h2", "", time.size(), time.front(), time.back());
	for (Int_t i = 0; i < final_arr.size() ; i++) {
		h2->SetBinContent(i, -1*final_arr.at(i)/(num_files - skipped_files.size()));
	}

	std::string title2 = "summed & averaged";

	h2->GetXaxis()->SetTitle(x_unit.c_str());
	h2->GetYaxis()->SetTitle("millivolt");
	h2->SetTitle(title2.c_str());
	h2->Draw();

	h2->GetYaxis()->SetTitleOffset(1.4);

	std::string plot_name2 = OutputDirectory + HistogramName + "summed_averaged_" + input + filetype_plot;
	c2->SaveAs(plot_name2.c_str());

	std::cout << "Skipped: " << skipped_files.size() << " Files" << std::endl;

	delete c1;
	delete h1;

	return 0;
}