#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

#include "global.h"
#include "uDGP.h"
#include "SED.h"

using namespace std;
using namespace Eigen;

using std::ifstream;
using std::string;

void process_cmdline(char *argv[], char **data_file, char **option_file, char **output_file) {
	
	*data_file	    = argv[1];
	*option_file	= argv[2];
	*output_file	= argv[3];
}

int main(int argc, char *argv[]){

	char *data_file, *option_file, *output_file;
	if (argc==4) {
		process_cmdline(argv, &data_file, &option_file, &output_file);
	} else {
		cout<<"Command line error!"<<endl;
		abort();
	}

	DataReader data_input(option_file);
	data_input.SetParameters();
    data_input.ReadData(data_file);

	uDGP *udgp_gd = new SED();
    udgp_gd->SetOutputFile(output_file);
    udgp_gd->SetData(&data_input);
    udgp_gd->SetMeasureMatrix();
    udgp_gd->Initialization();
    udgp_gd->GradientDescent();
    
    VectorXd smp_pos_out = udgp_gd->GetSamplePos();

    ofstream write_result(output_file, ios_base::trunc);
    for (int i=0; i<smp_pos_out.size(); i++) {
        write_result<<smp_pos_out(i)<<" ";
    }
    write_result<<"\n";
    write_result.close();
    
    delete udgp_gd;

	exit(0);
}
