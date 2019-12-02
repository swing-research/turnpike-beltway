#include "DataReader.h"
#include "global.h"
#include "omp.h"

DataReader::DataReader(char* option_file) {
	// read option files and reset the parameter values
	ifstream read_option;
	read_option.open(option_file);
	string opt_tmp, opt_name;
	double opt_val;
	int opt_idx=0;
	while(read_option>>opt_tmp){
		if (opt_idx==0) {opt_name=opt_tmp; opt_idx++;}
		else {opt_val=atof(opt_tmp.c_str()); opt_idx=0;}
		
		if(opt_idx==0) {options[opt_name]=opt_val;}
	}
	read_option.close();
}

void DataReader::SetParameters() {

	// reset the parameters
	
	num_thread = omp_get_max_threads();

    if (options.find("num_thread")!=options.end())   {num_thread=(int)options["num_thread"];}
	if (options.find("num_smp")!=options.end())   {num_smp=(int)options["num_smp"];}
	if (options.find("max_ite")!=options.end())   {max_ite=(int)options["max_ite"];}
	if (options.find("max_sg_ite")!=options.end())   {max_sg_ite=(int)options["max_sg_ite"];}
	if (options.find("min_space_unit")!=options.end()) {min_space_unit=options["min_space_unit"];}
	if (options.find("tau")!=options.end())     {tau=options["tau"];}
	if (options.find("sigma")!=options.end())     {sigma=options["sigma"];}
	if (options.find("perturb_std")!=options.end())  {perturb_std=options["perturb_std"];}
	if (options.find("sg_tol")!=options.end())  {sg_tol=options["sg_tol"];}
	if (options.find("step_ori")!=options.end())  {step_ori=options["step_ori"];}
	if (options.find("bkt_rate")!=options.end())  {bkt_rate=options["bkt_rate"];}
	if (options.find("step_thd")!=options.end())  {step_thd=options["step_thd"];}
	if (options.find("cvg_thd")!=options.end())   {cvg_thd=options["cvg_thd"];}

    cout<<"Minimum space unit:    "<<min_space_unit<<endl;
    cout<<"No. of threads:        "<<num_thread<<endl;
    cout<<"No. of samples:        "<<num_smp<<endl; 
    cout<<"Max iteration:         "<<max_ite<<endl;
    cout<<"Distance range:        "<<tau<<endl;
    cout<<"Noise std:             "<<sigma<<endl;
    cout<<"Backtracking rate:     "<<bkt_rate<<endl;
    cout<<"Stepsize ori:          "<<step_ori<<endl;
    cout<<"Stepsize threshold:    "<<step_thd<<endl;
    cout<<"Convergence threshold: "<<cvg_thd<<endl;
    cout<<endl;
}

void DataReader::ReadData(char* data_file) {
	
	// read data files and save it in a 2-dimensional vector mat
	ifstream read_data;
	
    read_data.open(data_file);
	istringstream istr;
	string str;
	int line_num=0;
	    
    while(getline(read_data, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        distribution_data.push_back(tmpvec);

	    tmpvec.clear();
	    istr.clear();
	    line_num++;
    }
	
	num_raw_uq_distance = line_num;
	read_data.close();

	cout<<"Reading data finished."<<endl;
}


vector< vector<double> > DataReader::getDistributionData() {return distribution_data;}

void DataReader::RemoveData() {distribution_data.clear();}
