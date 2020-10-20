#include <map>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>

using namespace std;
using std::string;
using std::cerr;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


class DataReader {
	
	private:
		map<string, double> options;
		vector< vector<double> > distribution_data;

	public:
		DataReader(char* option_file);
		vector< vector<double> > getDistributionData();
		
		void SetParameters();
        void ReadData(char* data_file);
		void RemoveData();
		
};
