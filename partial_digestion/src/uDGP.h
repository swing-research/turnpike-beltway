#ifndef UDGP_H
#define UDGP_H

#include <iostream>
#include <map>
#include <vector>

#include <random>

#include <Eigen/Dense>

#include "DataReader.h"
#include "global.h"

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

typedef std::pair<int,int> mypair;

class uDGP
{    
    protected:
    
        VectorXd valid_idx_pos;
        vector<int> valid_idx_vec;
        vector<int> valid_idx_vec_exclude;  // exclude the two anchor points
        
        vector< vector<double> > raw_distribution;
        vector<int> all_distance;   // discretized distance to be considered
        
        vector<int> all_distance_vec;    // note that this is distance, not squared distance, distance_grid
        vector<double> all_distance_true;   // real distances.
        
        vector<int> D_mat_distance_vec;

        VectorXd smp_pos_init;
        VectorXd smp_pos;
        
        map<int, double> all_distance_diff; // quantized distances between two different points
        
        map<int, double> all_distribution;
        map<int, double> est_distribution;
        vector<mypair> all_block_count;
        map<int, vector<int>> all_partition;
        
        map<int, double> D_mat_norm;
       
        map<int, int> idx_mapping;
        
        char *smp_pos_init_file;
        char *output_ite_file;
        
        int M;
        int md_p1;
        int num_pos;    // number of possible locations
        double step;
        double obj_val;
        double max_distance;
        double domain_sz;
        
        int anchor_one;
        map<int, double> anchor_one_seq;
        int anchor_two;
        map<int, double> anchor_two_seq;
        
        int num_thread_assign;
        
        int num_nev;
        int max_eigs_ite;
        double eigs_tol;
        double perturb_factor;
        
        
    public:
        uDGP();
        void SetOutputFile(char* output_file);
        void SetInitFile(char* init_file);
        void SetData(DataReader* data_reader);
        void SetMeasureMatrix();
        void Initialization();
        void GradientDescent();

        double NormalCdf(double x_val, double mu, double sigma);
        double GetObjFun();
        double ComputeEstDbt( VectorXd* smp_vec, int distance_val);
        void ComputeEstProj(VectorXd* smp_vec, VectorXd* smp_vec_proj, double mut_factor, int distance_val);
        VectorXd GetSamplePos();
        VectorXd ProjectOntoCvxSet(VectorXd smp_vec, int num_smp_proj);

        virtual double ComputeObjFun(VectorXd smp_vec);
        virtual void ComputeObjFunMuti(vector<int> all_distance_block, VectorXd* smp_vec_muti, vector<double>* obj_seq_pt, int val_idx );
        virtual VectorXd ComputeGradient(VectorXd smp_vec);
        virtual void ComputeGradientMuti( vector<int> all_distance_block, VectorXd* smp_vec_muti, VectorXd* smp_der_seq_pt, int val_idx );


        virtual ~uDGP();
        
};

#endif // UDGP_H
