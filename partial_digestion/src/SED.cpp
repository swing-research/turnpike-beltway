#include "SED.h"

// Could remove the columns that are all zeros and speed up the computation

SED::SED(){}

void SED::ComputeObjFunMuti( vector<int> all_distance_block, VectorXd* smp_vec_muti, vector<double>* obj_seq_pt, int val_idx ) {
    
    double obj_tmp = 0;
    for (int i=0; i<all_distance_block.size(); i++) {
        est_distribution[all_distance_block[i]] = ComputeEstDbt(smp_vec_muti, all_distance_block[i]);
        obj_tmp += pow( est_distribution[all_distance_block[i]] - all_distribution[all_distance_block[i]] , 2);
    }
    (*obj_seq_pt)[val_idx] = obj_tmp;
    
}

double SED::ComputeObjFun(VectorXd smp_vec) {
    
    vector<double> obj_seq(num_thread_assign, 0);
    
    if (num_thread_assign>1) {
        thread *multi_thread = new thread[num_thread_assign-1];
        for (int i=0; i<num_thread_assign-1; i++) {
            multi_thread[i] = thread(&SED::ComputeObjFunMuti, this,  all_partition[i], &smp_vec, &obj_seq, i);
        }
        
        ComputeObjFunMuti( all_partition[num_thread_assign-1], &smp_vec, &obj_seq, num_thread-1);
        
        for (int i=0; i<num_thread_assign-1; i++) {
            multi_thread[i].join();
        }
        
        delete [] multi_thread;
        
    } else {
        ComputeObjFunMuti( all_partition[num_thread_assign-1], &smp_vec, &obj_seq, num_thread-1);
    }
    
    double obj = 0;
    for (int i=0; i<num_thread_assign; i++) {
        obj += obj_seq[i];
    }

    return obj;
}

void SED::ComputeGradientMuti( vector<int> all_distance_block, VectorXd* smp_vec_muti, VectorXd* smp_der_seq_pt, int val_idx ) {

    VectorXd smp_der_seq_tmp = VectorXd::Zero(M);
    
    for (int i=0; i<all_distance_block.size(); i++) {

        double mut_factor = (2.0*( est_distribution[all_distance_block[i]] - all_distribution[all_distance_block[i]] ) );
        ComputeEstProj(smp_vec_muti, smp_der_seq_pt, mut_factor, all_distance_block[i]);
    }

}

VectorXd SED::ComputeGradient(VectorXd smp_vec) {

    vector<VectorXd> smp_der_seq(num_thread_assign, VectorXd::Zero(M));
    
    if (num_thread_assign>1) {
        thread *multi_thread = new thread[num_thread_assign-1];
        for (int i=0; i<num_thread_assign-1; i++) {
            multi_thread[i] = thread(&SED::ComputeGradientMuti, this,  all_partition[i], &smp_vec, &smp_der_seq[i], i);
        }
        
        ComputeGradientMuti( all_partition[num_thread_assign-1], &smp_vec, &smp_der_seq[num_thread_assign-1], num_thread_assign-1 );
        
        for (int i=0; i<num_thread_assign-1; i++) {
            multi_thread[i].join();
        }
        
        delete [] multi_thread;
        
    } else {
        ComputeGradientMuti( all_partition[num_thread_assign-1], &smp_vec, &smp_der_seq[num_thread_assign-1], num_thread_assign-1 );
    }

    VectorXd smp_der = VectorXd::Zero(M);
    for (int i=0; i<num_thread_assign; i++) {
        smp_der.noalias() += smp_der_seq[i];
    }
    return smp_der;
}

SED::~SED() {}
