#ifndef SED_H
#define SED_H

#include "uDGP.h"
#include <thread>


class SED:public uDGP
{
	private:
	    
	
	public:
		SED();

        
        virtual double ComputeObjFun(VectorXd smp_vect);
        virtual void ComputeObjFunMuti( vector<int> all_distance_block, VectorXd* smp_vec_muti, vector<double>* obj_seq_pt, int val_idx );
        virtual VectorXd ComputeGradient(VectorXd smp_vect);
        //virtual void ComputeGradientMuti( vector<int> all_distance_block, VectorXd* smp_vec_muti, vector<VectorXd>* smp_der_seq_pt, int val_idx );
        virtual void ComputeGradientMuti( vector<int> all_distance_block, VectorXd* smp_vec_muti, VectorXd* smp_der_seq_pt, int val_idx );
		
		virtual ~SED();
};

#endif // SED_H
