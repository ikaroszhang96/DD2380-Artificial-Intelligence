#ifndef _DUCKS_HMM_HPP_
#define _DUCKS_HMM_HPP_

#include "Action.hpp"
#include <vector>

namespace ducks {
    class HMM {
    public:
        HMM();
        HMM(const std::vector<EMovement> &obs);
        
        void forward_pass(std::vector<std::vector<double>> &alpha, std::vector<double> &c , const std::vector<EMovement> &obs);
        void back_pass(std::vector<double> &c, std::vector<std::vector<double>> &beta, const std::vector<EMovement> &obs);
        void re_estimate(std::vector<EMovement> obs, std::vector<std::vector<double> > gamma, std::vector<std::vector< std::vector<double> > > digamma);
        
        double probSeq(const std::vector<EMovement> &obs);
        void estimateModel(const std::vector<EMovement> &obs);
        std::vector<double> nextEmissions(const std::vector<EMovement> &obs);
        
        std::vector<std::vector<double>> A;
        std::vector<std::vector<double>> B;
        std::vector<double> pi;
  
    };
}

	


#endif

