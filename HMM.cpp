#include "HMM.hpp"
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>

#define N 3
#define M COUNT_MOVE
#define MAX_ITERS 10

namespace ducks{
    HMM::HMM():
        A(N, std::vector<double>(N, 0)),
        B(N, std::vector<double>(M, 0)),
        pi(N, 0)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0, 1);
        std::uniform_real_distribution<> pidist(0.30, 0.35);
        
        // Initialize pi
        double sum_row = 0;
        
        for (int i = 0; i < N; i++) {
            pi[i] = pidist(gen);
            sum_row += pi[i];
        }
        // Normalize pi
        for (int i = 0; i < N; i++)
            pi[i] /= sum_row;
        
        
        // Initialize A
        std::vector<double> sum_rows(N);
        
        for (int i = 0; i < N; i++) {
            sum_rows[i] = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    A[i][j] = 0.25;

                } else {
                    A[i][j] = 0.5;

                }
            }
        }

        // Initialize B
        for (int i = 0; i < N; i++) {
            sum_rows[i] = 0;
            for (int j = 0; j < M; j++) {
                B[i][j] = dist(gen);
                sum_rows[i] += B[i][j];
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                B[i][j] /= sum_rows[i];
            }
        }
    }
    
    HMM::HMM(const std::vector<EMovement> &obs) :
        A(N, std::vector<double>(N, 0)),
        B(N, std::vector<double>(M, 0)),
        pi(N, 0)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0, 1);
        std::uniform_real_distribution<> pidist(0.30, 0.35);
        
        // Initialize pi
        double sum_row = 0;
        
        for (int i = 0; i < N; i++) {
            pi[i] = pidist(gen);
            sum_row += pi[i];
        }
        // Normalize pi
        for (int i = 0; i < N; i++)
            pi[i] /= sum_row;
        
        
        // Initialize A
        std::vector<double> sum_rows(N);
        
        for (int i = 0; i < N; i++) {
            sum_rows[i] = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    A[i][j] = 0.25;
                    
                } else {
                    A[i][j] = 0.5;
                    
                }
            }
        }
        
        // Initialize B
        for (int i = 0; i < N; i++) {
            sum_rows[i] = 0;
            for (int j = 0; j < M; j++) {
                B[i][j] = dist(gen);
                sum_rows[i] += B[i][j];
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                B[i][j] /= sum_rows[i];
            }
        }
        
        // Estimate model parameters
        
        if (obs.size() > 1) {
            estimateModel(obs);
        }
    }

    // alpha-pass and beta-pass:
    void HMM::forward_pass( std::vector<std::vector<double>> &alpha, std::vector<double> &c , const std::vector<EMovement> &obs){
			int T = obs.size();
			std::vector<int> O(T);
        	for (int i = 0; i < T; i++) {
             O[i] = obs[i];
            }
            c[0] = 0;
            for(int i = 0; i<N; i++){
                alpha[0][i] = pi[i]* B[i][ O[0] ];
                c[0] += alpha[0][i];
            }
            c[0] = 1/c[0];
            for(int i = 0; i<N; i++)
                alpha[0][i] *= c[0];
            
            for(int t = 1; t<T; t++){
                int cur_emission = O[t];
                c[t] = 0;
                for(int curState = 0; curState<N; curState++){
                    alpha[t][curState] = 0;
                    for(int preState = 0; preState<N; preState++){
                        alpha[t][curState] += alpha[t-1][preState]* A[preState][curState] ;
                    }
                    alpha[t][curState] *=  B[curState][cur_emission];  //B[curState][cur_emission];
                    c[t] += alpha[t][curState];
                }
                c[t] = 1/c[t];
                for(int i = 0; i<N; i++)
                    alpha[t][i] *= c[t];
                
            }
	}

	void HMM::back_pass( std::vector<double> &c, std::vector<std::vector<double>> &beta, const std::vector<EMovement> &obs){
		int T = obs.size();
		std::vector<int> O(T);
        for (int i = 0; i < T; i++) {
             O[i] = obs[i];
        }
		for(int i = 0; i<N ;i++)
                beta[T-1][i] = c[T-1];
            // beta pass: beta_t(i)=sum_j( beta_t+1(j)* a_ij * b_j(O_t+1) )
            for(int t = T-2 ; t>=0; t--){
                int next_emission = O[t+1];
                
                for(int curState = 0; curState<N; curState++){
                    beta[t][curState] = 0;
                    for(int nextState = 0; nextState<N; nextState++)
                        beta[t][curState] += B[nextState][next_emission]* A[curState][nextState]* beta[t+1][nextState];
                    beta[t][curState] *= c[t];
                }
            }
	}

	//re-estimate model:
	void HMM::re_estimate( std::vector<EMovement> obs , std::vector<std::vector<double> > gamma, std::vector<std::vector< std::vector<double> > > digamma){
            //Re-estimate
		int T = obs.size();
		std::vector<int> O(T);
        for (int i = 0; i < T; i++) {
             O[i] = obs[i];
        }
            for(int i = 0; i<N ;i++)
                pi[i] = gamma[0][i];
            double numer = 0;
            double denom = 0;
            //update A:
            for(int i = 0; i<N ;i++){
                for(int j = 0; j<N ;j++){
                    numer = 0;
                    denom = 0;
                    for(int t = 0 ; t<T-1; t++){
                        numer+= digamma[t][i][j];
                        denom+= gamma[t][i];
                    }
                    A[i][j] = numer / denom;
                }
            }
            //update B:
            numer = 0;
            denom = 0;
            for(int i = 0; i<N; i++){
                for(int j = 0; j<M ;j++){
                    numer = 0;
                    denom = 0;
                    for(int t = 0; t<T; t++){
                        if(O[t] == j )
                            numer+= gamma[t][i];
                        denom+= gamma[t][i];
                    }
                    B[i][j] = numer/ denom;
                }
            }
	}
    
    double HMM::probSeq(const std::vector<EMovement> &obs){
        std::vector<double> alpha1(N);
        std::vector<double> alpha(N);
        int T = obs.size();
        
        std::vector<int> O(T);
        for (int i = 0; i < T; i++) {
            O[i] = obs[i];
        }
        
        for (int i = 0; i < N; i++) {
            alpha1[i] = pi[i] * B[i][O[0]];
        }
        
        for (int t = 1; t < T; t++) {
            for (int i = 0; i < N; i++) {
                alpha[i] = 0;
                for (int j = 0; j < N; j++) {
                    alpha[i] += alpha1[j] * A[j][i];
                }
                alpha[i] *= B[i][O[t]];
            }
            alpha1 = alpha;
        }
        
        double prob = 0;
        for (int i = 0; i < N; ++i) {
            prob += alpha[i];
        }
        
        return prob;
    }
    
    void HMM::estimateModel(const std::vector<EMovement> &obs) {
        int T = obs.size();
        
        int iters = 0;
        int maxIters = MAX_ITERS;
        double logProb = -999999990;
        double oldLogProb = -999999999;
        
        std::vector<int> O(T);
        for (int i = 0; i < T; i++) {
            O[i] = obs[i];
        }
        
        std::vector<double> c(T);
        std::vector<std::vector<double>> alpha(T, std::vector<double>(N));
        std::vector<std::vector<double>> beta(T, std::vector<double>(N));
        std::vector<std::vector<std::vector<double>>> digamma(T, std::vector<std::vector<double>>(N, std::vector<double>(N)));
        std::vector<std::vector<double>> gamma(T, std::vector<double>(N));
        
        while (iters < maxIters && logProb > oldLogProb){
             forward_pass(alpha,c, obs);
             back_pass(c,beta, obs);
            
            double denom;
            for(int t = 0; t<T-1 ; t++){
                denom = 0;
                int next_emission = O[t+1];
                for(int i = 0; i<N; i++){
                    for(int j = 0; j<N; j++){
                        denom+= alpha[t][i]*A[i][j]*B[j][next_emission]* beta[t+1][j];
                    }
                }
                
                for(int i = 0; i<N; i++){
                    for(int j = 0; j<N; j++){
                        digamma[t][i][j] = alpha[t][i]*A[i][j]*B[j][next_emission]* beta[t+1][j] / denom;
                    }
                }
            }
            
            for(int t = 0; t<T-1 ; t++){
                for(int i = 0; i<N; i++){
                    gamma[t][i] = 0;
                    for(int j = 0; j<N; j++){
                        gamma[t][i] += digamma[t][i][j];
                    }
                }
            }
            // for special case gamma_t-1(i)
            denom = 0;
            for(int i = 0; i<N; i++){
                denom += alpha[T-1][i];
            }
            for(int i = 0; i<N ;i++)
                gamma[T-1][i] = alpha[T-1][i] / denom;
            
 			// Re-estimate model:
            re_estimate(obs, gamma, digamma);
            
            // Compute log[P(O|lambda)]
            oldLogProb = logProb;
            logProb = 0;
            for (int i = 0; i < T; i++) {
                logProb += log(c[i]);
            }
            logProb = -logProb;
            
            iters += 1;
        }
    }

std::vector<double> HMM::nextEmissions(const std::vector<EMovement> &obs){
	int T = obs.size();
	std::vector<double> O(T);
        for (int i = 0; i < T; i++) {
            O[i] = obs[i];
        }
    std::vector<double> c(T);
    std::vector<std::vector<double>> alpha(T, std::vector<double>(N));
    c[0] = 0;
            for(int i = 0; i<N; i++){
                alpha[0][i] = pi[i]* B[i][ O[0] ];
                c[0] += alpha[0][i];
            }
            c[0] = 1/c[0];
            for(int i = 0; i<N; i++)
                alpha[0][i] *= c[0];
            
            for(int t = 1; t<T; t++){
                int cur_emission = O[t];
                c[t] = 0;
                for(int curState = 0; curState<N; curState++){
                    alpha[t][curState] = 0;
                    for(int preState = 0; preState<N; preState++){
                        alpha[t][curState] += alpha[t-1][preState]* A[preState][curState] ;
                    }
                    alpha[t][curState] *=  B[curState][cur_emission];  //B[curState][cur_emission];
                    c[t] += alpha[t][curState];
                }
                c[t] = 1/c[t];
                for(int i = 0; i<N; i++)
                    alpha[t][i] *= c[t];
                
            }
        std::vector<double> res(M, 0);
        for(int emission = 0; emission < M; emission++){
        	for(int i = 0; i<N; i++){
        		for(int j = 0; j<N; j++){
        			res[emission] += alpha[T-1][i]* A[i][j]*B[j][emission];

        		}
        	}

        }

        return res;

}


}
