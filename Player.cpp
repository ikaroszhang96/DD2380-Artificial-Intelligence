#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdio.h>
#include <algorithm>
#include "HMM.hpp"
#include "Constants.hpp"

namespace ducks
{
    Player::Player():
    hmms(COUNT_SPECIES, std::vector<HMM>())
{
    
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */
    // This line choose not to shoot
    //return cDontShoot;
    int nBirds = pState.getNumBirds();
    if(pState.getRound() == 0)
        return cDontShoot;
    
    for(int birdIndex = 0; birdIndex < nBirds; birdIndex++){
        Bird bird = pState.getBird(birdIndex);
        if(bird.getSeqLength() < 70 || !bird.isAlive()) continue;
        std::vector<EMovement> obs = getObservations70(bird);
        ESpecies birdSpecies = predict(bird);  //result of guess class
        if(birdSpecies == SPECIES_BLACK_STORK || birdSpecies == SPECIES_UNKNOWN)
        	continue;

        std::map<int, int> count_Movement;
        double max_max = 0;    //shooting threshold.
        for(auto birdModel : hmms[birdSpecies]){
        	std::vector<double> next_probs = birdModel.nextEmissions(obs);
        	int max_index = 0;
        	double max_prob = 0;
        	for(int i = 0; i<next_probs.size(); i++){
        		if( next_probs[i] > max_prob ){
        			max_prob = next_probs[i];
        			max_index = i;
        		}
        	}
            if(max_prob > max_max){
                max_max = max_prob;
            }
        	count_Movement[max_index]++;
        }
        if (max_max < 1.0e-280){
            continue;
        }
        double thres_move = 0.8;
        int max_times = 0;
        int optimal_MoveIndex = 0;
        std::map<int, int>::iterator iter = count_Movement.begin();
        while(iter != count_Movement.end() ){
        	if(max_times < iter->second ){
        		max_times = iter->second;
        		optimal_MoveIndex = iter->first;
        	}
        	iter++;
        }
        if( (double)( max_times/hmms[birdSpecies].size() ) > thres_move ){
        	return Action(birdIndex, EMovement(optimal_MoveIndex) );
        }
    }
    
    
    return cDontShoot;
    //This line would predict that bird 0 will move right and shoot at it
    //return Action(0, MOVE_RIGHT);
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */
    int nBirds = pState.getNumBirds();
    std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);
    if (pState.getRound() == 0)
    {
        for (int i = 0; i < nBirds; i++)
            lGuesses[i] = ESpecies(SPECIES_PIGEON);
    }
    else{
        for (int i = 0; i < nBirds; ++i) {
            Bird bird = pState.getBird(i);
            lGuesses[i] = predict(bird);
        }
    }

    return lGuesses;
}

void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */
    int nBirds = pSpecies.size();
    std::vector<EMovement> observations;
    
    for(int birdIndex = 0; birdIndex < nBirds; birdIndex++)
    {
        int species = pSpecies[birdIndex];
        if (pSpecies[birdIndex] == -1)
            continue;
        observations = getObservations(pState.getBird(birdIndex));
        
        hmms[species].push_back(HMM(observations));
    }
}

std::vector<EMovement> Player::getObservations(const Bird &b){
    int T = b.getSeqLength();
    int count = 0;
    std::vector<EMovement> sequence(T);
    for(int t = 0; t < T; t++){
        if(b.wasAlive(t)){
            sequence[t] = b.getObservation(t);
            count++;
        }
        else{
            break;
        }
    }
    std::vector<EMovement> sequencealive(count);
    
    for(int t = 0; t < count; t++)
        if(b.wasAlive(t))
            sequencealive[t] = sequence[t];
        
    return sequencealive;
}
    
    std::vector<EMovement> Player::getObservations70(const Bird &b){
        int T = b.getSeqLength();
        std::vector<EMovement> sequence(70);
        for(int t = 0 ; t < 70; t++){
            if(b.wasAlive(t)){
                sequence[t] = b.getObservation(T-70+t);
            }
        }
        return sequence;
    }



ESpecies Player::predict(const Bird &b) {   //average
    ESpecies species = SPECIES_PIGEON;
    std::vector<EMovement> movements = getObservations(b);
    double prob_species;
    int hmmsize;
    double max = 0;
    for (int i = 0; i < COUNT_SPECIES; ++i) {
        prob_species = 0;
        hmmsize = hmms[i].size();
        if (hmmsize > 1) {
            for (int j = 0; j < hmmsize ; ++j) {
                prob_species += hmms[i][j].probSeq(movements);
            }
            prob_species /= hmmsize;
            if (prob_species > max) {
                max = prob_species;
                species = (ESpecies)i;
            }
        }
    }
    if (max < 0) {
        species = SPECIES_BLACK_STORK;
    }
    return species;
}


/*
ESpecies Player::predict(const Bird &b) {   // top k in N
    ESpecies species = SPECIES_PIGEON;
    int k = 0;
    int num_allProb = 0;
    std::vector<EMovement> movements = getObservations(b);
    std::map<ESpecies, std::vector<double> > species_prob;
    std::map<ESpecies, int> count_prob;
    std::vector<double> all_prob;
    for (int i = 0; i < COUNT_SPECIES; ++i) {

            for (int j = 0; j < hmms[i].size() ; ++j) {
            	num_allProb += 1;
                species_prob[ (ESpecies)i ].push_back(hmms[i][j].probSeq(movements) );
                all_prob.push_back( hmms[i][j].probSeq(movements) );
            }

        
    }
    k = num_allProb/2;
    sort(all_prob.begin(), all_prob.end());
    std::reverse(all_prob.begin(), all_prob.end());
    for(int i = 0; i<k; i++){
    	std::map<ESpecies, std::vector<double> >::iterator iter = species_prob.begin();
    	while(iter != species_prob.end() ){
    		std::vector<double> tmp = iter->second;
    		if( find(tmp.begin(), tmp.end(), all_prob[i]) != tmp.end() ){
    			count_prob[iter->first]++;
    			break;
    		}
    		iter++;
    	}
    }
    int max_num = 0;
    std::map<ESpecies, int>::iterator itt= count_prob.begin();
    while(itt != count_prob.end()){
    	if(itt->second > max_num){
    		max_num = itt->second;
    		species = itt->first;
    	}
    	itt++;
    } 

    return species;
}
*/

} /*namespace ducks*/
