/*
 * Population.h
 *
 *  Created on: Apr 6, 2009
 *      Author: David Briemann
 */

#ifndef POPULATION_H_
#define POPULATION_H_

//#define RENEWAL_SWITCH 100

#include <set>
#include <string>
using namespace std;

#include "Protein.hpp"
#include "Conformation.hpp"

class Population {
	public:
		Population( int, Protein, float, float );
		virtual ~Population();

		void setFittest(void);
		Conformation * getFittest(void) const;
		Conformation * rouletteWheelSelect(void);
		void crossover(void);
		void dumpAll(void) const;
		bool isInsertable( const Conformation & );


	private:
		Conformation *individuals;
		Conformation *parent1;
		Conformation *parent2;
		Conformation child1;
		Conformation child2;
		Conformation *theFittest;

		set<string> setOfConformations;
		set<int> collisionSet;
		Protein protein;
		int size;
		float mutProb;
		float crossProb;
};

#endif /* POPULATION_H_ */
