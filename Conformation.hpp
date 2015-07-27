/*
 * Conformation.h
 *
 *  Created on: Apr 6, 2009
 *      Author: David Briemann
 */

#ifndef CONFORMATION_H_
#define CONFORMATION_H_

#include <set>
#include <string>
using namespace std;

#include "Protein.hpp"

//relative folding directions
#define LEFT -1
#define RIGHT 1
#define FORWARD 0

class Conformation {
	public:
		Conformation();
		//Conformation(const Conformation &);
		Conformation(const Conformation &, const Conformation &, set<int> *);
		Conformation(Protein *, set<int> *);
		virtual ~Conformation();

		const Conformation & operator=(const Conformation &);
		bool operator==(const Conformation &) const;
		bool operator!=(const Conformation &) const;

		static int point(short, short);
		static short extractX(int);
		static short extractY(int);
		static float randomFloat(void);

		static unsigned int energyEvalSteps;

		void calcAbsolutePositions(void);
		void generateRandomConformation(bool);
		void mutate(float);
		void calcValidity(void);
		string getConformationString(void) const;
		void calcFitness(void);
		int getFitness(void) const;
		int getLength(void) const;
		bool isValid(void) const;
		string getStatusString(void) const;
		int getGeneration(void) const;
		void olden(void);
		Protein * getProtein(void) const;
		void printAsciiPicture(void) const;
		int getAbsAt(int) const;

private:
		int generation;
		int fitness;
		Protein *protein;
		int length;
		char *encoding;
		int *absPositions;
		bool validState;
		set<int> *setOfPoints;
};

#endif /* CONFORMATION_H_ */
