/*
 * Population.cpp
 *
 *  Created on: Apr 6, 2009
 *      Author: David Briemann
 */

#include <cstdlib>
#include <iostream>
#include <set>
using namespace std;

#include "Population.hpp"
#include "Protein.hpp"

Population::Population(int size, Protein prot, float mutProb, float crossProb) {
	this->parent1 = NULL;
	this->parent2 = NULL;
	this->mutProb = mutProb;
	this->crossProb = crossProb;
	this->protein = prot;
	this->size = size;

	//init population
	this->individuals = new Conformation[this->size];
	Conformation temp;

	//COMM
	cout << "Generate Population: " << endl;

	int i = 0;
	while (i < this->size) {
		//create temp conformation
		temp = Conformation(&(this->protein), &(this->collisionSet));
		temp.calcFitness();

		//check if fitness is not 0 (and conformation is not already present in population)
		/*&& this->setOfConformations.insert(temp.getConformationString()).second*/
		if (temp.getFitness() != 0 ) {

			//add to individuals
			this->individuals[i] = temp;
			this->setOfConformations.insert(temp.getConformationString());

			//COMM
			cout << flush;
			cout << i << ".";
			i++;
		}
	}

	this->theFittest = &(this->individuals[0]);
	this->setFittest();
	cout << endl;
}

Population::~Population() {

	this->theFittest = NULL;

	this->parent1 = NULL;

	this->parent1 = NULL;

	//free some space
	delete[] this->individuals;
	this->individuals = NULL;
}

/*
 * checks if a conformation is not yet present in the individuals pool
 * if so returns true
 */
bool Population::isInsertable(const Conformation &candidate) {
	if ((this->setOfConformations.insert(candidate.getConformationString())).second) {
		return true;
	} else {
		return false;
	}
}

/*
 * determines and sets the fittest individual
 */
void Population::setFittest() {
	for (int i = 0; i < this->size; i++) {
		if (this->individuals[i].getFitness() < this->theFittest->getFitness()) {
			this->theFittest = &(this->individuals[i]);
		}
	}
}

/*
 * returns a pointer to the fittest individual
 */
Conformation * Population::getFittest() const {
	return this->theFittest;
}

/*
 * dumps all conformations to console
 */
void Population::dumpAll() const {
	for (int i = 0; i < this->size; i++) {
		cout << i << ": " + this->individuals[i].getStatusString() << endl;
	}
}

/*
 * spins the wheel and selects upon probability, which is determined by fitness
 * returns a pointer to a Conformation
 */
Conformation * Population::rouletteWheelSelect() {
	int totalFitness = 0;
	float randF = 0.0;
	float currentProb = 0.0;
	int i = 0;

	//calculate total fitness
	for (int i = 0; i < this->size; i++) {
		totalFitness += this->individuals[i].getFitness();
	}

	randF = Conformation::randomFloat();
	i = 0;

	currentProb = static_cast<float> (this->individuals[i].getFitness()) / static_cast<float> (totalFitness);

	//be aware of rounding failures? float!
	while ( (currentProb < randF) && (i < (this->size - 1))) {
		currentProb += static_cast<float> (this->individuals[i].getFitness()) / static_cast<float> (totalFitness);
		i++;
	}

	return &(this->individuals[i]);
}

void Population::crossover() {
	float randF = 0.0;

	//find parents
	this->parent1 = this->rouletteWheelSelect();
	randF = Conformation::randomFloat();
	if( this->crossProb < randF ) {
		return; //crossover rate
	}

	this->parent2 = this->rouletteWheelSelect();
	randF = Conformation::randomFloat();
	if( this->crossProb < randF ) {
		return; //crossover rate
	}

	//recombinate parents to children
	this->child1 = Conformation(*(this->parent1), *(this->parent2), &(this->collisionSet));
	this->child2 = Conformation(*(this->parent2), *(this->parent1), &(this->collisionSet));

	//mutate child1
	this->child1.mutate(this->mutProb);
	this->child1.calcValidity();

	if (this->child1.isValid()) { //if child is valid
		//update fitness of child1
		this->child1.calcFitness();

		//if conformation of child1 isnt present yet
		if (this->isInsertable(this->child1)) {
			//replace parent1 if child1 is fitter
			if (this->child1.getFitness() < this->parent1->getFitness()) {
				*(this->parent1) = this->child1;
			} else if (this->child1.getFitness() < this->parent2->getFitness()) { //replace parent2 if child1 is fitter
				*(this->parent2) = this->child1;
			}
		}
	} //else discard

	//mutate child2
	this->child2.mutate(this->mutProb);
	this->child2.calcValidity();

	if (this->child2.isValid()) { //if child 2 is valid
		//update fitness of child2
		this->child2.calcFitness();

		//if conformation of child2 isnt present yet
		if (this->isInsertable(this->child2)) {
			//replace parent2 if child2 is fitter
			if (this->child2.getFitness() < this->parent2->getFitness()) {
				*(this->parent2) = this->child2;
			} else if (this->child2.getFitness() < this->parent1->getFitness()) { //replace parent1 if child2 is fitter
				*(this->parent1) = this->child2;
			}
		}
	} //else discard


	if (this->parent1->getFitness() < this->theFittest->getFitness()) {
		this->theFittest = this->parent1;
		//call visualisation function
	/*	this->theFittest->printStatus();
		this->theFittest->printAsciiPicture();*/
	}

	if (this->parent2->getFitness() < this->theFittest->getFitness()) {
		this->theFittest = this->parent2;
		//call visualisation function
	/*	this->theFittest->printStatus();
		this->theFittest->printAsciiPicture();*/
	}
}
