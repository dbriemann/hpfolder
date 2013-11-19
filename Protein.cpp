/*
 * Protein.cpp
 *
 *  Created on: Apr 6, 2009
 *      Author: dlb
 */

#include "Protein.hpp"
#include <string>
using namespace std;


//default constructor
Protein::Protein() {
}

Protein::Protein( string seq ) {
	this->HPsequence = seq;
}

Protein::~Protein() {
}

char Protein::getNth(int i) const {
	return this->HPsequence[i];
}

int Protein::getLength() const {
	return this->HPsequence.length();
}
