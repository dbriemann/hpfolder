/*
 * Protein.h
 *
 *  Created on: Apr 6, 2009
 *      Author: dlb
 */

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include <string>
using namespace std;


class Protein {
	public:
		Protein();
		Protein( string );
		virtual ~Protein();

		char getNth(int) const;
		int getLength(void) const;

	private:
		string HPsequence;
};

#endif /* PROTEIN_H_ */
