/*
 * PointCollisionSet.h
 *
 *  Created on: May 5, 2009
 *      Author: dlb
 */

#ifndef POINTCOLLISIONSET_H_
#define POINTCOLLISIONSET_H_

class PointCollisionSet {
	public:
		PointCollisionSet();
		virtual ~PointCollisionSet();

	private:
		bool insert();
		void clear();

};

#endif /* POINTCOLLISIONSET_H_ */
