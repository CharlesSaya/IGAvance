#ifndef VOXEL_H
#define VOXEL_H

#include <cmath>
#include <iostream>

#include "src/Vec3.h"

class Voxel{

private :
	

public:
	std::vector<Vec3> vertices;
	float values[8];
	float id;

	Voxel(){}

	Voxel( std::vector<Vec3>  v ){
		vertices = v;
	}

	void evaluate( int id, Vec3 const & projectedPoint , Vec3 const & projectedNormal ){
		values[id] = Vec3::dot( ( vertices[id] - projectedPoint ), projectedNormal );
	
	}
};

#endif