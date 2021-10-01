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
    unsigned int id;

	Voxel(){}

	Voxel( std::vector<Vec3>  v ){
		vertices = v;
	}

	void evaluate( int id, Vec3 & projectedPoint , Vec3 & projectedNormal ){
        Vec3 temp = ( vertices[id] - projectedPoint );
		values[id] = Vec3::dot( temp , projectedNormal );

	}
};

#endif