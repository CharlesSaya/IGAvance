#ifndef GRID_H
#define GRID_H

#include <cmath>
#include <iostream>

#include "src/Vec3.h"
#include "src/Voxel.h"

class Grid{

private : 


public:
	std::vector<Voxel> voxels;

	
	Grid(){}

	Grid( float gridSize, Vec3 bbMin, Vec3 bbMax ){

		float stepX = ( bbMax[0] - bbMin[0] )/( gridSize -1.0);
		float stepY = ( bbMax[1] - bbMin[1] )/( gridSize -1.0);
		float stepZ = ( bbMax[2] - bbMin[2] )/( gridSize -1.0);

	    for ( int z = 0; z < (int)(gridSize) - 1; z++){
	        for ( int y = 0; y < (int)(gridSize) - 1; y++){
	            for ( int x = 0; x < (int)(gridSize) - 1; x++){
	            	std::vector<Vec3> corners;
	            	Vec3 point = Vec3( bbMin[0] + (float)x * stepX, bbMin[1] + (float)y * stepY, bbMin[2] + (float)z * stepZ );

	            	getCorners( corners, point, stepX, stepY, stepZ);
	            	Voxel voxel = Voxel(corners);

	            	voxels.push_back( voxel ); 

	            }
	        }
	    }

	}

	void getCorners( std::vector<Vec3> & corners, Vec3 &point, float stepX, float stepY, float stepZ ){

		corners.push_back( point );
		corners.push_back( Vec3( point[0] + stepX, point[1] , point[2] ) );
		corners.push_back( Vec3( point[0], point[1] + stepY, point[2] ) );
		corners.push_back( Vec3( point[0] + stepX, point[1] + stepY, point[2] ) );

		corners.push_back( Vec3( point[0] , point[1], point[2] + stepZ ) );
		corners.push_back( Vec3( point[0] + stepX, point[1] , point[2] + stepZ ) );
		corners.push_back( Vec3( point[0], point[1] + stepY, point[2] + stepZ ) );
		corners.push_back( Vec3( point[0] + stepX, point[1] + stepY, point[2] + stepZ ) );

	}


};

#endif