// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <math.h>   

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "src/Vec3.h"
#include "src/Grid.h"
#include "src/Voxel.h"
#include "src/Camera.h"
#include "src/jmkdtree.h"

std::vector< Vec3 > positions;
std::vector< Vec3 > normals;

std::vector< Vec3 > positions2;
std::vector< Vec3 > normals2;

std::vector< Vec3 > positions3;
std::vector< Vec3 > normals3;

std::vector< Vec3 > projectedPositions;
std::vector< Vec3 > projectedNormals;

std::vector<Vec3> dcVertices;
std::vector<unsigned int> indexes;
std::vector<Vec3> surfaceVertices;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 640;
static unsigned int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// ------------------------------------------------------------------------------------------------------------
// i/o and some stuff
// ------------------------------------------------------------------------------------------------------------
void loadPN (const std::string & filename , std::vector< Vec3 > & o_positions , std::vector< Vec3 > & o_normals ) {
    unsigned int surfelSize = 6;
    FILE * in = fopen (filename.c_str (), "rb");
    if (in == NULL) {
        std::cout << filename << " is not a valid PN file." << std::endl;
        return;
    }
    size_t READ_BUFFER_SIZE = 1000; // for example...
    float * pn = new float[surfelSize*READ_BUFFER_SIZE];
    o_positions.clear ();
    o_normals.clear ();
    while (!feof (in)) {
        unsigned numOfPoints = fread (pn, 4, surfelSize*READ_BUFFER_SIZE, in);
        for (unsigned int i = 0; i < numOfPoints; i += surfelSize) {
            o_positions.push_back (Vec3 (pn[i], pn[i+1], pn[i+2]));
            o_normals.push_back (Vec3 (pn[i+3], pn[i+4], pn[i+5]));
        }

        if (numOfPoints < surfelSize*READ_BUFFER_SIZE) break;
    }
    fclose (in);
    delete [] pn;
}
void savePN (const std::string & filename , std::vector< Vec3 > const & o_positions , std::vector< Vec3 > const & o_normals ) {
    if ( o_positions.size() != o_normals.size() ) {
        std::cout << "The pointset you are trying to save does not contain the same number of points and normals." << std::endl;
        return;
    }
    FILE * outfile = fopen (filename.c_str (), "wb");
    if (outfile == NULL) {
        std::cout << filename << " is not a valid PN file." << std::endl;
        return;
    }
    for(unsigned int pIt = 0 ; pIt < o_positions.size() ; ++pIt) {
        fwrite (&(o_positions[pIt]) , sizeof(float), 3, outfile);
        fwrite (&(o_normals[pIt]) , sizeof(float), 3, outfile);
    }
    fclose (outfile);
}
void scaleAndCenter( std::vector< Vec3 > & io_positions ) {
    Vec3 bboxMin( FLT_MAX , FLT_MAX , FLT_MAX );
    Vec3 bboxMax( FLT_MIN , FLT_MIN , FLT_MIN );
    for(unsigned int pIt = 0 ; pIt < io_positions.size() ; ++pIt) {
        for( unsigned int coord = 0 ; coord < 3 ; ++coord ) {
            bboxMin[coord] = std::min<float>( bboxMin[coord] , io_positions[pIt][coord] );
            bboxMax[coord] = std::max<float>( bboxMax[coord] , io_positions[pIt][coord] );
        }
    }
    Vec3 bboxCenter = (bboxMin + bboxMax) / 2.f;
    float bboxLongestAxis = std::max<float>( bboxMax[0]-bboxMin[0] , std::max<float>( bboxMax[1]-bboxMin[1] , bboxMax[2]-bboxMin[2] ) );
    for(unsigned int pIt = 0 ; pIt < io_positions.size() ; ++pIt) {
        io_positions[pIt] = (io_positions[pIt] - bboxCenter) / bboxLongestAxis;
    }
}

void applyRandomRigidTransformation( std::vector< Vec3 > & io_positions , std::vector< Vec3 > & io_normals ) {
    srand(time(NULL));
    Mat3 R = Mat3::RandRotation();
    Vec3 t = Vec3::Rand(1.f);
    for(unsigned int pIt = 0 ; pIt < io_positions.size() ; ++pIt) {
        io_positions[pIt] = R * io_positions[pIt] + t;
        io_normals[pIt] = R * io_normals[pIt];
    }
}

void subsample( std::vector< Vec3 > & i_positions , std::vector< Vec3 > & i_normals , float minimumAmount = 0.1f , float maximumAmount = 0.2f ) {
    std::vector< Vec3 > newPos , newNormals;
    std::vector< unsigned int > indices(i_positions.size());
    for( unsigned int i = 0 ; i < indices.size() ; ++i ) indices[i] = i;
    srand(time(NULL));
    std::random_shuffle(indices.begin() , indices.end());
    unsigned int newSize = indices.size() * (minimumAmount + (maximumAmount-minimumAmount)*(float)(rand()) / (float)(RAND_MAX));
    newPos.resize( newSize );
    newNormals.resize( newSize );
    for( unsigned int i = 0 ; i < newPos.size() ; ++i ) {
        newPos[i] = i_positions[ indices[i] ];
        newNormals[i] = i_normals[ indices[i] ];
    }
    i_positions = newPos;
    i_normals = newNormals;
}

bool save( const std::string & filename , std::vector< Vec3 > & vertices , std::vector< unsigned int > & triangles ) {
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open()) {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    unsigned int n_vertices = vertices.size() , n_triangles = triangles.size()/3;
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v ) {
        myfile << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for( unsigned int f = 0 ; f < n_triangles ; ++f ) {
        myfile << 3 << " " << triangles[3*f] << " " << triangles[3*f+1] << " " << triangles[3*f+2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}



// ------------------------------------------------------------------------------------------------------------
// rendering.
// ------------------------------------------------------------------------------------------------------------

void initLight () {
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    glCullFace (GL_BACK);
    glEnable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.f, 0.f, 0.f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
}



void drawTriangleMesh( std::vector< Vec3 > const & i_positions , std::vector< unsigned int > const & i_triangles ) {
    glBegin(GL_TRIANGLES);
    glDisable(GL_CULL_FACE);
    for(unsigned int i = 0 ; i < i_triangles.size(); i += 3) {
        Vec3 p0 = i_positions[i_triangles[i]];
        Vec3 p1 = i_positions[i_triangles[i + 1]];
        Vec3 p2 = i_positions[i_triangles[i + 2]];
        Vec3 n = Vec3::cross(p1-p0 , p2-p0);

        //Offset optionnel
        //p0[0] +=1;
        //p1[0] +=1;
        //p2[0] +=1;

        n.normalize();
        glNormal3f( n[0] , n[1] , n[2] );
        glVertex3f( p0[0] , p0[1] , p0[2] );
        glVertex3f( p1[0] , p1[1] , p1[2] );
        glVertex3f( p2[0] , p2[1] , p2[2] );
    }
    glEnd();
}


void drawPointSet( std::vector< Vec3 > const & i_positions , std::vector< Vec3 > const & i_normals ) {
    glBegin(GL_POINTS);
    for(unsigned int pIt = 0 ; pIt < i_positions.size() ; ++pIt) {
        glNormal3f( i_normals[pIt][0] , i_normals[pIt][1] , i_normals[pIt][2] );
        glVertex3f( i_positions[pIt][0], i_positions[pIt][1] , i_positions[pIt][2] );
    }
    glEnd();
}

void draw () {
    glPointSize(2); // for example...

    glColor3f(0.8,0.8,1);
    drawPointSet(positions , normals);

    glColor3f(0.,0.,1);

//    drawTriangleMesh( dcVertices, indexes );

     glColor3f(1.0,.0,.0);
     drawPointSet(positions3 , normals3);
//
//
    glColor3f(0.,1.0,0.0);
    drawPointSet(surfaceVertices , normals);

//    glColor3f(1.0,.0,0.0);
//    drawPointSet(positions3 , normals);

}

void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    glutPostRedisplay ();
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;

    case 'w':
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        if(polygonMode[0] != GL_FILL)
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        break;

    default:
        break;
    }
    idle ();
}

void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle ();
}

void motion (int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}

void reshape(int w, int h) {
    camera.resize (w, h);
}

void project( Vec3 inputPoint, Vec3 & outputPoint, Vec3 planPoint, Vec3 planNormal ){

    Vec3 pn = inputPoint - planPoint;
    
    outputPoint = inputPoint - ( Vec3::dot(pn, planNormal)/planNormal.length() * planNormal );
}


//void HPSS( Vec3 inputPoint, Vec3 & outputPoint, Vec3 & outputNormal,
//           std::vector<Vec3> const positions, std::vector<Vec3> const normals,
//           BasicANNkdTree const & kdTree, int kernelType= 0, unsigned int nbIterations = 10, unsigned int knn = 10 ){
//
//    ANNidxArray id_nearest_neighbors = new ANNidx[ knn ];
//    ANNdistArray square_distances_to_neighbors = new ANNdist[ knn ];
//
//    // iterate to project input point on centroid plan
//    while ( nbIterations > 0 ){
//        float weitghSum = 0.0;
//        float h = 1.0;
//
//        Vec3 avgPoint = Vec3( 0., 0., 0. );
//        Vec3 avgNormal = Vec3( 0., 0., 0. );
//
//        // get k-nearest and iterate on them to calculate the centroid
//        kdTree.knearest( inputPoint, knn, id_nearest_neighbors, square_distances_to_neighbors );
//
//        for( unsigned int i = 0; i < knn; i++ ){
//
//            Vec3 projectedPoint;
//
//            int index = id_nearest_neighbors[i];
//            // project input point on each neighbor's plan
//            project( inputPoint, projectedPoint, positions[index], normals[index] );
//
//            Vec3 diff = inputPoint - positions[index] ;
//
//            float weight;
//            if ( kernelType == 0 )
//                weight = exp( - (diff.length() * diff.length() ) / ( h * h ) ) ;
//            else if( kernelType == 1 )
//                weight = 1.0;
//            else
//                weight = 1.0;
//
//            avgPoint  += weight * projectedPoint;
//            avgNormal += weight * normals[index];
//            weitghSum += weight;
//        }
//
//        outputPoint  = avgPoint / weitghSum ;
//        outputNormal = avgNormal / weitghSum;
//        outputNormal.normalize();
//        Vec3 temp = Vec3( inputPoint );
//
//        // project input point on centroid's plan
//        project( temp, inputPoint, outputPoint, outputNormal );
//
//        nbIterations--;
//
//    }
//
//    delete [] id_nearest_neighbors;
//    delete [] square_distances_to_neighbors;
//
//}

Vec3 project_orth( const Vec3 & x, const Vec3 & c , const Vec3 & n ){
    Vec3 cx = x-c;
    return x - Vec3::dot(cx, n)*n;
}


void HPSS(Vec3 input_point,
          Vec3 &output_point, Vec3 &output_normal,
          const std::vector<Vec3> &positions, const std::vector<Vec3> &normals,
          const BasicANNkdTree &kd_tree,
          unsigned int nb_iter = 10, uint knn = 10) {

    ANNidxArray id_nearest_neighbors = new ANNidx[knn];
    ANNdistArray square_distances_to_neighbors = new ANNdist[knn];
    Vec3 x_k = input_point, c, n;

    for (unsigned int it = 0; it < nb_iter; it++) {
        kd_tree.knearest(x_k, knn, id_nearest_neighbors, square_distances_to_neighbors);
        double sqrradius = square_distances_to_neighbors[knn-1];
        float total_w = 0.f;
        c = Vec3(0.f, 0.f, 0.f);
        n = Vec3(0.f, 0.f, 0.f);
        for (unsigned int k = 0; k < knn; k++) {
            int j = id_nearest_neighbors[k];
            Vec3 p_tilde = project_orth(x_k, positions[j], normals[j]);
            float w = std::exp(-(x_k - p_tilde).squareLength() / sqrradius);
            total_w += w;
            c += w * p_tilde;
            n += w * normals[j];
        }
        c /= total_w;
        n.normalize();

        x_k = project_orth(x_k, c, n);
    }

    output_point = x_k;
    output_normal = n;

    delete[] id_nearest_neighbors;
    delete[] square_distances_to_neighbors;
}

void AABB( Vec3 &bbMin, Vec3 &bbMax, std::vector<Vec3> const positions ){
    for( unsigned int i =0; i<positions.size();i++){
        Vec3 point = positions[i];

        if( point[0] < bbMin[0] )
            bbMin[0] = point[0];

        else if( point[0] > bbMax[0])
            bbMax[0] = point[0];

        else if( point[1] < bbMin[1])
            bbMin[1] = point[1];

        else if( point[1] > bbMax[1])
            bbMax[1] = point[1];

        else if( point[2] < bbMin[2])
            bbMin[2] = point[2];

        else if( point[2] > bbMax[2])
            bbMax[2] = point[2];
    }

}

    
void dualContouring( std::vector<Vec3> & test, std::vector<Vec3> & test2, std::vector<Vec3> const positions, std::vector<Vec3> const normals, BasicANNkdTree const & kdTree, int kernelType ){
    float gridSize = 32;
    Vec3 bbMin = Vec3(FLT_MAX, FLT_MAX, FLT_MAX), bbMax = Vec3( -FLT_MAX, -FLT_MAX, -FLT_MAX);
    AABB( bbMin, bbMax, positions);
    Grid grid = Grid( gridSize, bbMin, bbMax ) ;


    for ( unsigned int i = 0; i < grid.voxels.size() ; i++ ) {
        Voxel voxel = grid.voxels[i];

        for (unsigned int c = 0; c < 8; c++) {
            Vec3 projectedPoint, projectedNormal;

            HPSS(voxel.vertices[c], projectedPoint, projectedNormal, positions, normals, kdTree );
            voxel.evaluate(c, projectedPoint, projectedNormal);
            positions3.push_back( voxel.vertices[c] );
            normals3.push_back( voxel.vertices[c] );
        }
    }

    for  (Voxel voxel : grid.voxels) {
        Vec3 min = voxel.vertices[0];
        Vec3 max = voxel.vertices[6];
        bool out = false;
        for (int i = 0; i < 7; i++){
            for (int j = 1; i < 8; j++) {
                if (sign(voxel.values[i]) != sign(voxel.values[j])) {
                    surfaceVertices.push_back(Vec3(max + min) / 2.0);
                    out = true;
                    break;
                }
            }
            if (out)
                break;
        }
    }
    std::cout << surfaceVertices.size() << std::endl;

    for (Voxel voxel : grid.voxels ){
        Vec3 min = voxel.vertices[0];
        // X AXIS
        if (min[1]>bbMin[1] && min[2]>bbMin[2]) {
            if ((voxel.values[0]) != (voxel.values[1])){
                Vec3 v1, v2, v3, v4;
                Vec3 n1, n2, n3, n4;

                v1 = surfaceVertices[voxel.id];
                v2 = surfaceVertices[voxel.id  - gridSize];
                v3 = surfaceVertices[voxel.id  - gridSize * gridSize - gridSize];
                v4 = surfaceVertices[voxel.id  - gridSize * gridSize];

//                HPSS(surfaceVertices[voxel.id], v1, n1, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize], v2, n2, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize * gridSize - gridSize], v3, n3, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize * gridSize], v4, n4, positions, normals, kdTree);

                dcVertices.push_back(v1);
                dcVertices.push_back(v2);
                dcVertices.push_back(v3);
                dcVertices.push_back(v4);

                unsigned int count = dcVertices.size();
                indexes.insert( indexes.end(), { count-4, count-3, count-2, count-4, count-2, count-1 });

            }
        }
        // Y AXIS
        if (min[0]>bbMin[0] && min[2]>bbMin[2]) {
            if ((voxel.values[0]) != (voxel.values[3])){
//                Vec3 v1, v2, v3, v4;
//                Vec3 n1, n2, n3, n4;

                Vec3 v1 = surfaceVertices[voxel.id];
                Vec3 v2 = surfaceVertices[voxel.id - 1];
                Vec3 v3 = surfaceVertices[voxel.id - gridSize * gridSize - 1];
                Vec3 v4 = surfaceVertices[voxel.id - gridSize * gridSize];

//                HPSS(surfaceVertices[voxel.id], v1, n1, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - 1], v2, n2, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize * gridSize - 1], v3, n3, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize * gridSize], v4, n4, positions, normals, kdTree);

                dcVertices.push_back(v1);
                dcVertices.push_back(v2);
                dcVertices.push_back(v3);
                dcVertices.push_back(v4);

                unsigned int count = dcVertices.size();
                indexes.insert( indexes.end(), { count-4, count-3, count-2, count-4, count-2, count-1 });

            }
        }
        // Z AXIS
        if (min[0]>bbMin[0] && min[1]>bbMin[1]) {
            if ((voxel.values[0]) != (voxel.values[4])){
//                Vec3 v1, v2, v3, v4;
//                Vec3 n1, n2, n3, n4;
//
                Vec3 v1 = surfaceVertices[voxel.id];
                Vec3 v2 = surfaceVertices[voxel.id - 1];
                Vec3 v3 = surfaceVertices[voxel.id - gridSize -1];
                Vec3 v4 = surfaceVertices[voxel.id - gridSize];

//                HPSS(surfaceVertices[voxel.id], v1, n1, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - 1], v2, n2, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize - 1], v3, n3, positions, normals, kdTree);
//                HPSS(surfaceVertices[voxel.id  - gridSize], v4, n4, positions, normals, kdTree);

                dcVertices.push_back(v1);
                dcVertices.push_back(v2);
                dcVertices.push_back(v3);
                dcVertices.push_back(v4);

                unsigned int count = dcVertices.size();
                indexes.insert( indexes.end(), { count-4, count-3, count-2, count-4, count-2, count-1 });

            }
        }

    }
    std::cout << "end"<< std::endl;
}

int main (int argc, char ** argv) {
    if (argc > 2) {
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("tp point processing");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);


    {
        // Load a first pointset, and build a kd-tree:
        loadPN("pointsets/igea.pn" , positions , normals);

        BasicANNkdTree kdtree;
        kdtree.build(positions);

        // Create a second pointset that is artificial, and project it on pointset1 using MLS techniques:

        dualContouring(positions2, normals2, positions, normals, kdtree, 0);
        save("tp.off", dcVertices, indexes);
    }

    glutMainLoop ();
    return EXIT_SUCCESS;
}
