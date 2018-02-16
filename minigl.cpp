/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
 vec3 color;
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

MGLpoly_mode MGLmode;

struct Vertex {
	vec4 position;
	vec3 color;
};

struct Triangle {
	Vertex A,B,C;
};

vector<Vertex> vertexList;
vector<Triangle> triangles;
/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
}

float triArea(float fia,float fib,float fic,float fja,float fjb,float fjc) {
	return fia*(fjb-fjc) + fja*(fic-fib) +(fib*fjc-(fjb*fic));
};

void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data) {
	float fia,fib,fic,fja,fjb,fjc; 
	fia = ((tri.A.position[0]+1)*width)/2 - 0.5;
	fib = ((tri.B.position[0]+1)*width)/2 - 0.5;
	fic = ((tri.C.position[0]+1)*width)/2 - 0.5;
	fja = ((tri.A.position[1]+1)*width)/2 - 0.5;
	fjb = ((tri.B.position[1]+1)*width)/2 - 0.5;
	fjc = ((tri.C.position[1]+1)*width)/2 - 0.5;
	
	float Area = triArea(fia,fib,fic,fja,fjb,fjc);
	
	
}
/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	MGLmode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(MGLmode == MGL_TRIANGLES){
		for(int i = 0; i < ((vertexList.size()/3)*3); i+=3) {
			Triangle t;
			t.A = vertexList[i];
			t.B = vertexList[i+1];
			t.C = vertexList[i+2];
			triangles.push_back(t);
		} 
		}else if (MGLmode == MGL_QUADS) {
			for(int i = 0; i < ((vertexList.size()/4)*4); i+=4) {
				Triangle t1;
				Triangle t2;
				t1.A = vertexList[i];
				t1.B = vertexList[i+1];
				t1.C = vertexList[i+2];
				t2.A = vertexList[i];
				t2.B = vertexList[i+3];
				t2.C = vertexList[i+2];
				triangles.push_back(t1);
				triangles.push_back(t2);
			}
		
	}
	
	triangles.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	Vertex vertex;
	vertex.position[0] = x;
	vertex.position[1] = y;
	vertex.position[2] = 0;
	vertex.color = color;
	vertexList.push_back(vertex);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	Vertex vertex;
	vertex.position[0] = x;
	vertex.position[1] = y;
	vertex.position[2] = z;
	vertex.position[3] = 1;
	vertex.color = color;
	vertexList.push_back(vertex);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	color[0] = red;
	color[1] = green;
	color[2] = blue;
}
