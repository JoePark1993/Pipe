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



MGLmatrix_mode currentMatrixMode;
mat4 currentMatrix = {{1,0,0,0,
					   0,1,0,0,
					   0,0,1,0,
					   0,0,0,1}};
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
vector<mat4> projectionMat(1);
vector<mat4> modelviewMat(1);
/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}



float triArea(float fia,float fib,float fic,float fja,float fjb,float fjc) {
	return fia*(fjb-fjc) + fja*(fic-fib) +((fib*fjc)-(fjb*fic));
};

mat4& top_of_stack(){
	if (currentMatrixMode == MGL_PROJECTION) {
		return projectionMat.back();
	} else {
		return modelviewMat.back();
	}
};

void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data) {
	float fia,fib,fic,fja,fjb,fjc; 
	fia = ((tri.A.position[0]/tri.A.position[3]+1)*width)/2 - 0.5;
	fib = ((tri.B.position[0]/tri.B.position[3]+1)*width)/2 - 0.5;
	fic = ((tri.C.position[0]/tri.C.position[3]+1)*width)/2 - 0.5;
	fja = ((tri.A.position[1]/tri.A.position[3]+1)*height)/2 - 0.5;
	fjb = ((tri.B.position[1]/tri.B.position[3]+1)*height)/2 - 0.5;
	fjc = ((tri.C.position[1]/tri.C.position[3]+1)*height)/2 - 0.5;
	float Area = triArea(fia,fib,fic,fja,fjb,fjc);
  	cout << Area << endl ;	
	for(int i = 0; i < width; i++) {
 	  for(int j = 0; j< height; j++) {
           // float fi = ((i+1)*width)/2 -0.5;
           // float fj = ((j+1)*height)/2 - 0.5;
            float alpha = triArea(i,fib,fic,j,fjb,fjc)/Area;
	    float beta = triArea(fia,i,fic,fja,j,fjc)/Area;
            float gamma = triArea(fia,fib,i,fja,fjb,j)/Area;
            if(alpha >=0.0 && alpha <= 1.0 && beta >= 0.0 && beta <= 1.0 && gamma <= 1 && gamma >= 0) {
              data[i+j*width] = Make_Pixel(tri.A.color[0]*255,tri.B.color[1]*255,tri.C.color[2]*255);
            }
          } 
	}
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
  for(unsigned int i = 0; i < width; i++) {
    for(unsigned int j = 0; j < height; j++) {
      data[i+j*width] = Make_Pixel(0,0,0);
    }
  }
  for(auto tri : triangles) {
     Rasterize_Triangle(tri,width,height,data); 
  }
  triangles.clear();
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
int size = ((vertexList.size()/3)*3);
int size2 = ((vertexList.size()/4)*4);
	if(MGLmode == MGL_TRIANGLES){
		for(int i = 0; i < size; i+=3) {
			Triangle t;
			t.A = vertexList[i];
			t.B = vertexList[i+1];
			t.C = vertexList[i+2];
			triangles.push_back(t);
		} 
	} else if (MGLmode == MGL_QUADS) {
			for(int i = 0; i < size2; i+=4) {
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
	
	vertexList.clear();
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
        vertex.position[3] = 1;
	vertex.color = color;
	vertex.position = projectionMat.back() * modelviewMat.back() * vertex.position;
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
    vertex.position = projectionMat.back() * modelviewMat.back() * vertex.position;
	vertexList.push_back(vertex);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  currentMatrixMode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	mat4& ref = top_of_stack();
	if(currentMatrixMode == MGL_PROJECTION) {
		projectionMat.push_back(ref);
	} else if(currentMatrixMode == MGL_MODELVIEW) {
		modelviewMat.push_back(ref);
	}

	
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(currentMatrixMode == MGL_PROJECTION) {
		projectionMat.pop_back();
	} else if(currentMatrixMode == MGL_MODELVIEW) {
		modelviewMat.pop_back();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	mat4& ref = top_of_stack();
	ref = {{1,0,0,0,
		    0,1,0,0,
		    0,0,1,0,
		    0,0,0,1}};
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
		mat4 matrix1;
	for ( int i = 0; i < 16; i++) {
		matrix1.values[i] = matrix[i];
	}
	mat4& curr = top_of_stack();
	curr = matrix1;
	
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
	mat4 matrix1;
	for ( int i = 0; i < 16; i++) {
		matrix1.values[i] = matrix[i];
	}
	mat4& curr = top_of_stack();
	curr = curr * matrix1;

}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	mat4 trans = {{1,0,0,0,
				   0,1,0,0,
				   0,0,1,0,
				   x,y,z,1}};
	mat4& curr = top_of_stack();
	curr = curr * trans;
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
	MGLfloat rad = (angle*M_PI)/180;
	MGLfloat c = cos(rad);
	MGLfloat s = sin(rad);
	vec3 joseph;
	joseph[0] = x;
	joseph[1] = y;
	joseph[2] = z;
	joseph = joseph.normalized();
	x = joseph[0];
	y = joseph[1];
	z = joseph[2];
	MGLfloat e = (1-c);
	mat4 rot = {{x*x*e+c, y*x*e + z*s, x*z*e - y*s, 0,
			     x*y*e-z*s,y*y*e+c, y*z*e + x*s,0,
			     x*z*e+y*s,y*z*e-x*s,z*z*e+c,0,
			     0,0,0,1}};
			     
	mat4& curr = top_of_stack();
	curr = curr * rot;
	 
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat4 scale = {{x,0,0,0,
				   0,y,0,0,
				   0,0,z,0,
				   0,0,0,1}};
	mat4& curr = top_of_stack();
	curr = curr * scale;
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

  float A = (right + left)/(right - left);
  float B = (top + bottom)/(top - bottom);
  float C = -(far + near)/(far-near);
  float D = -(2*far*near)/(far-near);

  float a = (2*near)/(right-left);
  float b = (2*near)/(top-bottom);
  mat4 frust = {{ a,0.0f,0.0f,0.0f,0.0f,b,0.0f,0.0f,A,B,C,-1.0f,0.0f,0.0f,D,0.0f}};
  mat4& test = top_of_stack();
  test = frust * test;
 

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
  float tx = -(right+left)/(right-left);
  float ty = -(top+bottom)/(top-bottom);
  float tz = -(far+near)/(far-near);
  float a = (2/(right-left));
  float b = (2/(top-bottom));
  float c = -(2/(far-near));
  mat4 ortho = {{a,0.0f,0.0f,0.0f,0.0f,b,0.0f,0.0f,0.0f,0.0f,c,0.0f,tx,ty,tz,1.0f}};
  mat4& test = top_of_stack();
 test = ortho * test;


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
