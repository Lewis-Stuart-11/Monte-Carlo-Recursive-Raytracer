//////////////////////////////////////////////////////////////////////
//
//	University of Leeds
//	COMP 5812M Foundations of Modelling & Rendering
//	User Interface for Coursework
//
//	September, 2020
//
//  -----------------------------
//  Raytrace Render Widget
//  -----------------------------
//
//	Provides a widget that displays a fixed image
//	Assumes that the image will be edited (somehow) when Render() is called
//  
////////////////////////////////////////////////////////////////////////

// include guard
#ifndef _RAYTRACE_RENDER_WIDGET_H
#define _RAYTRACE_RENDER_WIDGET_H

// include the relevant QT headers
#include <QOpenGLWidget>
#include <QMouseEvent>

// and include all of our own headers that we need
#include "TexturedObject.h"
#include "RenderParameters.h"
#include "Cartesian3.h"

#include <limits>
#include <cmath>
#include <math.h>
#include <stdlib.h>

// A ray that is used for retrieving light values in a scene
class Ray
    { 
    public:

	// Origin point of the ray
    Cartesian3 pointOrigin;

	// Direction ray is traveling in
    Cartesian3 directionVector;
    
	// The lighting intensity of the ray (used to determine whether to discard ray)
	float albedo;
	}; 


// Stores all properties needed for a point on face
class FaceProperties
	{

	public:

	// ID of the face in question
	int FaceID;

	// Point on the face
	Cartesian3 PointPosition;

	// The face normal vector
	Cartesian3 Normal;

	// The texture ID to use for texturing
	int textureID;

	// UV coordinates
	float textureImageCoordS;
	float textureImageCoordT;

	// Material values
	Material faceMaterial;

	// Face colour
	RGBAValue RGBColour;

	// Chance of a reflection 
	float impulseChance;

	// Albedo for a relection (light intensity)
	float impulseAlbedo;
	
	};

// Stores a point on a triangle face 
// (smaller calss than Face Properties but basically performs the same function)
class TriangleWithPoint
	{
	public:

	// Face and triangle ID 
	int FaceID;
	int TriangleID;

	// Used for barycentric interpolation
	float alpha;
	float beta;
	float gamma;
	
	// Position of point on the triangle
	Cartesian3 PointPosition;

	};


// A light source in the scene
class Luminaire
	{

	public:

	// Light for position and colour
	float LightPosition[4];
	float lightColor[4];

	// Stores whether the light source is infinite
	bool isInfinite;

	};

// class Material is defined in TextureObject.cpp

// class for a render widget with arcball linked to an external arcball widget
class RaytraceRenderWidget : public QOpenGLWidget										
	{ // class RaytraceRenderWidget
	Q_OBJECT
	private:	
	// the geometric object to be rendered
	TexturedObject *texturedObject;	

	// the render parameters to use
	RenderParameters *renderParameters;

	// An image to use as a framebuffer
	RGBAImage frameBuffer;

	// Matrices for converting WCS to NDCS
	std::array< std::array<float,4>, 4> vertexMatrix;

	std::array< std::array<float,4>, 4> normalMatrix;

	std::array< std::array<float,4>, 4> lightingMatrix;

	// NDCS coordinates for scene vertices and normals
	std::vector<Cartesian3> NDCSVertices;

	std::vector<Cartesian3> NDCSNormals;

	// Stores all light sources in the scene
	std::vector<Luminaire> allLighting;

	// ##### Ray handling methods ####

	// Routine that generates the image
    void Raytrace();

	// Traces the ray through the scene and returns the lighting value if intersections are made
	float * TracePath(Ray ray);

	// ###### Lighting Methods ######

	// Gets the indirect light for a face in the scene
	float *CalculateIndirectLight(FaceProperties face, Cartesian3 LightOut, Luminaire Light, float currentRayAlbedo);

	// Calculates the direct light for a face in the scene
	float *CalculateDirectLight(FaceProperties face, Cartesian3 LightOut, Luminaire Light);

	// Calculates the BRDF (Specular and diffuse) for a point on the face
	float *GetBRDF(FaceProperties face, Cartesian3 LightOut, Cartesian3 LightIn);

	// Gets a random reflection vector to use in scattering
	Cartesian3 GetMonteCarloVector(Cartesian3 normal);

	// Configures lighting objects in the scene
	void ConfigureLuminairs();

	// ##### Face calculation methods #####

	// Gets the closest triangle that intersects with the current ray
	TriangleWithPoint GetClosestTriangle(Ray R, float cutOffDistanceAlongRay);

	// Gets all required properties for the point on the face using interpolation
	FaceProperties GetSurfaceElementFromTriangle(TriangleWithPoint triangle);

	// ##### Coordinate handling methods #####

	// Converts vertices and normals to NDCS   
	void ConvertWCSToNDCSForScene();

	// Handles matrices to convert to NDCS
	void ConfigureVertexMatrix();

	void ConfigureNormalMatrix();

	void ConfigureLightingMatrix();

	// Converts element to NDCS
	Cartesian3 TransformElement(Cartesian3 element, std::array< std::array<float,4>, 4> &matrix);

	// Multiplies homogenrous coordiates with matrix
    Homogeneous4 multiplyMatrixAndHomogeneoursCoords(Homogeneous4 position, std::array< std::array<float,4>, 4> &matrix);

	// ##### Modified OpenGL methods ######

    void multiplyMatrices(std::array< std::array<float,4>, 4> matrix1, std::array< std::array<float,4>, 4> &matrix2);

    void Translatef(float xTranslate, float yTranslate, float zTranslate, std::array< std::array<float,4>, 4> &matrix);

    void MultMatrixf(const float *columnMajorCoordinates, std::array< std::array<float,4>, 4> &matrix);

    void Ortho(float left, float right, float bottom, float top, float zNear, float zFar, std::array< std::array<float,4>, 4> &matrix);

	public:
	// constructor
	RaytraceRenderWidget
			(
	 		// the geometric object to show
			TexturedObject 		*newTexturedObject,
			// the render parameters to use
			RenderParameters 	*newRenderParameters,
			// parent widget in visual hierarchy
			QWidget 			*parent
			);
	
	// destructor
	~RaytraceRenderWidget();
			
	protected:

	// called when OpenGL context is set up
	void initializeGL();
	// called every time the widget is resized
	void resizeGL(int w, int h);
	// called every time the widget needs painting
	void paintGL();

	// mouse-handling
	virtual void mousePressEvent(QMouseEvent *event);
	virtual void mouseMoveEvent(QMouseEvent *event);
	virtual void mouseReleaseEvent(QMouseEvent *event);

	// these signals are needed to support shared arcball control
	public:
	signals:
	// these are general purpose signals, which scale the drag to 
	// the notional unit sphere and pass it to the controller for handling
	void BeginScaledDrag(int whichButton, float x, float y);
	// note that Continue & End assume the button has already been set
	void ContinueScaledDrag(float x, float y);
	void EndScaledDrag(float x, float y);
	}; // class RaytraceRenderWidget

#endif
