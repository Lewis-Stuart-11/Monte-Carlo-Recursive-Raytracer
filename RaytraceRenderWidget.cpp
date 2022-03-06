//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  Raytrace Render Widget
//  -----------------------------
//
//	Provides a widget that displays a fixed image
//	Assumes that the image will be edited (somehow) when Render() is called
//  
////////////////////////////////////////////////////////////////////////

#include <math.h>

// include the header file
#include "RaytraceRenderWidget.h"

// constructor
RaytraceRenderWidget::RaytraceRenderWidget
        (   
        // the geometric object to show
        TexturedObject      *newTexturedObject,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // parent widget in visual hierarchy
        QWidget             *parent
        )
    // the : indicates variable instantiation rather than arbitrary code
    // it is considered good style to use it where possible
    : 
    // start by calling inherited constructor with parent widget's pointer
    QOpenGLWidget(parent),
    // then store the pointers that were passed in
    texturedObject(newTexturedObject),
    renderParameters(newRenderParameters)
    { // constructor
    // leaves nothing to put into the constructor body
    } // constructor    

// destructor
RaytraceRenderWidget::~RaytraceRenderWidget()
    { // destructor
    // empty (for now)
    // all of our pointers are to data owned by another class
    // so we have no responsibility for destruction
    // and OpenGL cleanup is taken care of by Qt
    } // destructor                                                                 

// called when OpenGL context is set up
void RaytraceRenderWidget::initializeGL()
    { // RaytraceRenderWidget::initializeGL()
	// this should remain empty
    } // RaytraceRenderWidget::initializeGL()

// called every time the widget is resized
void RaytraceRenderWidget::resizeGL(int w, int h)
    { // RaytraceRenderWidget::resizeGL()
    // resize the render image
    frameBuffer.Resize(w, h);

    } // RaytraceRenderWidget::resizeGL()
    
// called every time the widget needs painting
void RaytraceRenderWidget::paintGL()
    { // RaytraceRenderWidget::paintGL()
    // set background colour to white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    // If lighting is enabled, then configure all light sources in the scene (default is 1)
    if(renderParameters -> useLighting)
        ConfigureLuminairs();

    // Convert all vertex and normal positions to NDCS for correct intersection testing 
    ConvertWCSToNDCSForScene();
    
    // Perform ray tracing
    Raytrace();

    // Display the image
    glDrawPixels(frameBuffer.width, frameBuffer.height, GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer.block);
    } // RaytraceRenderWidget::paintGL()

// Generates an image by shooting a series of rays through the image pixels into the scene and calculating
// the light to use for each pixel
void RaytraceRenderWidget::Raytrace()
    { // RaytraceRenderWidget::Raytrace()

    // Sets the default extinction and minimum albedo values, these are only used if recursive scattering
    // is enabled. Feel free to change these if required
    if(renderParameters->useRecursiveScattering){
        renderParameters->extinctionFactor = 0.05f;
        renderParameters->albedoThreshold = 0.03f;
    }

    Cartesian3 eyePosition, pixelPosition;

    Ray newRay;

    RGBAValue pixelColor; 

    float NDCSRow, NDCSCol;

    // Default background colour to use if no faces intercept with the current ray
    float background[4] = {0.8f, 0.8f, 0.6f, 1.0f};

    // Loops through each pixel on the screen and shoots a ray through it
    for (int row = 0; row < frameBuffer.height; row++){
        for (int col = 0; col < frameBuffer.width; col++){
            
            // Sets the total light to 0
            float totalLight[4] = {0.0f, 0.0f, 0.0f, 0.0f};

            // Converts the pixel positions from DCS to NDCS, this is so that 
            // the ray is in the correct coordinate system as the rest of the vertices
            NDCSRow = ((float)row / (float)(frameBuffer.height/2)) - 1.0f;
            NDCSCol = ((float)col / (float)(frameBuffer.width/2)) - 1.0f;

            // Creates Cartesian3 position for the pixel
            pixelPosition = Cartesian3(NDCSCol, NDCSRow, 0);

            // This assumes that orthogonal projection is used, as rays are shot 
            // 'straight through' the current pixel in the Z direction. Persepctice 
            // projection could used by setting a default position for the eye for 
            // every pixel
            eyePosition = Cartesian3(NDCSCol, NDCSRow, -2.0f);

            newRay.pointOrigin = eyePosition;

            // Ray has direction (0,0,2)
            newRay.directionVector = pixelPosition-eyePosition;

            // Default albedo is 1.0f, as the ray is new and has not bounced
            newRay.albedo = 1.0f;

            // Loops through the number of samples
            for(int sample=0; sample < renderParameters->numberOfSamples; sample++){

                // Ray is traced through the scene and the calculated lighting values
                // are returned
                float *returnedLightingIntensity = TracePath(newRay);

                // In the trace ray method, if no intersections are found, then a flag is 
                // returned in the 5 index of the lighting array. This indicates that 
                // the default background colour should be used. The reason TracePath 
                // simply does not return the background colour if no ray is found is
                // because this would result in lighting errors when using scattering
                // (no intersections should return 0 light)
                if(*(returnedLightingIntensity+4) == 1.0f)
                    for(int i=0; i<4; i++)
                        totalLight[i] += background[i];
                
                // If intersections are found, then add to total light
                else
                    for(int i=0; i<4; i++)
                        totalLight[i] += *(returnedLightingIntensity+i);
                
            }

            // Convert lighting array to RGBA value and divide by the number of samples
            pixelColor.red = (totalLight[0]*255)/(float)renderParameters->numberOfSamples;
            pixelColor.green = (totalLight[1]*255)/(float)renderParameters->numberOfSamples;
            pixelColor.blue = (totalLight[2]*255)/(float)renderParameters->numberOfSamples;
            pixelColor.alpha = (totalLight[3]*255)/(float)renderParameters->numberOfSamples;

            // Add new pixel to frame buffer
            frameBuffer.block[(row*frameBuffer.width) + col] = pixelColor;
        }
    }
} // RaytraceRenderWidget::Raytrace()

// Traces the ray through the scene and returns the lighting value if intersections are made
float * RaytraceRenderWidget::TracePath(Ray ray){
    float textureColourAlbedo[4], lightingColour[4];

    RGBAValue textureColour;

    // Sets the 5th array value to 1.0f to indicate that the background colour should be used
    float backgroundColour[5] = {0.0f, 0.0f, 0.0f, 0.0f, 1.0f};

    // Default colour to use if lighting is not enabled
    float defaultColour[4] =  {0.7f, 0.7f, 0.7f, 1.0f};

    // Stores the colour values to return
    static float finalColour[5];

    finalColour[4] = 0.0f;

    // If the ray albedo is too low, then return no light
    if(ray.albedo < renderParameters->albedoThreshold){
        for(int i=0; i<4; i++)
            finalColour[i] = backgroundColour[i];
        return finalColour;
    }

    for(int i=0; i<4; i++)
        finalColour[i] = defaultColour[i];

    // Get the closest triangle along the ray, the maximum distance from the ray to find 
    // a face is set to infinite (no bounding box)
    TriangleWithPoint closestTriangle = GetClosestTriangle(ray, std::numeric_limits<float>::max());

    // If the face Id of the returned face is -1, then it means that no intersections have been found, and the 
    // background colour should be returned
    if(closestTriangle.FaceID == -1){
        for(int i=0; i<5; i++)
            finalColour[i] = backgroundColour[i];
        return finalColour;
    }

    // Surface properties are found for that point on the face
    FaceProperties faceProperty = GetSurfaceElementFromTriangle(closestTriangle);

    // If lighting is enabled
    if(renderParameters -> useLighting){

        // Set default lighting values to the face emission properties
        float lightingColorIntensities[4] = {faceProperty.faceMaterial.Emission[0], faceProperty.faceMaterial.Emission[1], 
                                            faceProperty.faceMaterial.Emission[2], faceProperty.faceMaterial.Emission[3]};

        float *returnedLightingIntensity;

        // Loop through each light source in the scene
        for(int light=0; light < allLighting.size(); light++){
            
            // Get the light values for each light source and add to the total light
            returnedLightingIntensity = CalculateDirectLight(faceProperty, (-1)*ray.directionVector, allLighting[light]);

            for(int i=0; i<4; i++)
                lightingColorIntensities[i] += *(returnedLightingIntensity + i);

        }

        // Get the indirect light values and add to total light
        returnedLightingIntensity = CalculateIndirectLight(faceProperty, (-1)*ray.directionVector, allLighting[0], ray.albedo);

        for(int i=0; i<4; i++)
            lightingColorIntensities[i] += *(returnedLightingIntensity + i);
        
        // Make sure that light values are not larger than 1 and add to final light colour
        for(int i=0; i<4; i++){
            lightingColour[i] = std::min(lightingColorIntensities[i], (float) 1.);
            finalColour[i] = lightingColour[i];
        }
    }

    // If rendering is enabled
    if(renderParameters->texturedRendering){ 

        int texelColumn, texelRow, texturePixelLookup;

        // If the scene assets should be used (values stored in .scene file) then use these textures
        if(renderParameters->useSceneAssets){
            if(isnan(faceProperty.textureImageCoordS))
                faceProperty.textureImageCoordS = 0.0f;

            if(isnan(faceProperty.textureImageCoordT))
                faceProperty.textureImageCoordT = 0.0f;

            // If texture ID is equal to -1, then this face should not be textured
            if(faceProperty.textureID == -1)
                for(int i=0; i<4; i++)
                    textureColourAlbedo[i] = defaultColour[i];
            else{
                // Get texture row and column values in texture image
                texelColumn = (int)round(faceProperty.textureImageCoordT * texturedObject->textures[faceProperty.textureID].height);
                texelRow = (int)round(faceProperty.textureImageCoordS * texturedObject->textures[faceProperty.textureID].width);

                texturePixelLookup = (texelColumn * texturedObject->textures[faceProperty.textureID].width) + texelRow;

                // Get texture colour
                textureColour = texturedObject->textures[faceProperty.textureID].block[texturePixelLookup];

                // Texture RGBA value is converted into lighting array
                textureColourAlbedo[0] = (float)textureColour.red/255;
                textureColourAlbedo[1] = (float)textureColour.green/255;
                textureColourAlbedo[2] = (float)textureColour.blue/255;
                textureColourAlbedo[3] = (float)textureColour.alpha/255;
            }
        }

        // Uses the loaded object texture
        else{

            // Get texture row and column values in texture image
            texelColumn = (int)round(faceProperty.textureImageCoordT * texturedObject->texture.height);
            texelRow = (int)round(faceProperty.textureImageCoordS * texturedObject->texture.width);

            texturePixelLookup = (texelColumn * texturedObject->texture.width) + texelRow;

            // Get texture colour
            textureColour = texturedObject->texture.block[texturePixelLookup];

            // Texture RGBA value is converted into lighting array
            textureColourAlbedo[0] = (float)textureColour.red/255;
            textureColourAlbedo[1] = (float)textureColour.green/255;
            textureColourAlbedo[2] = (float)textureColour.blue/255;
            textureColourAlbedo[3] = (float)textureColour.alpha/255;
        }
        
        // If texture modulation is used, then combine lighting and texturing
        if(renderParameters -> textureModulation && renderParameters -> useLighting)
            for(int i=0; i<4; i++)
                finalColour[i] = textureColourAlbedo[i] * lightingColour[i];
        
        // Otherwise just use the texture values for the final colour
        else
            for(int i=0; i<4; i++)
                finalColour[i] = textureColourAlbedo[i];
    }

    // Return calculated final colour
    return finalColour;
}


// Calculates the direct light for a face in the scene
float * RaytraceRenderWidget::CalculateDirectLight(FaceProperties face, Cartesian3 LightOut, Luminaire Light){

    // Default light is set to 0
    static float totalDirectLight[4];
    for(int i = 0; i < 4; i++)
        totalDirectLight[i] = 0.0f;

    Cartesian3 LightPos(Light.LightPosition[0], Light.LightPosition[1], Light.LightPosition[2]);

    // The light direction vector is calcualted as the direction from the position on the face to the current light source
    Cartesian3 LightIn = LightPos - face.PointPosition;

    // If shadows are used, then shoot a ray from the position on the face towards the light, if an intersection
    // is made then the light is blocked, and return 0 for the light. The maximum length to traverse along the ray
    // is the distance from the light to the position on the face, this is so that if an face is behind the light
    // then it will not count as blocking the light
    if(renderParameters->useShadows){
        Ray shadowCheckRay;

        shadowCheckRay.pointOrigin = face.PointPosition;
        shadowCheckRay.directionVector = LightIn;

        if(GetClosestTriangle(shadowCheckRay, LightIn.length()).FaceID != -1)
            return totalDirectLight;
    }

    // Calculate the loss in lighting intensity as the distance squared from the light to the position on the face
    float distanceSquared;
    if(Light.isInfinite)
        distanceSquared = 1;
    else
        distanceSquared = LightIn.dot(LightIn); 

    // Get the total BRDF and add the lighting values to total light 
    // (divided by the distance squared and multiplied by light colour)
    float *totalBRDF = GetBRDF(face, LightOut, LightIn);
    for(int i = 0; i < 4; i++)
        totalDirectLight[i] = *(totalBRDF + i) * renderParameters->lightColor[i] / distanceSquared;

    return totalDirectLight;
}

// Gets the indirect light for a face in the scene
float * RaytraceRenderWidget::CalculateIndirectLight(FaceProperties face, Cartesian3 LightOut, Luminaire Light, float currentRayAlbedo){

    static float totalLight[5];

    // If recursive scattering is not used, then return the simple Ambient material calculation
    if(!renderParameters->useRecursiveScattering){
        for(int i = 0; i < 4; i++)
            totalLight[i] = (float)renderParameters->lightColor[i] * face.faceMaterial.Ambient[i];
        return totalLight;
    }

    for(int i = 0; i < 5; i++)
        totalLight[i] = 0.0f; 

    float *albedoValues, *inLight; 
    float maxAlbedo = 0.0f;

    // Chooses a random value between 0-1 (with 3 decimal places)- if the value is lower than the extinction 
    // factor then return no light (kill ray) 
    if((rand() % 1000)/1000.0f < renderParameters->extinctionFactor)
        return totalLight;

    // If reflections are enabled, then choose a random value between 0-1 (with 3 decimal places). If this
    // is less than the impulse reflection chance, then use reflections
    bool useReflections = false;
    if(renderParameters->useReflections)
        if((rand() % 1000)/1000.0f < face.impulseChance)
            useReflections = true;

    Cartesian3 LightIn;

    if(useReflections){
        
        // Use perfect reflection as the ray to use for indirect light
        LightIn = (2.0f * face.Normal) - LightOut;

        static float reflectiveAlbedoValues[4];

        for(int i=0; i < 4; i++)
            reflectiveAlbedoValues[i] = face.impulseAlbedo;

        // Get the albedo for a refelection
        albedoValues = reflectiveAlbedoValues;
    }

    else{
        // Get the ray to use as indirect light from the Monte Carlo vector generator
        LightIn = GetMonteCarloVector(face.Normal);

        // Get the BRDF values for that vector
        albedoValues = GetBRDF(face, LightOut, LightIn);
    }
    
    // Find the largest value (either r,g or b) for that light, and use that as the largest albedo light value 
    // returned
    for(int i=0; i < 3; i++)
         if(*(albedoValues + i) > maxAlbedo)
            maxAlbedo = *(albedoValues + i);


    // Create a new ray and shoot it from the face in the direction of the monte carlo vector to 
    // get the incoming indirect light
    Ray newRay;
    newRay.pointOrigin = face.PointPosition;
    newRay.directionVector = LightIn;
    newRay.albedo = maxAlbedo * currentRayAlbedo;

    inLight = TracePath(newRay);

    // Use this as the total incoming light and multiply by the BRDF
    for(int i=0; i<4; i++){
        totalLight[i] = (float) *(inLight+i) * *(albedoValues+i);
    }

    // Makes sure that background flag is not set to 1 (nasty bug)
    inLight[4] = 0.0f;

    return totalLight;
}

// Calculates the BRDF (Specular and diffuse) for a point on the face
float * RaytraceRenderWidget::GetBRDF(FaceProperties face, Cartesian3 LightOut, Cartesian3 LightIn){

    // Default light is set to 0
    static float totalBRDF[4];
    for(int i=0; i < 4; i++)
        totalBRDF[i] = 0.0f;

    float specularLightOfSurface[4], diffuseLightOfSurface[4], specularReflectionWithShine, diffuseReflection;

    // Reflection vector is found as the sum of the vector from the light to point, and the vector
    // from the pixel to the point
    Cartesian3 reflectionVector = LightIn + LightOut;

    // Specular value is found with shineness
    specularReflectionWithShine = pow(std::max(reflectionVector.unit().dot(face.Normal.unit()), (float)0.), face.faceMaterial.Shineness);

    // Diffuse values are calculated
    diffuseReflection = LightIn.unit().dot(face.Normal.unit());

    // Add these values to total light, multiplying with the material values
    for(int i=0; i < 4; i++){
        totalBRDF[i] += face.faceMaterial.Specular[i] * specularReflectionWithShine;
        totalBRDF[i] += face.faceMaterial.Diffuse[i] * std::max(diffuseReflection, (float)0.);
        totalBRDF[i] = std::min(totalBRDF[i], (float) 1.);
    }

    return totalBRDF;
}

// Gets a random reflection vector to use in scattering
Cartesian3 RaytraceRenderWidget::GetMonteCarloVector(Cartesian3 normal){

    Cartesian3 newVector;
    int accuracy = 2000;
    float pi = 3.1415;
    float angle;

    // Loops until a valid vector is generated
    while(true){

        // Generates random x, y and z values for the direction vector between -1 and 1
        newVector.x = (float)((rand() % accuracy) - (accuracy/2))/accuracy;
        newVector.y = (float)((rand() % accuracy) - (accuracy/2))/accuracy;
        newVector.z = (float)((rand() % accuracy) - (accuracy/2))/accuracy;

        // Calculates the angle between the new vector and the normal
        angle = acos(normal.dot(newVector)/(normal.length() * newVector.length()));

        // If the angle is between 90 and 270 degrees, then the vector goes into the face,
        // which does not make sense (unless material is translucent)
        if(angle > pi/2 && angle < pi*(3/2))
            continue;

        // Checks that length is valid for an acceptable reflection vector to be returned
        // (returned vectors should not show bias in diagonals)
        if(newVector.length() > 1.0f || newVector.length() < 0.1f)
            continue;
        

        // Returns unit vector if the vector is valid (angle and length are acceptable)
        return newVector.unit();
    }

    return Cartesian3(0.0f, 0.0f, 0.0f);
}


// Gets the closest triangle that intersects with the current ray
TriangleWithPoint RaytraceRenderWidget::GetClosestTriangle(Ray R, float cutOffDistanceAlongRay){

    // The value of T indicates how far the face is along the ray from the starting position.
    // This is used to measure which face is closest to the origin position of the ray (smaller 
    // value of T indicates that face is closer than current face), hence this is set to infinity
    float currentSmallestTDistance = std::numeric_limits<float>::max();

    // Creates a defualt face with FaceID = -1. This indicates that no faces were found
    TriangleWithPoint clostestTriangle;
    clostestTriangle.FaceID = -1;

    // If object is not shown then return
    if (!renderParameters->showObject)
        return clostestTriangle;

    // Value ensures that faces behind the ray origin are not considered (T must be positive)
    // and reduces chances of self intersection. This is an alternative to shifting the ray along
    // a face normal when calculating shadows, etc...
    float e = R.directionVector.length()*0.002;

    // ##### The next loop is copied from the Texture object script for reading faces #####

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
        
    // If object scaling is requested, apply it as well 
    if (renderParameters->scaleObject)
        scale /= texturedObject->objectSize;

    // Loops through each face
    for (int face = 0; face < texturedObject->faceVertices.size(); face++)
        { 

        for (unsigned int triangle = 0; triangle < texturedObject->faceVertices[face].size() - 2; triangle++)
            { 
            
            std::vector<Cartesian3> currentVertices;

            // Loops through each vertex on the face and adds to the current list of vertices.
            // Vertex positions are taken in the NDCS coordinate system (calculated earlier)
            for (unsigned int vertex = 0; vertex < 3; vertex++)
                { 

                int faceVertex = 0;

                if (vertex != 0)
                    faceVertex = triangle + vertex;

                Cartesian3 currentVertex;

                currentVertex.x = scale * NDCSVertices[ texturedObject->faceVertices   [face][faceVertex]].x;
                currentVertex.y = scale * NDCSVertices[ texturedObject->faceVertices   [face][faceVertex]].y;
                currentVertex.z = scale * NDCSVertices[ texturedObject->faceVertices   [face][faceVertex]].z;
                
                currentVertices.push_back(currentVertex);

            } 

            // Points
            Cartesian3 pointO, pointP, pointQ, pointR;

            pointP = currentVertices[0];
            pointQ = currentVertices[1];
            pointR = currentVertices[2];

            // The following code performs the intersection calculation specified in Lecture 14

            // Vectors
            Cartesian3 vectU, vectV, vectN, vectW, vectO;

            // Prime points
            Cartesian3 primePointO, primePointP, primePointQ, primePointR;

            // Plane points
            Cartesian3 planePointO, planePointP, planePointQ, planePointR;

            // Calculates vectors for basis
            vectU = pointQ - pointP;
            vectV = pointR - pointP;

            vectU = vectU.unit();

            vectN = (vectU.cross(vectV)).unit();
            vectW = (vectN.cross(vectU)).unit();

            vectO = R.pointOrigin - pointP;

            // Calculates the point along the ray that intersects with the basis plane
            float scalarT = -(vectO.dot(vectN))/(R.directionVector.dot(vectN));

            pointO = R.pointOrigin + (R.directionVector*scalarT);

            // Point is converted to plane basis
            primePointO = pointO - pointP;
            primePointP = pointP - pointP;
            primePointQ = pointQ - pointP;
            primePointR = pointR - pointP;

            planePointO.x = primePointO.dot(vectU);
            planePointO.y = primePointO.dot(vectW);
            planePointO.z = primePointO.dot(vectN);

            planePointP.x = primePointP.dot(vectU);
            planePointP.y = primePointP.dot(vectW);
            planePointP.z = primePointP.dot(vectN);

            planePointQ.x = primePointQ.dot(vectU);
            planePointQ.y = primePointQ.dot(vectW);
            planePointQ.z = primePointQ.dot(vectN);

            planePointR.x = primePointR.dot(vectU);
            planePointR.y = primePointR.dot(vectW);
            planePointR.z = primePointR.dot(vectN);

            Cartesian3 planePointPQ, planePointQR, planePointRP;

            planePointPQ = planePointQ - planePointP;
            planePointQR = planePointR - planePointQ;
            planePointRP = planePointP - planePointR;

            Cartesian3 normalPQ, normalQR, normalRP;

            normalPQ.x = -planePointPQ.y;
            normalPQ.y = planePointPQ.x;
            normalPQ.z = 0.0;
            
            normalQR.x = -planePointQR.y;
            normalQR.y = planePointQR.x;
            normalQR.z = 0.0;

            normalRP.x = -planePointRP.y;
            normalRP.y = planePointRP.x;
            normalRP.z = 0.0;

            // Performs half-plane test and find alpha, beta and gamma for point on the face
            float lineConstantPQ, lineConstantQR, lineConstantRP;

            lineConstantPQ = normalPQ.dot(planePointP);
            lineConstantQR = normalQR.dot(planePointQ);
            lineConstantRP = normalRP.dot(planePointR);

            float distanceP, distanceQ, distanceR;

            distanceP = normalQR.dot(planePointP) - lineConstantQR;
            distanceQ = normalRP.dot(planePointQ) - lineConstantRP;
            distanceR = normalPQ.dot(planePointR) - lineConstantPQ;

            float alpha, beta, gamma;

            alpha = (normalQR.dot(planePointO) - lineConstantQR) / distanceP;
            beta = (normalRP.dot(planePointO) - lineConstantRP) / distanceQ;
            gamma = (normalPQ.dot(planePointO) - lineConstantPQ) / distanceR;

            // Half-plane test
            if ((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0))
                continue;
            
            // Checks if the current face is closest to the ray origin by comparing scalar T values
            if(scalarT >= currentSmallestTDistance)
                continue;

            // Face must be along the positive direction of the ray
            if(scalarT < e)
                continue;

            // Checks the distance of the point along the ray is not larger than the length cut-off (used for shadows)
            if((R.directionVector*scalarT).length() > cutOffDistanceAlongRay)
                continue;

            // New smallest value of T is updated
            currentSmallestTDistance = scalarT;

            // A new triangle object is created which stores the face ID. 
            TriangleWithPoint newTriangle;

            newTriangle.TriangleID = triangle;
            newTriangle.FaceID = face;

            // Alpha, beta and gamma are also stored as these will be need for barycentric interpolation later when 
            // retrieving other values (such as normal along the face)
            newTriangle.alpha = alpha;
            newTriangle.beta = beta;
            newTriangle.gamma = gamma;

            // Point is stored as this is need when calculating lighting etc...
            newTriangle.PointPosition = pointO;

            // Triangle to return is set to the new triangle
            clostestTriangle = newTriangle;
            
        } 
    } 

    // Returns the current triangle, if no intersections are found, then this will return the default triangle
    return clostestTriangle;
}

// Gets all required properties for the point on the face using interpolation
FaceProperties RaytraceRenderWidget::GetSurfaceElementFromTriangle(TriangleWithPoint triangle){

    std::vector<float> triangleUCoord;
    std::vector<float> triangleVCoord;
    std::vector<float *> triangleColour;
    std::vector<Cartesian3> triangleNormals;

    int face = triangle.FaceID;
    float alpha = triangle.alpha;
    float beta = triangle.beta;
    float gamma = triangle.gamma;

    // Loops through each vertex
    for (unsigned int vertex = 0; vertex < 3; vertex++){
        int faceVertex = 0;

        if (vertex != 0)
            faceVertex = triangle.TriangleID + vertex;

        // Gets the UV coordinates for the vertex on the face
        triangleUCoord.push_back(texturedObject->textureCoords   [texturedObject->faceTexCoords  [face][faceVertex]  ].x);
        triangleVCoord.push_back(texturedObject->textureCoords   [texturedObject->faceTexCoords  [face][faceVertex]  ].y);
        
        // Gets the normal for the vertex on the face
        Cartesian3 currentNormal;

        currentNormal.x = NDCSNormals     [texturedObject->faceNormals    [face][faceVertex]  ].x;
        currentNormal.y = NDCSNormals     [texturedObject->faceNormals    [face][faceVertex]  ].y;
        currentNormal.z = NDCSNormals     [texturedObject->faceNormals    [face][faceVertex]  ].z;

        triangleNormals.push_back(currentNormal);

        // If UV is mapped to RGB, then add this colour to the triangle colour
        if(renderParameters->mapUVWToRGB){
            float *colourPointer = (float *) &(texturedObject->textureCoords[texturedObject->faceTexCoords[face][faceVertex]]);
            triangleColour.push_back(colourPointer);
        }
        
    }

    // New face property to stor these values
    FaceProperties newFaceProperty;

    newFaceProperty.FaceID = face;
    newFaceProperty.PointPosition = triangle.PointPosition;

    // Gets the reflection chance and albedo if reflections are enabled
    if(renderParameters->useReflections){
        newFaceProperty.impulseChance = texturedObject->faceImpulseReflectionsChance[face];
        newFaceProperty.impulseAlbedo = texturedObject->faceImpulseReflectionAlbedo[face];
    }
    
    // Interpolates UV coordinates along triangle
    newFaceProperty.textureImageCoordS = alpha * triangleUCoord[0] + beta * triangleUCoord[1] + gamma * triangleUCoord[2];
    newFaceProperty.textureImageCoordT = alpha * triangleVCoord[0] + beta * triangleVCoord[1] + gamma * triangleVCoord[2];
    
    // Interpolates normal along triangle
    newFaceProperty.Normal.x =  alpha * triangleNormals[0].x + beta * triangleNormals[1].x + gamma * triangleNormals[2].x;
    newFaceProperty.Normal.y =  alpha * triangleNormals[0].y + beta * triangleNormals[1].y + gamma * triangleNormals[2].y;
    newFaceProperty.Normal.z =  alpha * triangleNormals[0].z + beta * triangleNormals[1].z + gamma * triangleNormals[2].z;
    
    // If mapping to UV is used, then get the interpolated material values and set the face colour to the 
    // calculated RGB values
    if(renderParameters->mapUVWToRGB){

        for(int i=0; i<3; i++){
            newFaceProperty.faceMaterial.Ambient[i] = alpha * triangleColour[0][i] + beta * triangleColour[1][i]  + gamma * triangleColour[2][i];
            newFaceProperty.faceMaterial.Diffuse[i] = alpha * triangleColour[0][i] + beta * triangleColour[1][i]  + gamma * triangleColour[2][i];
            newFaceProperty.faceMaterial.Specular[i] = alpha * triangleColour[0][i] + beta * triangleColour[1][i]  + gamma * triangleColour[2][i];
            newFaceProperty.faceMaterial.Emission[i] = renderParameters->emissive;
        }

        newFaceProperty.faceMaterial.Ambient[3] = 1.0f; 
        newFaceProperty.faceMaterial.Diffuse[3] = 1.0f; 
        newFaceProperty.faceMaterial.Emission[3] = 1.0f; 
        newFaceProperty.faceMaterial.Specular[3] = 1.0f; 

        newFaceProperty.RGBColour.red = (alpha * triangleColour[0][0] + beta * triangleColour[1][0]  + gamma * triangleColour[2][0])*255;
        newFaceProperty.RGBColour.green = (alpha * triangleColour[0][1] + beta * triangleColour[1][1]  + gamma * triangleColour[2][1])*255;
        newFaceProperty.RGBColour.blue = (alpha * triangleColour[0][2] + beta * triangleColour[1][2]  + gamma * triangleColour[2][2])*255;

        newFaceProperty.faceMaterial.Shineness = renderParameters->specularExponent;

    }

    // If UV to RGB is not used
    else{

        // If the scene assets are used, then the material is set to the material stored in the textured object
        // which is read in from the scene file
        if(renderParameters->useSceneAssets)
            newFaceProperty.faceMaterial = texturedObject->faceMaterials[face];

        // Otherwise, if an object is rendered, use the material values defined by the window sliders
        else{
            for(int i=0; i<3; i++){
                newFaceProperty.faceMaterial.Ambient[i] = renderParameters->ambient; 
                newFaceProperty.faceMaterial.Diffuse[i] = renderParameters->diffuse;
                newFaceProperty.faceMaterial.Emission[i] = renderParameters->emissive;
                newFaceProperty.faceMaterial.Specular[i] = renderParameters->specular; 
            }

            newFaceProperty.faceMaterial.Ambient[3] = 1.0f; 
            newFaceProperty.faceMaterial.Diffuse[3] = 1.0f; 
            newFaceProperty.faceMaterial.Emission[3] = 1.0f; 
            newFaceProperty.faceMaterial.Specular[3] = 1.0f;

            newFaceProperty.faceMaterial.Shineness = renderParameters->specularExponent; 
        } 

        // Sets object colour to the default value
        newFaceProperty.RGBColour = RGBAValue(0.7f*255, 0.7f*255, 0.7f*255, 255.f);
    }

    // If scene assets are used, set the texture to use when rendering as the texture ID stored for that face
    if(renderParameters->useSceneAssets)
        newFaceProperty.textureID = texturedObject->faceTextures[face];

    return newFaceProperty;
}

// Configures lighting objects in the scene
void RaytraceRenderWidget::ConfigureLuminairs(){

    Luminaire currentLighting;

    allLighting.clear();

    // Configures lighting matrix to be in NDCS and the correct position
    ConfigureLightingMatrix();

    // All lighting is not infinite
    currentLighting.isInfinite = false;

    float coords[4];

    // If the scene is used, then the coordinates for the lighting position in the scene file is used, otherwise
    // set this to the default position defined by the render parameters
    if(!renderParameters->useSceneAssets){
        coords[0] = renderParameters->lightPosition[0];
        coords[1] = renderParameters->lightPosition[1];
        coords[2] = renderParameters->lightPosition[2];
    }
    else{
        coords[0] = texturedObject->lightPosition.x;
        coords[1] = texturedObject->lightPosition.y;
        coords[2] = texturedObject->lightPosition.z;
    }

    coords[3] = 1.0f;
    
    // Calculates the new lighting position (multiplies position with lighting matrix)
    for(int i=0; i < 3; i++){
        float newCoord = 0;
        for(int j=0; j < 4; j++){
            newCoord += lightingMatrix[i][j] * coords[j];
        }
        currentLighting.LightPosition[i] = newCoord;
        currentLighting.lightColor[i] = renderParameters->lightColor[i];
    }

    currentLighting.LightPosition[3] = renderParameters->lightPosition[3];
    currentLighting.lightColor[3] = renderParameters->lightColor[3];
    
    // Adds this to all lighting
    allLighting.push_back(currentLighting);
}

// Converts vertices and normals to NDCS   
void RaytraceRenderWidget::ConvertWCSToNDCSForScene(){

    // Configures the vertex and normal matrices
    ConfigureVertexMatrix();
    ConfigureNormalMatrix();

    NDCSVertices.resize(texturedObject->vertices.size());
    NDCSNormals.resize(texturedObject->normals.size());

    // Transforms vertices and normals
    for(int vertexID=0; vertexID < texturedObject->vertices.size(); vertexID++)
        NDCSVertices[vertexID] = TransformElement(texturedObject->vertices[vertexID], vertexMatrix);

    for(int normalID=0; normalID < texturedObject->normals.size(); normalID++)
        NDCSNormals[normalID] = TransformElement(texturedObject->normals[normalID], normalMatrix);
    
}

// Configures the vertex matrix to transform vertices to NDCS. These use the same calls that 
// the render widget uses with OpenGL 
void RaytraceRenderWidget::ConfigureVertexMatrix(){

    // Creates a new identity matrix
    std::array< std::array<float,4>, 4> newIdentityMatrix{{
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}
    }};

    vertexMatrix = newIdentityMatrix;

    // translate by the visual translation
    Translatef(renderParameters->xTranslate, renderParameters->yTranslate, 0.0f, vertexMatrix);

    // apply rotation matrix from arcball
    MultMatrixf(renderParameters->rotationMatrix.columnMajor().coordinates, vertexMatrix);

    // compute the aspect ratio of the widget
    float aspectRatio = (float) frameBuffer.width / (float) frameBuffer.height;
    
    // we want to capture a sphere of radius 1.0 without distortion
    // so we set the ortho projection based on whether the window is portrait (> 1.0) or landscape
    // portrait ratio is wider, so make bottom & top -1.0 & 1.0
    if (aspectRatio > 1.0)
        Ortho(-aspectRatio, aspectRatio, -1.0, 1.0, -1.0, 1.0, vertexMatrix);
    // otherwise, make left & right -1.0 & 1.0
    else
        Ortho(-1.0, 1.0, -1.0/aspectRatio, 1.0/aspectRatio, -1.0, 1.0, vertexMatrix);

    if (renderParameters->showObject){
        
        // Scale defaults to the zoom setting
        float scale = renderParameters->zoomScale;
        
        // if object scaling is requested, apply it as well 
        if (renderParameters->scaleObject)
            scale /= texturedObject->objectSize;
            
        // apply the translation to the centre of the object if requested
        if (renderParameters->centreObject)
            Translatef(-texturedObject->centreOfGravity.x * scale, -texturedObject->centreOfGravity.y * scale, -texturedObject->centreOfGravity.z * scale, vertexMatrix);
    }
}

// Configures lighting matrix (using similar calls to OpenGl)
void RaytraceRenderWidget::ConfigureLightingMatrix(){
    std::array< std::array<float,4>, 4> newIdentityMatrix{{
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}
    }};

    lightingMatrix = newIdentityMatrix;

    MultMatrixf(renderParameters->lightMatrix.columnMajor().coordinates, lightingMatrix);

    Translatef(renderParameters->xTranslate, renderParameters->yTranslate, 0.0f, lightingMatrix);

}

// Configures normal matrix (using similar calls to OpenGl) 
void RaytraceRenderWidget::ConfigureNormalMatrix(){
    std::array< std::array<float,4>, 4> newIdentityMatrix{{
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}
    }};

    normalMatrix = newIdentityMatrix;

    MultMatrixf(renderParameters->rotationMatrix.columnMajor().coordinates, normalMatrix);
}

/*  ###### The following methods are copied from the previous coursework and have the purpose of  ######
    ###### converting all vertices in the scene to DCS so that raytracing can occur correctly. ######*/

Cartesian3 RaytraceRenderWidget::TransformElement(Cartesian3 element, std::array< std::array<float,4>, 4> &matrix){

    Homogeneous4 homoElement;
    
    homoElement.x = element.x;
    homoElement.y = element.y;
    homoElement.z = element.z;
    homoElement.w = 1.;

    // Converts OCS to VCS (position is multiplied by the view-model matrix)
    Homogeneous4 currentPosition = multiplyMatrixAndHomogeneoursCoords(homoElement, matrix);

    Cartesian3 newCartesianCoords;

    // Divide by W to get NDCS
    newCartesianCoords.x = currentPosition.x / currentPosition.w;
    newCartesianCoords.y = currentPosition.y / currentPosition.w;
    newCartesianCoords.z = currentPosition.z / currentPosition.w;

    return newCartesianCoords;
}


// Multiplies a matrix with a homogeneous coordinate
Homogeneous4 RaytraceRenderWidget::multiplyMatrixAndHomogeneoursCoords(Homogeneous4 position, std::array< std::array<float,4>, 4> &matrix){

    // Converts Homogeneous4 object to array for easy traversal
    float coords[4] = {position.x, position.y, position.z, position.w};

    float newCoords[4];

    // Each coordinate is multiplied by each entry in the respective matrix row
    for(int i=0; i < 4; i++){
        float newCoord = 0;
        for(int j=0; j < 4; j++){
            newCoord += matrix[i][j] * coords[j];
        }
        newCoords[i] = newCoord;
    }

    // Position is updated with the new coordinates and returned
    position.x = newCoords[0];
    position.y = newCoords[1];
    position.z = newCoords[2];
    position.w = newCoords[3];

    return position;
}

// Multiplies two 4x4 matrices together and returns the product
void RaytraceRenderWidget::multiplyMatrices(std::array< std::array<float,4>, 4> matrix1, std::array< std::array<float,4>, 4> &matrix2){
    std::array< std::array<float,4>, 4> multipliedMatrix{{
        {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}
    }};

    // Performs a 4x4 matrix multiplication between the two arguments
    for(int i=0; i < 4; i++){
        for(int j=0; j < 4; j++){
            for(int k=0; k < 4; k++){
                multipliedMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    matrix2 = multipliedMatrix;
}


// translate the matrix
void RaytraceRenderWidget::Translatef(float xTranslate, float yTranslate, float zTranslate, std::array< std::array<float,4>, 4> &matrix){ // Translatef()

    // Sets translation matrix values
    std::array< std::array<float,4>, 4> translateMatrix{{
        {1,0,0,xTranslate}, {0,1,0,yTranslate}, {0,0,1,zTranslate}, {0,0,0,1}
    }};

    multiplyMatrices(translateMatrix, matrix);
} // Translatef()


// multiply by a known matrix in column-major format
void RaytraceRenderWidget::MultMatrixf(const float *columnMajorCoordinates, std::array< std::array<float,4>, 4> &matrix){ // MultMatrixf()

    std::array< std::array<float,4>, 4> formattedMatrix;

    // Converts column major coordinate matrix to a 4x4 float vector (which is the used format for matrices)
    for(int i=0; i < 4; i++){
        for(int j=0; j < 4; j++){
            formattedMatrix[i][j] = columnMajorCoordinates[(j*4) + i];
        }
    }

    multiplyMatrices(formattedMatrix, matrix);
} // MultMatrixf()


// sets an orthographic projection matrix
void RaytraceRenderWidget::Ortho(float left, float right, float bottom, float top, float zNear, float zFar, std::array< std::array<float,4>, 4> &matrix){ // Ortho()

    // Sets default ortho matrix values
    std::array< std::array<float,4>, 4> orthoMatrix{{
        {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,1}
    }};

    // Sets ortho matrix values according to the OpenGL specifications
    orthoMatrix[0][0] = 2 / (right - left);
    orthoMatrix[0][3] = -((right + left) / (right - left));

    orthoMatrix[1][1] = 2 / (top - bottom);
    orthoMatrix[1][3] = -((top + bottom) / (top - bottom));

    orthoMatrix[2][2] = -(2 / (zFar - zNear));
    orthoMatrix[2][3] = -((zFar + zNear) / (zFar - zNear));

    multiplyMatrices(orthoMatrix, matrix);
} // Ortho()


// mouse-handling
void RaytraceRenderWidget::mousePressEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mousePressEvent()
    // store the button for future reference
    int whichButton = event->button();
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0 * event->x() - size) / size;
    float y = (size - 2.0 * event->y() ) / size;

    
    // and we want to force mouse buttons to allow shift-click to be the same as right-click
    int modifiers = event->modifiers();
    
    // shift-click (any) counts as right click
    if (modifiers & Qt::ShiftModifier)
        whichButton = Qt::RightButton;
    
    // send signal to the controller for detailed processing
    emit BeginScaledDrag(whichButton, x,y);
    } // RaytraceRenderWidget::mousePressEvent()
    
void RaytraceRenderWidget::mouseMoveEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseMoveEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0 * event->x() - size) / size;
    float y = (size - 2.0 * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit ContinueScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseMoveEvent()
    
void RaytraceRenderWidget::mouseReleaseEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseReleaseEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0 * event->x() - size) / size;
    float y = (size - 2.0 * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit EndScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseReleaseEvent()
