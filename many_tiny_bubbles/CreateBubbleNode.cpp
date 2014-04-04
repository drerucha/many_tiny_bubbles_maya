#include "CreateBubbleNode.h"
#define MNoVersionString
#define MNoPluginEntry
#include <maya/MFnPlugin.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Windows.h>
#include <atlstr.h>
#include <maya/MFnStringData.h>
#include <time.h>
#include <maya/MStatus.h>
#include <maya/MObject.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MPlug.h>
#include <maya/MVector.h>
#include <maya/MDoubleArray.h>
#include <maya/MVectorArray.h>
#include <maya/MPlugArray.h>
#include <maya/MPointArray.h>

#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnContainerNode.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>

#include <maya/MDynamicsUtil.h>
#include <maya/MGlobal.h>
#include <sstream>
#include "mac_grid.h"
#include "constants.h"


extern MString mllPath;

MStatus returnStatus;

#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
	}


//MObject CreateBubbleNode::containerSizeX;
//MObject CreateBubbleNode::containerSizeY;
//MObject CreateBubbleNode::containerSizeZ;
//MObject CreateBubbleNode::containerResolutionX;
//MObject CreateBubbleNode::containerResolutionY;
//MObject CreateBubbleNode::containerResolutionZ;
MObject CreateBubbleNode::viscosity;
MObject CreateBubbleNode::density;
MObject CreateBubbleNode::scatterFreq;
MObject CreateBubbleNode::scatterCoef;
MObject CreateBubbleNode::bubbleSize;
MObject CreateBubbleNode::time;

MObject CreateBubbleNode::outputMesh;
MTypeId CreateBubbleNode::id( 0x00128 );



MACGrid mGrid;

void* CreateBubbleNode::creator()
{
	return new CreateBubbleNode;
}


MStatus CreateBubbleNode::initialize()
{


	//Constants::setContainerSize(9,5,1);////////////////////////
	//mGrid.reset();
	//MFnNumericAttribute defaultContainerSizeX;
	//MFnNumericAttribute defaultContainerSizeY;
	//MFnNumericAttribute defaultContainerSizeZ;
	//MFnNumericAttribute defaultContainerResolutionX;
	//MFnNumericAttribute defaultContainerResolutionY;
	//MFnNumericAttribute defaultContainerResolutionZ;
	MFnNumericAttribute defaultViscosity;
	MFnNumericAttribute defaultDensity;
	MFnNumericAttribute defaultScatterFreq;
	MFnNumericAttribute defaultScatterCoef;
	MFnNumericAttribute defaultBubbleSize;
	MFnUnitAttribute defaultTime;
	MFnTypedAttribute defaultGeometry;



	CreateBubbleNode::viscosity = defaultViscosity.create( "viscosity", "v", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultViscosity.setMax(10.0);
	defaultViscosity.setMin(0.1);

	CreateBubbleNode::density = defaultDensity.create( "density", "d", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultDensity.setMax(10.0);
	defaultDensity.setMin(0.1);

	CreateBubbleNode::scatterFreq = defaultScatterFreq.create( "scaterfreq", "f", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultScatterFreq.setMax(10.0);
	defaultScatterFreq.setMin(0.1);

	CreateBubbleNode::scatterCoef = defaultScatterCoef.create( "scatercoef", "c", MFnNumericData::kDouble, 0.8, &returnStatus );
	defaultScatterCoef.setMax(10.0);
	defaultScatterCoef.setMin(0.1);

	CreateBubbleNode::bubbleSize = defaultBubbleSize.create( "bubblesize", "bs", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultBubbleSize.setMax(10.0);
	defaultBubbleSize.setMin(0.1);


	CreateBubbleNode::time = defaultTime.create( "time", "tm", MFnUnitAttribute::kTime, 2.0, &returnStatus );

	CreateBubbleNode::outputMesh = defaultGeometry.create( "outputMesh", "out", MFnData::kMesh, &returnStatus ); 
	//defaultGeometry.setStorable(false);
	//defaultGeometry.setStorable(false);



	//adding the attribute
	//returnStatus = addAttribute(CreateBubbleNode::containerSizeX);
	//returnStatus = addAttribute(CreateBubbleNode::containerSizeY);
	//returnStatus = addAttribute(CreateBubbleNode::containerSizeZ);
	//returnStatus = addAttribute(CreateBubbleNode::containerResolutionX);
	//returnStatus = addAttribute(CreateBubbleNode::containerResolutionY);
	//returnStatus = addAttribute(CreateBubbleNode::containerResolutionZ);
	returnStatus = addAttribute(CreateBubbleNode::viscosity);
	returnStatus = addAttribute(CreateBubbleNode::density);
	returnStatus = addAttribute(CreateBubbleNode::scatterFreq);
	returnStatus = addAttribute(CreateBubbleNode::scatterCoef);
	returnStatus = addAttribute(CreateBubbleNode::bubbleSize);
	returnStatus = addAttribute(CreateBubbleNode::time);
	returnStatus = addAttribute(CreateBubbleNode::outputMesh);

	//assigning Dependencies
	returnStatus = attributeAffects(CreateBubbleNode::time,        CreateBubbleNode::outputMesh);
	returnStatus = attributeAffects(CreateBubbleNode::viscosity,   CreateBubbleNode::outputMesh);
	returnStatus = attributeAffects(CreateBubbleNode::density,     CreateBubbleNode::outputMesh);
	returnStatus = attributeAffects(CreateBubbleNode::scatterFreq, CreateBubbleNode::outputMesh);
	returnStatus = attributeAffects(CreateBubbleNode::scatterCoef, CreateBubbleNode::outputMesh);
	returnStatus = attributeAffects(CreateBubbleNode::bubbleSize,  CreateBubbleNode::outputMesh);



	return MS::kSuccess;
}


MStatus CreateBubbleNode::compute(const MPlug& plug, MDataBlock& data)
{
	MStatus returnStatus;
	
	if (plug == outputMesh) 
	{
		MDataHandle timeData = data.inputValue(time, &returnStatus ); 
		MTime times = timeData.asTime();


		MDataHandle viscosityData = data.inputValue(viscosity, &returnStatus );
		double viscositys = viscosityData.asDouble();
		MDataHandle densityData = data.inputValue(density, &returnStatus );
		double densitys = densityData.asDouble();
		MDataHandle scatterFreqData = data.inputValue(scatterFreq, &returnStatus );
		double scatterFreqs = scatterFreqData.asDouble();
		MDataHandle scatterCoefData = data.inputValue(scatterCoef, &returnStatus );
		double scatterCoefs = scatterCoefData.asDouble();
		MDataHandle bubbleSizeData = data.inputValue(bubbleSize, &returnStatus );
		double bubbleSizes = 0.05;// set bubble size

		//currently, this part is called in each frame, but it should be called only when the container size or resolution is changed 
		double containSizeX,containSizeY,containSizeZ;
		int resolutionX, resolutionY, resolutionZ;
		MIntArray resolutionArray;
		MGlobal::executeCommand("getAttr fluidShape1.dimensionsW", containSizeX);
		MGlobal::executeCommand("getAttr fluidShape1.dimensionsH", containSizeY);
		MGlobal::executeCommand("getAttr fluidShape1.dimensionsD", containSizeZ);
		MGlobal::executeCommand("getAttr fluidShape1.resolution", resolutionArray);
		resolutionX = resolutionArray[0];
		resolutionY = resolutionArray[1];
		resolutionZ = resolutionArray[2];

		double gridSize = containSizeX / resolutionX; 


		Constants::setContainerDim(resolutionX, resolutionY, resolutionZ);//resolution cannot less than 3 
		Constants::setGridSize(gridSize);
		Constants::setBubbleRadius(bubbleSizes);



		MDataHandle outputHandle = data.outputValue(outputMesh, &returnStatus);
		MFnMeshData dataCreator;
		MObject newOutputData = dataCreator.create(&returnStatus);

		returnStatus = createBubble(times, newOutputData, containSizeX, containSizeY, containSizeZ, viscositys, densitys, scatterFreqs, scatterCoefs, bubbleSizes, plug, data );
		outputHandle.set(newOutputData);
		data.setClean( plug );


	} 
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}


MStatus CreateBubbleNode::createBubble(const MTime &time, MObject& outData, double &containerSizeX, double &containerSizeY, double &containerSizeZ, double &viscosity, double &density, double &scatterFreq, double &scatterCoef, double &bubbleSize, const MPlug& plug, MDataBlock& block )
{
	int bubbleRadiusCount = 5;
	double bubbleRadiusMin = 0.01;
	double bubbleRadiusMax = 0.1;
	double timeStep = 0.1;	
	double bubbleBreakFreq = 0.1;
	/////////////////////////////below should be only set once
	mGrid.setDensity(density);
	mGrid.setViscosity(viscosity);
	mGrid.setScatterFreq(scatterFreq);
	mGrid.setScatterCoef(scatterCoef);
	mGrid.setBubbleBreakFreq(bubbleBreakFreq);
	mGrid.setTimeStep(timeStep);//time step
	mGrid.setBubbleRadius(bubbleRadiusMin, bubbleRadiusMax, bubbleRadiusCount);
	////////////////////////////////////////////////


	int	frame = (int)time.as( MTime::kFilm );
	if (frame == 0) 
		frame = 1;

	stringstream ss;
	string tmp;


	//**delete existing bubbles**//
	for(int i = 1; i <= bubbleRadiusCount ; ++i)
	{
		string particleExist = "particleExists bubbleParticle";

		ss<<i;
		ss>>tmp;
		particleExist += tmp + ";";
		ss.clear();

		char* particleExistChar=(char *)particleExist.c_str();
		MString strParticleExist = particleExistChar;

		int exist;
		MGlobal::executeCommand(strParticleExist, exist);
		if(exist)
		{
			string deleteParticle = "select -r bubbleParticle";
			deleteParticle += tmp + ";";
			char* deleteParticleChar=(char *)deleteParticle.c_str();
			MString strDeleteParticle = deleteParticleChar;

			MGlobal::executeCommand(strDeleteParticle);
			MGlobal::executeCommand("doDelete");
		}

	}



	mGrid.doSimulation(frame);
	for(int radiusIndex = 1; radiusIndex <= bubbleRadiusCount ; ++radiusIndex)
	{
		int pointNum;
		float* ptPos = getParticlePositions(frame, radiusIndex-1, containerSizeX, containerSizeY, containerSizeZ, viscosity, density, scatterFreq, scatterCoef, &pointNum);

		string index = "";
		ss<<radiusIndex;
		ss>>index;
		ss.clear();
		char* indexChar=(char *)index.c_str();
		MString strIndex = indexChar;

		string particleAttr = "";
		if(pointNum > 0)
			particleAttr += "particle ";
		for(int i = 0; i < pointNum ; ++i)
		{
			particleAttr += "-p ";

			ss<<ptPos[i * 3];
			ss>>tmp;
			particleAttr += tmp + " ";
			ss.clear();
			
			ss<<ptPos[i * 3 + 1];
			ss>>tmp;
			particleAttr += tmp + " ";
			ss.clear();
				
			ss<<ptPos[i * 3 + 2];
			ss>>tmp;
			particleAttr += tmp + " ";
			ss.clear();
		}
		particleAttr += "-c 1 -n bubbleParticle";

		char* particleAttrChar=(char *)particleAttr.c_str();
		MString strparticleAttr = particleAttrChar;
		MGlobal::executeCommand(strparticleAttr + strIndex);//create particles

		//string particleRenderAttr = "setAttr bubbleParticle";
		//particleRenderAttr += tmp;
		//particleRenderAttr += "Shape.particleRenderType 4";
		//char* particleRenderAttrChar=(char *)particleRenderAttr.c_str();
		//MString strParticleRenderAttr = particleRenderAttrChar;
		MGlobal::executeCommand("setAttr bubbleParticle" + strIndex + "Shape.particleRenderType 4");//render particle to sphere

		//MGlobal::executeCommand("setAttr bubbleParticle" + 1 + "Shape.particleRenderType 4");//render particle to sphere
	
		string enableRadiusSet = (string)"addAttr -is true -ln radius -at " +'\"' + (string)"float" +'\"' + (string)" -min 0 -max 10 -dv 0.5 bubbleParticle" + index + "Shape";
		char* enableRadiusSetChar=(char *)enableRadiusSet.c_str();
		MString strEnableRadiusSet = enableRadiusSetChar;
		MGlobal::executeCommand(strEnableRadiusSet);
	}

	//**set bubble size**//
	for(int i = 1; i <= bubbleRadiusCount ; ++i)
	{
		ss<<i;
		ss>>tmp;
		ss.clear();

		char* particleGroupCountChar=(char *)tmp.c_str();
		MString strParticleGroupCountChar = particleGroupCountChar;

		string bSize;
		ss<< mGrid.getBubbleRadius(i-1);
		ss>> bSize;
		ss.clear();
		char* bSizeChar=(char *)bSize.c_str();
		MString strBSize = bSizeChar;
		MGlobal::executeCommand("setAttr bubbleParticle"+ strParticleGroupCountChar +"Shape.radius " + strBSize);
	}

	//**dummy cube mesh**//
	MFnMesh	meshFS;
	double cubeSize = 0.001f * bubbleSize;
	MFloatPointArray points;
	const int numFaces = 6;
	int numVertices	= 8;
	const int numFaceConnects = 24;

	MFloatPoint vtx_1( -cubeSize, -cubeSize, -cubeSize );
	MFloatPoint vtx_2(  cubeSize, -cubeSize, -cubeSize );
	MFloatPoint vtx_3(  cubeSize, -cubeSize,  cubeSize );
	MFloatPoint vtx_4( -cubeSize, -cubeSize,  cubeSize );
	MFloatPoint vtx_5( -cubeSize,  cubeSize, -cubeSize );
	MFloatPoint vtx_6( -cubeSize,  cubeSize,  cubeSize );
	MFloatPoint vtx_7(  cubeSize,  cubeSize,  cubeSize );
	MFloatPoint vtx_8(  cubeSize,  cubeSize, -cubeSize );
	points.append( vtx_1 );
	points.append( vtx_2 );
	points.append( vtx_3 );
	points.append( vtx_4 );
	points.append( vtx_5 );
	points.append( vtx_6 );
	points.append( vtx_7 );
	points.append( vtx_8 );
	int face_counts[numFaces] = { 4, 4, 4, 4, 4, 4 };
	MIntArray faceCounts( face_counts, numFaces );
	int face_connects[ numFaceConnects ] = {0,1,2,3,4,5,6,7,3,2,6,5,0,3,5,4,0,4,7,1,1,7,6,2};
	MIntArray faceConnects( face_connects, numFaceConnects );
	MObject newMesh = meshFS.create(numVertices, numFaces, points, faceCounts, faceConnects, outData);



	return MS::kSuccess;
}



float* CreateBubbleNode::getParticlePositions(int frame, int radiusIndex, double containerSizeX, double containerSizeY, double containerSizeZ, double viscosity, double density, double scatterFreq, double scatterCoef, int* ptNum)
{

	int size;
	float* postest =  mGrid.getBubblePosition(radiusIndex, &size);// index
	float* pos = new float[size * 3];
	for(int i = 0 ; i < size ; ++i)
	{
		pos[i * 3]     = postest[i * 3] - containerSizeX / 2.0f;
		pos[i * 3 + 1] = postest[i * 3 + 1];
		pos[i * 3 + 2] = postest[i * 3 + 2] - containerSizeZ / 2.0f;
	}
	*ptNum = size;

	return pos;
}


