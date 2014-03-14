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





void* CreateBubbleNode::creator()
{
	return new CreateBubbleNode;
}


MStatus CreateBubbleNode::initialize()
{

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

	//creating the attribut
	//CreateBubbleNode::containerSizeX = defaultContainerSizeX.create( "sizex", "sx", MFnNumericData::kInt, 10, &returnStatus );
	//defaultContainerSizeX.setMax(30.0);
	//defaultContainerSizeX.setMin(5);

	//CreateBubbleNode::containerSizeY = defaultContainerSizeY.create( "sizey", "sy", MFnNumericData::kInt, 10, &returnStatus );
	//defaultContainerSizeY.setMax(30.0);
	//defaultContainerSizeY.setMin(5);

	//CreateBubbleNode::containerSizeZ = defaultContainerSizeZ.create( "sizez", "sz", MFnNumericData::kInt, 10, &returnStatus );
	//defaultContainerSizeZ.setMax(30.0);
	//defaultContainerSizeZ.setMin(5);

	//CreateBubbleNode::containerResolutionX = defaultContainerResolutionX.create( "reslx", "rx", MFnNumericData::kInt, 10, &returnStatus );
	//defaultContainerResolutionX.setMax(30.0);
	//defaultContainerResolutionX.setMin(5);

	//CreateBubbleNode::containerResolutionY = defaultContainerResolutionY.create( "resly", "ry", MFnNumericData::kInt, 10, &returnStatus );
	//defaultContainerResolutionY.setMax(30.0);
	//defaultContainerResolutionY.setMin(5);

	//CreateBubbleNode::containerResolutionZ = defaultContainerResolutionZ.create( "reslz", "rz", MFnNumericData::kInt, 10, &returnStatus );
	//defaultContainerResolutionZ.setMax(30.0);
	//defaultContainerResolutionZ.setMin(5);

	CreateBubbleNode::viscosity = defaultViscosity.create( "viscosity", "v", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultViscosity.setMax(10.0);
	defaultViscosity.setMin(0.1);

	CreateBubbleNode::density = defaultDensity.create( "density", "d", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultDensity.setMax(10.0);
	defaultDensity.setMin(0.1);

	CreateBubbleNode::scatterFreq = defaultScatterFreq.create( "scaterfreq", "f", MFnNumericData::kDouble, 1.0, &returnStatus );
	defaultScatterFreq.setMax(10.0);
	defaultScatterFreq.setMin(0.1);

	CreateBubbleNode::scatterCoef = defaultScatterCoef.create( "scatercoef", "c", MFnNumericData::kDouble, 1.0, &returnStatus );
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
		double bubbleSizes = bubbleSizeData.asDouble();

		double containSizeX,containSizeY,containSizeZ;
		double resolutionX, resolutionY, resolutionZ;
		MDoubleArray resolutionArray;
		MGlobal::executeCommand("getAttr fluidShape1.dimensionsW", containSizeX);
		MGlobal::executeCommand("getAttr fluidShape1.dimensionsH", containSizeY);
		MGlobal::executeCommand("getAttr fluidShape1.dimensionsD", containSizeZ);
		MGlobal::executeCommand("getAttr fluidShape1.resolution", resolutionArray);
		resolutionX = resolutionArray[0];
		resolutionY = resolutionArray[1];
		resolutionZ = resolutionArray[2];



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
	int	frame = (int)time.as( MTime::kFilm );
	if (frame == 0) 
		frame = 1;

	stringstream ss;
	string tmp;
	string particleAttr = "";



	int exist;
	MGlobal::executeCommand("particleExists bubbleParticle1;",exist);
	if(exist)
	{
		MGlobal::executeCommand("select -r bubbleParticle1;");
		MGlobal::executeCommand("doDelete");
	}

	int pointNum;
	float* ptPos = getParticlePositions(frame, viscosity, density, scatterFreq, scatterCoef, &pointNum);

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
	particleAttr += "-c 1 -n bubbleParticle1";

	char* particleAttrChar=(char *)particleAttr.c_str();
	MString strparticleAttr = particleAttrChar;
	MGlobal::executeCommand(strparticleAttr);//create particles
	MGlobal::executeCommand("setAttr bubbleParticle1Shape.particleRenderType 4");//render particle to sphere
	
	string enableRadiusSet = (string)"addAttr -is true -ln radius -at " +'\"' + (string)"float" +'\"' + (string)" -min 0 -max 10 -dv 0.5 bubbleParticle1Shape";
	char* enableRadiusSetChar=(char *)enableRadiusSet.c_str();
	MString strEnableRadiusSet = enableRadiusSetChar;
	MGlobal::executeCommand(strEnableRadiusSet);

	//**set bubble size**//
	string bSize;
	ss<<bubbleSize;
	ss>>bSize;
	ss.clear();
	char* bSizeChar=(char *)bSize.c_str();
	MString strBSize = bSizeChar;
	MGlobal::executeCommand("setAttr bubbleParticle1Shape.radius " + strBSize);


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



float* CreateBubbleNode::getParticlePositions(int frame, double viscosity, double density, double scatterFreq, double scatterCoef, int* ptNum)
{
	float* pos = new float[frame * 3];


	for(int i = 0 ; i < frame ; ++i)
	{
		pos[i * 3] = viscosity * (frame / 1000 + 1) * cos(2 * 3.1415926 / (float)frame * i);
		pos[i * 3 + 1] = 0;
		pos[i * 3 + 2] = viscosity * (frame / 1000 + 1) * sin(2 * 3.1415926 / (float)frame * i);
	}

	*ptNum = frame;
	return pos;
}


//MStatus LSystemNode::computeMesh( const MPlug& plug, MDataBlock& block )
////
////	Description:
////
////		To illustrate a custom use of the arrayMapper node, we'll define this
////	alternative compute that will map particles to mesh vertex positions. If
////	array lengths do not match, we'll wrap around the arrays.
////
////	How to use:
////
////		1) Create a poly surface.
////		2) connectAttr polyShape1.outMesh particleAttr1.computeNode
////		3) connectAttr particleShape1.count particleAttr1.particleCount
////		4) connectAttr particleAttr1.outPosition particleShape1.rampPosition (*)
////
////		* You may choose to drive any vector based per particle attribute.
////		  rampPosition was selected here as an example.
////
//{
//	MStatus status = MS::kSuccess;
//
//	// Verify that computeNode has a mesh connected to computeNode
//	//
//	MPlug compPlug( thisMObject(), computeNode );
//	MPlugArray conns;
//	compPlug.connectedTo( conns, true, false, &status );
//	if( conns.length() <= 0 )
//	{
//		return MS::kFailure;
//	}
//	MPlug conn = conns[0];
//	MObject compNode = conn.node();
//	MFnMesh meshFn( compNode, &status );
//	if( status != MS::kSuccess )
//	{
//		return MS::kFailure;
//	}
//
//	// Retrieve the mesh vertices
//	//
//	MPointArray points;
//	meshFn.getPoints( points );
//	unsigned int nPoints = points.length();
//
//	// Retrieve the current particle count.
//	//
//	// NOTE: Due to the order in which the particle system requests
//	// various pieces of data, some attributes are requested prior
//	// to the actual emission of particles (eg. rampPosition), whereas
//	// other attributes are requested after particles have been emitted.
//	//
//	// If the driven PP attribute on the particleShape is requested prior
//	// to the emission of particles, this compute() method will be called
//	// before any particles have been emitted. As a result, the effect
//	// will be lagged by one frame.
//	//
//	unsigned int nParticles = 0;
//	int nSignedPart = block.inputValue( particleCount ).asInt();
//	if( nSignedPart > 0 )
//	{
//		nParticles = nSignedPart;
//	}
//
//	// Get pointer to destination attribute: outPositionPP
//	//
//	MFnVectorArrayData dataVectorArrayFn;
//	MVectorArray outPosArray;
//	MObject posD = block.outputValue( outPositionPP ).data();
//	//const char* typeStr = posD.apiTypeStr();
//	status = dataVectorArrayFn.setObject( posD );
//	if( status == MS::kSuccess )
//	{
//		outPosArray = dataVectorArrayFn.array();
//	}
//
//	outPosArray.setLength( nParticles );
//	for( unsigned int i = 0; i < nParticles; i++ )
//	{
//		unsigned int index = i;
//		if( nParticles > nPoints )
//		{
//			index = i % nPoints;
//		}
//		MPoint point = points[index];
//		MVector pos( point.x, point.y, point.z );
//		outPosArray[i] = pos;
//	}
//	dataVectorArrayFn.set( outPosArray );
//	return MS::kSuccess;
//}
