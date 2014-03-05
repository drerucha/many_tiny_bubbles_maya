#include "LSystemNode.h"
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

extern MString mllPath;

MStatus returnStatus;

#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
	}


MObject LSystemNode::time;
MObject LSystemNode::outputMesh;
MObject LSystemNode::angle;
MObject LSystemNode::randAngle;
MObject LSystemNode::stepSize;
MObject LSystemNode::grammar;
MTypeId LSystemNode::id( 0x00128 );

void* LSystemNode::creator()
{
	return new LSystemNode;
}


MStatus LSystemNode::initialize()
{
	MFnNumericAttribute defaultAngle;
	MFnNumericAttribute defaultStepSize;
	MFnNumericAttribute defaultRandAngle;
	MFnTypedAttribute defaultGrammar;
	MFnUnitAttribute defaultTime;
	MFnTypedAttribute defaultGeometry;

	//creating the attribut
	LSystemNode::time = defaultTime.create( "time", "tm", MFnUnitAttribute::kTime, 2.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating LSystemNode time attribute\n");

	LSystemNode::outputMesh = defaultGeometry.create( "outputMesh", "out", MFnData::kMesh, &returnStatus ); 
	McheckErr(returnStatus, "ERROR creating LSystemNode output attribute\n");
	defaultGeometry.setStorable(false);
	defaultGeometry.setStorable(false);

	LSystemNode::angle = defaultAngle.create( "angle", "ag", MFnNumericData::kDouble, 30, &returnStatus );
	McheckErr(returnStatus, "ERROR creating LSystemNode angle attribute\n");
	defaultAngle.setMax(180.0);
	defaultAngle.setMin(0.0);

	LSystemNode::randAngle = defaultRandAngle.create( "randAngle", "ra", MFnNumericData::kDouble, 10, &returnStatus );
	McheckErr(returnStatus, "ERROR creating LSystemNode angle attribute\n");
	defaultRandAngle.setMax(30.0);
	defaultRandAngle.setMin(0.0);

	LSystemNode::stepSize = defaultStepSize.create( "stepSize", "ss", MFnNumericData::kDouble, 1.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating LSystemNode stepSize attribute\n");
	defaultStepSize.setMax(10.0);
	defaultStepSize.setMin(0.1);

	MFnStringData fnStringData;
	MObject defaultString;
	defaultString = fnStringData.create( "simple1.txt" );
	LSystemNode::grammar = defaultGrammar.create( "grammarFile", "gm", MFnData::kString, defaultString, &returnStatus ); 
	McheckErr(returnStatus, "ERROR creating LSystemNode grammar attribute\n");



	//adding the attribute
	returnStatus = addAttribute(LSystemNode::time);
	McheckErr(returnStatus, "ERROR adding time attribute\n");

	returnStatus = addAttribute(LSystemNode::angle);
	McheckErr(returnStatus, "ERROR adding angle attribute\n");

	returnStatus = addAttribute(LSystemNode::randAngle);
	McheckErr(returnStatus, "ERROR adding angle attribute\n");

	returnStatus = addAttribute(LSystemNode::stepSize);
	McheckErr(returnStatus, "ERROR adding stepSize attribute\n");

	returnStatus = addAttribute(LSystemNode::grammar);
	McheckErr(returnStatus, "ERROR adding grammar attribute\n");

	returnStatus = addAttribute(LSystemNode::outputMesh);
	McheckErr(returnStatus, "ERROR adding outputMesh attribute\n");




	//assigning Dependencies
	returnStatus = attributeAffects(LSystemNode::time,
								    LSystemNode::outputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(LSystemNode::angle,
								    LSystemNode::outputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(LSystemNode::randAngle,
								    LSystemNode::outputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(LSystemNode::stepSize,
								    LSystemNode::outputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(LSystemNode::grammar,
								    LSystemNode::outputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	return MS::kSuccess;
}


MStatus LSystemNode::compute(const MPlug& plug, MDataBlock& data)
{
	MStatus returnStatus;

	if (plug == outputMesh) {
		/* Get time */
		MDataHandle timeData = data.inputValue(time, &returnStatus ); 
		McheckErr(returnStatus, "Error getting time data handle\n");
		MTime times = timeData.asTime();

		MDataHandle grammarData = data.inputValue(grammar, &returnStatus );	
		McheckErr(returnStatus, "Error getting time data handle\n");
		MString grammars = grammarData.asString();

		MDataHandle angleData = data.inputValue(angle, &returnStatus );
		McheckErr(returnStatus, "Error getting time data handle\n");
		double angles = angleData.asDouble();

		MDataHandle randAngleData = data.inputValue(randAngle, &returnStatus );
		McheckErr(returnStatus, "Error getting time data handle\n");
		double randAngles = randAngleData.asDouble();

		MDataHandle stepSizeData = data.inputValue(stepSize, &returnStatus );
		McheckErr(returnStatus, "Error getting time data handle\n");
		double steps = stepSizeData.asDouble();



		/* Get output object */

		MDataHandle outputHandle = data.outputValue(outputMesh, &returnStatus);
		McheckErr(returnStatus, "ERROR getting polygon data handle\n");

		MFnMeshData dataCreator;
		MObject newOutputData = dataCreator.create(&returnStatus);
		McheckErr(returnStatus, "ERROR creating outputData");

		createMesh(times, grammars, angles, steps, newOutputData, randAngles);
		McheckErr(returnStatus, "ERROR creating new Plant");

		outputHandle.set(newOutputData);
		data.setClean( plug );
	} else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}

MObject LSystemNode::createMesh(const MTime &time, MString &grammarFiles, double &angles, double &steps, MObject &outData, double &randAngles)
{
	int	frame;

	string grammar="F\nF->F[+F]F[-F]F";
	string folderPath;
	string path = mllPath.asChar();

	int slashCount = 0;
	for(int j = 1; j <=path.length(); j++)
	{
		if(path.at(path.length()-j) == '/')
		{
			//slashCount++;

			//if(slashCount == 2)
				folderPath = path.substr(0, path.length() - j + 1);
				break;
		}	
	}
	folderPath = folderPath + "plants/" + grammarFiles.asChar();

	ifstream inFile(folderPath, ios::in);
	istreambuf_iterator<char> begs(inFile), ends;
	string grammarstr(begs, ends);
	inFile.close();

	if(grammarstr!="")
		grammar = grammarstr;




	int iter = 2;

	frame = (int)time.as( MTime::kFilm );
	if (frame == 0) 
		frame = 1;
	iter = frame;

	//double randNum;
	//unsigned seed = (unsigned)time(NULL);
	//srand(seed);
	//for(int i = 0 ; i < 10 ;i++)
	//{
	//	randNum = rand() % 10;
	//}

	//create the branches
	LSystem plant;
	plant.setDefaultStep(steps);
	plant.setDefaultAngle(angles);
	plant.setRandAngleRange(randAngles);
	plant.loadProgramFromString(grammar);

	vector<LSystem::Branch> branches;
	plant.process(iter, branches);

	double cylinderRadius = 0.1;
	int numVertices;
	int numFaces;
	MPointArray points;
	MIntArray faceConnects;
	MIntArray faceCounts;

	for(int i = 0; i < branches.size(); i++)
	{	
		vec3 startPos;
		vec3 endPos;

		startPos = branches[i].first;
		endPos = branches[i].second;
		double sPosX = startPos[0];
		double sPosY = startPos[2];
		double sPosZ = startPos[1];
		double ePosX = endPos[0];
		double ePosY = endPos[2];
		double ePosZ = endPos[1];

		MPoint startMPoint = MPoint(sPosX, sPosY, sPosZ, 1.0);
		MPoint endMPoint = MPoint(ePosX, ePosY, ePosZ, 1.0);
		


		CylinderMesh newCylinder(startMPoint, endMPoint, cylinderRadius);		
		newCylinder.appendToMesh(points, faceCounts, faceConnects);
	}


	MFnMesh	meshFS;
	numVertices = points.length();
	numFaces = faceCounts.length();
	MObject newMesh = meshFS.create(numVertices, numFaces,points, faceCounts, faceConnects, outData);

	return newMesh;
}
