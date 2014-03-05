#include <maya/MTime.h>
#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>

#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MFnMeshData.h>

#include <maya/MIOStream.h>
#include "LSystem.h"
#include "cylinder.h"
#include <string>
#include <Windows.h>
class LSystemNode : public MPxNode
{
public:
					LSystemNode() {};
	virtual 		~LSystemNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	static MObject	time;
	static MObject	outputMesh;
	static MObject	angle;
	static MObject	stepSize;
	static MObject	grammar;
	static MTypeId	id;
	static MObject  randAngle;
	//string WstringToString(const wstring);
protected:
	//MObject createMesh(const MTime& time, MObject& outData, MStatus& stat);
	MObject createMesh(const MTime& time, MString &grammarFiles, double &angles, double &steps, MObject& outData, double &randAngles);
};