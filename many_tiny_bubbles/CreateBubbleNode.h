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
#include "mac_grid.h"
class CreateBubbleNode : public MPxNode
{
public:
					CreateBubbleNode() {};
	virtual 		~CreateBubbleNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	static MObject time;
	static MObject outputMesh;
	static MObject containerSizeX;
	static MObject containerSizeY;
	static MObject containerSizeZ;
	static MObject containerResolutionX;
	static MObject containerResolutionY;
	static MObject containerResolutionZ;
	static MObject viscosity;
	static MObject density;
	static MObject scatterFreq;
	static MObject scatterCoef;
	static MObject bubbleSize;
	static MTypeId	id;


	//string WstringToString(const wstring);
protected:
	MStatus	computeMesh(const MTime& time, const MPlug& plug, MDataBlock& block );
	MStatus	createBubble(const MTime& time, MObject& outData, double &containerSizeX, double &containerSizeY, double &containerSizeZ, double &viscosity, 
						 double &density, double &scatterFreq, double &scatterCoef, double &bubbleSize, const MPlug& plug, MDataBlock& block );
	float* getParticlePositions(int frame, double viscosity, double density, double scatterFreq, double scatterCoef, int* ptNum);
	MACGrid mGrid;
};