#include "CreateBubbleCmd.h"
//#include "LSystem.h"
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MGlobal.h>
#include <list>
#include <sstream>

using namespace std;
CreateBubbleCmd::CreateBubbleCmd() : MPxCommand()
{
}

CreateBubbleCmd::~CreateBubbleCmd() 
{
}

const char *containerSizeXFlag = "-sx", *containerSizeXLongFlag = "-sizex";//short flag can't be longer than 3 char, and long flag can't be less than 4 char
const char *containerSizeYFlag = "-sy", *containerSizeYLongFlag = "-sizey";
const char *containerSizeZFlag = "-sz", *containerSizeZLongFlag = "-sizez";
const char *containerResXFlag  = "-rx", *containerResXLongFlag = "-reslx";
const char *containerResYFlag  = "-ry", *containerResYLongFlag = "-resly";
const char *containerResZFlag  = "-rz", *containerResZLongFlag = "-reslz";
const char *bubbleNumFlag      = "-n", *bubbleNumLongFlag     = "-bubblenum";
const char *fluidViscosityFlag = "-v", *fluidViscosityLongFlag = "-viscosity";
const char *fluidDensityFlag   = "-d", *fluidDensityLongFlag = "-density";
const char *atmosphereFlag     = "-a", *atmosphereLongFlag = "-atmosphere";
const char *scaterFreqFlag     = "-f", *scaterFreqLongFlag = "-scaterfreq";
const char *scaterCoefFlag     = "-c", *scaterCoefLongFlag = "-scatercoef";
const char *bblBreakFreqFlag   = "-bf", *bblBreakFreqLongFlag = "-breakfreq";


//
MSyntax CreateBubbleCmd::syntax()
{
    MSyntax syntax;

	syntax.addFlag(containerSizeXFlag, containerSizeXLongFlag, MSyntax::kLong);
	syntax.addFlag(containerSizeYFlag, containerSizeYLongFlag, MSyntax::kLong);
	syntax.addFlag(containerSizeZFlag, containerSizeZLongFlag, MSyntax::kLong);
	syntax.addFlag(containerResXFlag,  containerResXLongFlag,  MSyntax::kLong);
	syntax.addFlag(containerResYFlag,  containerResYLongFlag,  MSyntax::kLong);
	syntax.addFlag(containerResZFlag,  containerResZLongFlag,  MSyntax::kLong);
	syntax.addFlag(bubbleNumFlag,      bubbleNumLongFlag,      MSyntax::kLong);

	syntax.addFlag(fluidViscosityFlag, fluidViscosityLongFlag, MSyntax::kDouble);
    syntax.addFlag(fluidDensityFlag,   fluidDensityLongFlag,   MSyntax::kDouble);
	syntax.addFlag(atmosphereFlag,     atmosphereLongFlag,     MSyntax::kDouble);
    syntax.addFlag(scaterFreqFlag,     scaterFreqLongFlag,     MSyntax::kDouble);
	syntax.addFlag(scaterCoefFlag,     scaterCoefLongFlag,     MSyntax::kDouble);
    syntax.addFlag(bblBreakFreqFlag,   bblBreakFreqLongFlag,   MSyntax::kDouble);

	return syntax;
}

MStatus CreateBubbleCmd::doIt( const MArgList& args )
{
	
	MString name;
	MString id;

	MGlobal::displayInfo("Implement me");
	
	MString grammar=" ";
	int bubbleNum = 1;
	int containerSizeX = 5;
	int containerSizeY = 5;
	int containerSizeZ = 5;
	int containerResolutionX = 5;
	int containerResolutionY = 5;
	int containerResolutionZ = 5;
	double containerTransX = 0.0;
	double containerTransY = 0.0;
	double containerTransZ = 0.0;
	double fluidViscosity = 1.0;
	double fluidDensity = 1.0;
	double atmosphere = 1.0;
	double scaterFreq = 1.0;
	double scaterCoef = 1.0;
	double bblBreakFreq = 1.0;

	MArgDatabase argData(syntax(), args );//need to include <maya/MArgDatabase.h>
#pragma region getAttributes

	if(argData.isFlagSet(containerSizeXFlag))
		argData.getFlagArgument(containerSizeXFlag, 0, containerSizeX);

	if(argData.isFlagSet(containerSizeYFlag))
		argData.getFlagArgument(containerSizeYFlag, 0, containerSizeY);

	if(argData.isFlagSet(containerSizeZFlag))
		argData.getFlagArgument(containerSizeZFlag, 0, containerSizeZ);

	if(argData.isFlagSet(containerResXFlag))
		argData.getFlagArgument(containerResXFlag, 0, containerResolutionX);

	if(argData.isFlagSet(containerResYFlag))
		argData.getFlagArgument(containerResYFlag, 0, containerResolutionY);

	if(argData.isFlagSet(containerResZFlag))
		argData.getFlagArgument(containerResZFlag, 0, containerResolutionZ);

	if(argData.isFlagSet(bubbleNumFlag))
		argData.getFlagArgument(bubbleNumFlag, 0, bubbleNum);

	if(argData.isFlagSet(fluidViscosityFlag))
		argData.getFlagArgument(fluidViscosityFlag, 0, fluidViscosity);

	if(argData.isFlagSet(fluidDensityFlag))
		argData.getFlagArgument(fluidDensityFlag, 0, fluidDensity);

	if(argData.isFlagSet(atmosphereFlag))
		argData.getFlagArgument(atmosphereFlag, 0, atmosphere);

	if(argData.isFlagSet(scaterFreqFlag))
		argData.getFlagArgument(scaterFreqFlag, 0, scaterFreq);

	if(argData.isFlagSet(scaterCoefFlag))
		argData.getFlagArgument(scaterCoefFlag, 0, scaterCoef);

	if(argData.isFlagSet(bblBreakFreqFlag))
		argData.getFlagArgument(bblBreakFreqFlag, 0, bblBreakFreq);
#pragma endregion

	containerTransX = 0.0f;
	containerTransY = (double)containerSizeY / 2.0f;
	containerTransZ = 0.0f;



#pragma region get container attributes
	stringstream ss;
	string tmp;
	string container;
	string containerTran;
	
	ss<<containerResolutionX;
	ss>>tmp;
	container += tmp + " ";
	ss.clear();

	ss<<containerResolutionY;
	ss>>tmp;
	container += tmp + " ";
	ss.clear();

	ss<<containerResolutionZ;
	ss>>tmp;
	container += tmp + " ";
	ss.clear();

	ss<<containerSizeX;
	ss>>tmp;
	container += tmp + " ";
	ss.clear();

	ss<<containerSizeY;
	ss>>tmp;
	container += tmp + " ";
	ss.clear();

	ss<<containerSizeZ;
	ss>>tmp;
	container += tmp + " ";
	ss.clear();

	ss<<containerTransX;
	ss>>tmp;
	containerTran += tmp + " ";
	ss.clear();

	ss<<containerTransY;
	ss>>tmp;
	containerTran += tmp + " ";
	ss.clear();

	ss<<containerTransZ;
	ss>>tmp;
	containerTran += tmp + " ";
	ss.clear();

	char* containerChar=(char *)container.c_str();
	MString strContainer = containerChar;
	char* containerTranChar=(char *)containerTran.c_str();
	MString strContainerTran = containerTranChar;
#pragma endregion

	MGlobal::executeCommand("create3DFluid " + strContainer);//the first three parameters are for resolution, the last three are for size
	MGlobal::executeCommand("move -r " + strContainerTran);

	//MGlobal::executeCommand("setAttr fluidShape1.surfaceRender 1");//暫時先不用
	//MGlobal::executeCommand("setAttr fluidShape1.softSurface 1");
	//MGlobal::executeCommand("setAttr fluidShape1.densityBuoyancy -1.5");
	//MGlobal::executeCommand("setAttr fluidShape1.surfaceThreshold 0.133734");


	//MGlobal::executeCommand("fluidEmitter -pos 0 0 0 -type omni -der 1 -her 1 -fer 1 -fdr 2 -r 100.0 -cye none -cyi 1 -mxd 1 -mnd 0");
	//MGlobal::executeCommand("emitter -pos 0 0 0 -type omni -r 100 -sro 0 -nuv 0 -cye none -cyi 1 -spd 1 -srn 0 -nsp 1 -tsp 0 -mxd 0 -mnd 0 -dx 1 -dy 0 -dz 0 -sp 0");
	//MGlobal::executeCommand("particle");
	//MGlobal::executeCommand("connectDynamic -em emitter1 particle1");
	//MGlobal::executeCommand("setAttr emitter1.rate 1");//bubble count
	
	//MGlobal::executeCommand("playbackOptions -ast 1 -aet 100");//Start and end of the animation:
	MGlobal::executeCommand("playbackOptions -min 1 -max 20");//Start and end of playback time range:

	//MGlobal::executeCommand("currentTime 500");
	//int a = MGlobal::executeCommand("particle -count -q particle1");
	//int b = 0;
	//MGlobal::executeCommand("changeSelectMode -component");
	//MGlobal::executeCommand("hilite particle1") ;
	//MGlobal::executeCommand("select -r particle1.pt[1]") ;
	//MGlobal::executeCommand("particle -e -or 1 -at position -vv 1 -0.633177 -1.849741 particleShape1");
	//MGlobal::executeCommand("particle -p 1.40844 0 -1.014435 -p 1.301524 0 1.368541 -c 1");


	//MFnParticleSystem



    return MStatus::kSuccess;
}

