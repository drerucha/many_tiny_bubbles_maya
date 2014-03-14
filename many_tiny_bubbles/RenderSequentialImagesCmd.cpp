#include "RenderSequentialImagesCmd.h"
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MGlobal.h>
#include <list>
#include <sstream>

using namespace std;
RenderSequentialImagesCmd::RenderSequentialImagesCmd() : MPxCommand()
{
}

RenderSequentialImagesCmd::~RenderSequentialImagesCmd() 
{
}

const char *startFrameFlag = "-s", *startFrameLongFlag = "-startframe";//short flag can't be longer than 3 char, and long flag can't be less than 4 char
const char *endFrameFlag = "-e", *endFrameLongFlag = "-endframe";



//
MSyntax RenderSequentialImagesCmd::syntax()
{
    MSyntax syntax;

	syntax.addFlag(startFrameFlag, startFrameLongFlag, MSyntax::kLong);
	syntax.addFlag(endFrameFlag, endFrameLongFlag, MSyntax::kLong);


	return syntax;
}

MStatus RenderSequentialImagesCmd::doIt( const MArgList& args )
{
	MString name;
	MString id;

	MGlobal::displayInfo("Render Sequential Images");
	
	int start = 1;
	int end = 2;


	MArgDatabase argData(syntax(), args );//need to include <maya/MArgDatabase.h>

	if(argData.isFlagSet(startFrameFlag))
		argData.getFlagArgument(startFrameFlag, 0, start);

	if(argData.isFlagSet(endFrameFlag))
		argData.getFlagArgument(endFrameFlag, 0, end);


	if(start < 1)
		start = 1;
	if(end < start)
		end = start;

	//**should add warning message**//
	//TODO





	//////////////////////for rendering////////////////////////

	string renderType = (string)"setAttr " + '\"' + (string)"defaultRenderGlobals.currentRenderer" + '\"' + (string)" -type " + '\"' + (string)"string" + '\"' +" "+ '\"'+ (string)"mentalRay"+ '\"';
	char* renderTypeChar=(char *)renderType.c_str();
	MString strRenderType = renderTypeChar;
	MGlobal::executeCommand(strRenderType);//set render type to mental ray
	MGlobal::executeCommand("setAttr defaultRenderGlobals.outFormatControl 0");
	MGlobal::executeCommand("setAttr defaultRenderGlobals.animation 1");
	MGlobal::executeCommand("setAttr defaultRenderGlobals.putFrameBeforeExt 1");
	MGlobal::executeCommand("setAttr defaultRenderGlobals.periodInExt 1");
	MGlobal::executeCommand("setAttr defaultRenderGlobals.extensionPadding 4");//Frame padding
	MGlobal::executeCommand("setAttr defaultRenderGlobals.imageFormat 8");//jpg 8; bmp 20
	MGlobal::executeCommand("setAttr defaultRenderGlobals.startFrame 1");//set the start frame
	MGlobal::executeCommand("setAttr defaultRenderGlobals.endFrame 20");//set the end frame
	MGlobal::executeCommand("setAttr defaultResolution.width 1280");
	MGlobal::executeCommand("setAttr defaultResolution.height 720");

	//renderWindowSaveImageCallback "renderView" "D:/Course_Advanced Animation/Authuring tool/Final Project/teste001" "JPEG";
	for(int i = 1 ; i <= end - start + 1 ; ++i)
	{
		stringstream ss;
		string frameCount;
		string frameIndex;

		ss<<i;
		ss>>frameCount;
		ss.clear();

		ss<<start + i - 1;
		ss>>frameIndex;
		ss.clear();

		string saveRender = (string)"renderWindowSaveImageCallback " + '\"' + (string)"renderView" + '\"' + " " + '\"' + (string)"D:/Course_Advanced Animation/Authuring tool/Final Project/test" + frameCount + '\"' + " "+ '\"'+ (string)"JPEG"+ '\"';
		char* saveRenderChar=(char *)saveRender.c_str();
		MString strSaveRender = saveRenderChar;

		char* frameIndexChar=(char *)frameIndex.c_str();
		MString strFrameIndex = frameIndexChar;

		MGlobal::executeCommand("currentTime " + strFrameIndex);
		MGlobal::executeCommand("renderIntoNewWindow render");
		MGlobal::executeCommand(strSaveRender);
	}
	return MStatus::kSuccess;
}

