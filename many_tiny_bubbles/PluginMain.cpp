#include <maya/MPxCommand.h>
#include <maya/MFnPlugin.h>
#include <maya/MIOStream.h>
#include <maya/MString.h>
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSimple.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MDGModifier.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MStringArray.h>
#include <list>

#include "CreateBubbleCmd.h"
#include "LSystemNode.h"

MString mllPath;

MStatus initializePlugin( MObject obj )
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj, "MyPlugin_CreateBubble", "1.11", "Any");

    // Register Command
    status = plugin.registerCommand( "CreateBubbleCmd", CreateBubbleCmd::creator );
    if (!status) {
        status.perror("registerCommand");
        return status;
    }

	status = plugin.registerNode("LSystemNode", LSystemNode::id,LSystemNode::creator, LSystemNode::initialize);
			if (!status) {
		status.perror("registerNode");
		return status;
	}
	mllPath = plugin.loadPath();

	MString command ="source \"" + plugin.loadPath() + "/gui.mel\";";
	MGlobal::executeCommand(command);

    return status;
}

MStatus uninitializePlugin( MObject obj)
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj );
	MGlobal::executeCommand( " if (`menu -exists myMenu`) {deleteUI myMenu;}");

    status = plugin.deregisterCommand( "CreateBubbleCmd" );
    if (!status) {
	    status.perror("deregisterCommand");
	    return status;
    }
	
	status = plugin.deregisterNode(LSystemNode::id);
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

    return status;
}


