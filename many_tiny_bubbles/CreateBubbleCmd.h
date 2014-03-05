#ifndef CreateBubbleCmd_H_
#define CreateBubbleCmd_H_

#include <maya/MPxCommand.h>
#include <string>
#include <maya/MDGModifier.h>
#include <maya/MPlug.h>


class CreateBubbleCmd : public MPxCommand
{
public:
    CreateBubbleCmd();
    virtual ~CreateBubbleCmd();
    static void* creator() { return new CreateBubbleCmd(); }
    MStatus doIt( const MArgList& args );
	static MSyntax syntax();
};

#endif