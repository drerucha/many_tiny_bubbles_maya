#ifndef RenderSequentialImagesCmd_H_
#define RenderSequentialImagesCmd_H_

#include <maya/MPxCommand.h>
#include <string>
#include <maya/MDGModifier.h>
#include <maya/MPlug.h>


class RenderSequentialImagesCmd : public MPxCommand
{
public:
    RenderSequentialImagesCmd();
    virtual ~RenderSequentialImagesCmd();
    static void* creator() { return new RenderSequentialImagesCmd(); }
    MStatus doIt( const MArgList& args );
	static MSyntax syntax();
};

#endif