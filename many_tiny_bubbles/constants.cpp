// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
int theDim[3] = {50, 50, 50};
#else
int theDim[3] = {2, 2, 1};//12, 12, 4
#endif

double theCellSize = 0.5;
//double bubbleRadius = 0.03;

// TODO: add constants for dispersed bubble simulation

Constants::Constants()
{
}

Constants::~Constants()
{
}

void Constants::setContainerDim(int x, int y, int z)
{
	theDim[0] = x;
	theDim[1] = y;
	theDim[2] = z;
}

void Constants::setGridSize(double sx)
{
	theCellSize = sx;
}

void Constants::setBubbleRadius(double radius)
{
	//bubbleRadius = radius;
}