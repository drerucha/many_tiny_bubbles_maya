// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b


// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern int theDim[3];
extern double theCellSize;
extern double bubbleRadius;

// TODO: add constants for dispersed bubble simulation

class Constants
{
	public:
		Constants( void );
		~Constants( void );

		static void setContainerDim(int x, int y, int z);
		static void setGridSize(double sx);
		static void setBubbleRadius(double radius);
	private:
};

#endif