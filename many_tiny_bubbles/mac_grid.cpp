// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
//#include "open_gl_headers.h"
//#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>
#include <ctime>


#define PI 3.1415926 
// Globals:
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define UNIT_X vec3( 1.0f, 0.0f, 0.0f )
#define UNIT_Y vec3( 0.0f, 1.0f, 0.0f )
#define UNIT_Z vec3( 0.0f, 0.0f, 1.0f )

double FLUID_DENSITY = 1.0f;
const double INITIAL_TEMPERATURE = 0.0f;
const double PARTICLE_MASS = 1.0f;
double FLUID_VISCOSITY = 1.5f;
double SCATTER_FREQUENCY = 1.0f;
double SCATTER_COEFFICIENT = 0.5f;
double BUBBLE_BREAK_FREQ = 0.5f;
double TIME_STEP = 0.1;
int currentFrame = 0;




MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
	containerSize[0] = 10;
	containerSize[1] = 10;
	containerSize[2] = 1;
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
	if (&orig == this)
	{
		return *this;
	}
	mU = orig.mU;
	mV = orig.mV;
	mW = orig.mW;
	mP = orig.mP;
	mD = orig.mD;
	mT = orig.mT;   
	mTemp = orig.mTemp;
	mConfForceX = orig.mConfForceX;
	mConfForceY = orig.mConfForceY;
	mConfForceZ = orig.mConfForceZ;
	mFractionField = orig.mFractionField;
	mScatterOdd = orig.mScatterOdd;

	//bubbleData = 

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
	mU.initialize();
	mV.initialize();
	mW.initialize();
	mP.initialize();
	mD.initialize();
	mT.initialize( INITIAL_TEMPERATURE );
   	mTemp.initialize();
	mConfForceX.initialize();
	mConfForceY.initialize();
	mConfForceZ.initialize();
	mFractionField.initialize(1.0f); // set default fraction field 
	mScatterOdd.initialize();
	srand(time(NULL));
	setUpAMatrix();
	bubblePosList.clear();
	bubbleVelList.clear();

}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: set initial values for density, temperature, and velocity

	mFractionField.initialize(1.0f); // set default fraction field 

	mT(theDim[MACGrid::X] * theCellSize / 2.0f, theCellSize, theDim[MACGrid::Z] * theCellSize / 2.0f) = 280.0f;
	mD(theDim[MACGrid::X] * theCellSize / 2.0f, theCellSize, theDim[MACGrid::Z] * theCellSize / 2.0f) = 1.0f;
	//mU(0, 0, 0) = 1.0;
	//mV( 4, 1, 0 ) = 1.0f;
	//mV( 4, 8, 0 ) = -1.0f;

}


////////////////New added function////////////////////

void MACGrid::generateBubbles()
{
	//bubblePos.push_back(vec3(5,1,5));
	//bubbleVel.push_back(vec3(0,1,0));



	for(int k = 0; k < bubbleRadiusList.size(); k++)
	{
		bubblePos bubblePosition;
		bubbleVel bubbleVelocity;

		for(int i = 0; i < 180 ; i += 100)
		{
			for(int j = 0; j < 360 ; j += 100)
			{
				double ran_num = (rand()%10)+1; //theoretically, the bigger bubble is less occured
				if(ran_num > 5)
				{
					double theta = i / 180.0f * PI;
					double phi = j / 180.0f * PI;

					vec3 center(theDim[MACGrid::X] * theCellSize / 2.0f, theCellSize, theDim[MACGrid::Z] * theCellSize / 2.0f);
					double radius = 0.1f;
					double x = theDim[MACGrid::X] * theCellSize / 2.0f + radius * sin(theta) * cos(phi);
					double y = theCellSize + radius * cos(theta);
					double z = theDim[MACGrid::Z] * theCellSize / 2.0f + radius * sin(theta) * sin(phi);

					
					if(bubblePosList.size() == bubbleRadiusList.size())
						bubblePosList[k].push_back(vec3(x, y, z));
					else
						bubblePosition.push_back(vec3(x, y, z));

					if(bubbleVelList.size() == bubbleRadiusList.size())
						bubbleVelList[k].push_back(vec3(0,1,0));
					else
						bubbleVelocity.push_back(vec3(0,1,0));

				}
			}
		}
		if(bubblePosList.size() != bubbleRadiusList.size())
			bubblePosList.push_back(bubblePosition);
		if(bubbleVelList.size() != bubbleRadiusList.size())
			bubbleVelList.push_back(bubbleVelocity);
	}


}

void MACGrid::computeScatterOdd()
{
	FOR_EACH_CELL
	{

	}
}

void MACGrid::computeFractionField()
{
	int j = 0;
	for(std::vector<bubblePos>::iterator iterRadius = bubblePosList.begin(); iterRadius != bubblePosList.end(); ++iterRadius, ++j)
	{	
		int i = 0;
		for(std::vector<vec3>::iterator iter = bubblePosList[j].begin(); /*i < bubblePosList[j].size()*/ iter != bubblePosList[j].end(); ++iter, ++i)
		{
			vec3 position = bubblePosList[j][i];
			int position_grid_X = (int)(position[0] / theCellSize);
			int position_grid_Y = (int)(position[1] / theCellSize);
			int position_grid_Z = (int)(position[2] / theCellSize);

			mFractionField(position_grid_X, position_grid_Y, position_grid_Z) -= 4.0f / 3.0f * PI * pow(bubbleRadiusList[j] / theCellSize, 3);
		}
		int a = 0;
	}
}

void MACGrid::advectBubbles( double dt )
{
	int j = 0;
	for(std::vector<bubblePos>::iterator iterRadius = bubblePosList.begin(); iterRadius != bubblePosList.end(); ++iterRadius, ++j)
	{
		int i = 0;
		for(std::vector<vec3>::iterator iterPos = bubblePosList[j].begin(), iterVel = bubbleVelList[j].begin(); iterPos != bubblePosList[j].end(); ++i)
		{


			vec3 position = bubblePosList[j][i];
			int position_grid_X = (int)(position[0] / theCellSize);
			int position_grid_Y = (int)(position[1] / theCellSize);
			int position_grid_Z = (int)(position[2] / theCellSize);

			//vec3 velocity = getVelocity(position);
			//vec3 velocity(0,5,0);

			//bubbleVel[i] += vec3(0,1,0);
			//vec3 velocity = bubbleVel[i];
			//velocity = vec3(0,1,0);
			vec3 velocity = vec3(0, bubbleVelList[j][i].Length(), 0);


			double ran_num = (rand()%100) / 100.0f;
			//SCATTER_FREQUENCY = 0.1;
			double velocityMag = velocity.Length(); 
			double fractionField = mFractionField(position_grid_X, position_grid_Y, position_grid_Z);
			double scatterOdd = SCATTER_FREQUENCY * (1 - fractionField) * velocityMag * velocityMag;//between 0~1 
			if(scatterOdd > ran_num) // alter the direction
			{
				//SCATTER_COEFFICIENT = 0.9;
				double x = 2 * ran_num * SCATTER_COEFFICIENT - SCATTER_COEFFICIENT + 1;
				double y = 2 * ran_num + SCATTER_COEFFICIENT - 1;
				double cosTheta;
				double theta;
				if(x == 0)
					theta = PI / 2.0f;
				else
				{
					cosTheta = y / x;
					theta = acos(cosTheta);
				}

				//if(x >=0 && y >= 0)
					//theta = acos(cosTheta);
				//else if(x < 0 && y >= 0)
				//	theta = acos(cosTheta);
				//else if(x >= 0 && y < 0)
				//	theta = acos(cosTheta) + PI;
				//else if(x < 0 && y < 0)
				//	theta = acos(cosTheta) + PI;
 				if(theta != 0)
				{
					//double ran_num2 = (((rand()%200) + 1) - 100) / 100.0f;
					double rotateAxisX = 1;
					double rotateAxisY = 0;
					double rotateAxisZ = 0;
					if(velocity[2] != 0)
					{
						rotateAxisX = (((rand()%200) + 1) - 100) / 100.0f;
						rotateAxisY = (((rand()%200) + 1) - 100) / 100.0f;
						while(rotateAxisX == 0 && rotateAxisY == 0)
						{
							rotateAxisX = (((rand()%200) + 1) - 100) / 100.0f;
							rotateAxisY = (((rand()%200) + 1) - 100) / 100.0f;
						}

						rotateAxisZ = -(rotateAxisX * velocity[0] + rotateAxisY * velocity[1]) / velocity[2];
					}
					else if(velocity[1] != 0)
					{
						rotateAxisX = (((rand()%200) + 1) - 100) / 100.0f;
						rotateAxisZ = (((rand()%200) + 1) - 100) / 100.0f;
						while(rotateAxisX == 0 && rotateAxisZ == 0)
						{
							rotateAxisX = (((rand()%200) + 1) - 100) / 100.0f;
							rotateAxisZ = (((rand()%200) + 1) - 100) / 100.0f;
						}
						rotateAxisY = -(rotateAxisX * velocity[0] + rotateAxisZ * velocity[2]) / velocity[1];
					}
					else if(velocity[0] != 0)
					{
						rotateAxisY = (((rand()%200) + 1) - 100) / 100.0f;
						rotateAxisZ = (((rand()%200) + 1) - 100) / 100.0f;
						while(rotateAxisY == 0 && rotateAxisZ == 0)
						{
							rotateAxisY = (((rand()%200) + 1) - 100) / 100.0f;
							rotateAxisZ = (((rand()%200) + 1) - 100) / 100.0f;
						}
						rotateAxisX = -(rotateAxisY * velocity[1] + rotateAxisZ * velocity[2]) / velocity[0];
					}
					double length = sqrt(rotateAxisX * rotateAxisX + rotateAxisY * rotateAxisY + rotateAxisZ * rotateAxisZ);
					rotateAxisX = rotateAxisX / length;
					rotateAxisY = rotateAxisY / length;
					rotateAxisZ = rotateAxisZ / length;

					double ss = cos(theta / 2.0f); 
					double xx = sin(theta / 2.0f) * rotateAxisX;
					double yy = sin(theta / 2.0f) * rotateAxisY;
					double zz = sin(theta / 2.0f) * rotateAxisZ;

					double newX = (1 - 2*yy*yy - 2*zz*zz)*velocity[0] + (2*xx*yy - 2*ss*zz)*velocity[1] + (2*xx*zz + 2*ss*yy)*velocity[2];
					double newY = (2*xx*yy + 2*ss*zz)*velocity[0] + (1 - 2*xx*xx - 2*zz*zz)*velocity[1] + (2*yy*zz - 2*ss*xx)*velocity[2];
					double newZ = (2*xx*zz - 2*ss*yy)*velocity[0] + (2*yy*zz + 2*ss*xx)*velocity[1] + (1 - 2*xx*xx - 2*yy*yy)*velocity[2];
					velocity[0] = newX;
					velocity[1] = newY;
					velocity[2] = newZ;
				}

			}

			bubbleVelList[j][i] = velocity + vec3(0, 1, 0) * (1 - ((rand()%10) + 1) / 100.0f);

			//add the terminal speed effect
			if(bubbleVelList[j][i].Length() > 5)//need to add the bubble radius parameter  5 is the terminal speed
			{
				bubbleVelList[j][i] = bubbleVelList[j][i] / bubbleVelList[j][i].Length() * 5;
			}

			position += dt * bubbleVelList[j][i];



			//Particle goes outside the container
			if(position[0] < 0 || position[0] > theDim[MACGrid::X] * theCellSize
			|| position[1] < 0 || position[1] >= theDim[MACGrid::Y] * theCellSize - 0.001
			|| position[2] < 0 || position[2] > theDim[MACGrid::Z] * theCellSize)
			{
				if(iterPos == bubblePosList[j].begin())
				{
					vec3 testA = *iterPos;
					bubblePosList[j].erase(iterPos);
					iterPos = bubblePosList[j].begin();
					vec3 testB = *iterPos;
					bubbleVelList[j].erase(iterVel);
					iterVel = bubbleVelList[j].begin();
					i--;
				}
				else
				{
					--(iterPos = bubblePosList[j].erase(iterPos));
					--(iterVel = bubbleVelList[j].erase(iterVel));
					i--;
					++iterPos; 
					++iterVel;
				}
			}
			else
			{
				bubblePosList[j][i] = position;


				//BUBBLE BREAK
				if(BUBBLE_BREAK_FREQ > (rand()%100) / 100.0f && j > 0)
				{

					vec3 newPosition1 = (position + vec3(bubbleRadiusList[j], 0, 0));
					vec3 newPosition2 = (position - vec3(bubbleRadiusList[j], 0, 0));
					if(newPosition1[0] < 0 || newPosition1[0] > theDim[MACGrid::X] * theCellSize)
						newPosition1 = position;
					if(newPosition2[0] < 0 || newPosition2[0] > theDim[MACGrid::X] * theCellSize)
						newPosition2 = position;	

					bubblePosList[j-1].push_back(newPosition1);
					bubbleVelList[j-1].push_back(bubbleVelList[j][i]);
					bubblePosList[j-1].push_back(newPosition2);
					bubbleVelList[j-1].push_back(bubbleVelList[j][i]);


					if(iterPos == bubblePosList[j].begin())
					{
						vec3 testA = *iterPos;
						bubblePosList[j].erase(iterPos);
						iterPos = bubblePosList[j].begin();
						vec3 testB = *iterPos;
						bubbleVelList[j].erase(iterVel);
						iterVel = bubbleVelList[j].begin();
						i--;
					}
					else
					{
						--(iterPos = bubblePosList[j].erase(iterPos));
						--(iterVel = bubbleVelList[j].erase(iterVel));
						i--;
						++iterPos; 
						++iterVel;
					}
				}
				else
				{
					++iterPos; 
					++iterVel;
				}
			}

		}
	}
}

void MACGrid::doSimulation(int frame)
{
	if(frame < currentFrame)
	{
		//RESET
		reset();
		currentFrame = 0;
	}

	//Do simulation
	for(int i = 0; i < frame - currentFrame; i++)
	{
		// Step0: Gather user forces
		generateBubbles();
		updateSources();

		// Step1: Calculate new velocities
		//advectVelocity(TIME_STEP);
		//addExternalForces(TIME_STEP);
		//project(TIME_STEP);

		// Step2: Calculate new temperature
		//advectTemperature(TIME_STEP);

		// Step3: Calculate new density 
		//advectDensity(TIME_STEP);
		computeFractionField();
		advectBubbles(TIME_STEP);
	}

	currentFrame = frame;
}

float* MACGrid::getBubblePosition(int index, int* size)
{
	float* bubblePositionArray = new float[bubblePosList[index].size() * 3];
	*size = bubblePosList[index].size();
	for(int i = 0; i < bubblePosList[index].size(); i++)
	{
		vec3 position = bubblePosList[index][i];
		bubblePositionArray[i * 3] = position[0];
		bubblePositionArray[i * 3 + 1] = position[1];
		bubblePositionArray[i * 3 + 2] = position[2];
	}

	return bubblePositionArray;
}

void MACGrid::setViscosity(double viscosity)
{
	FLUID_VISCOSITY = viscosity;
}

void MACGrid::setDensity(double density)
{
	FLUID_DENSITY = density;
}

void MACGrid::setScatterFreq(double freq)
{
	SCATTER_FREQUENCY = freq;
}

void MACGrid::setScatterCoef(double coef)
{
	SCATTER_COEFFICIENT = coef;
}

void MACGrid::setSourcePos(std::vector<vec3> pos)
{
	sourcePos.clear();
	sourcePos = pos;
}

void MACGrid::setTimeStep(double dt)
{
	TIME_STEP = dt;
}

void MACGrid::setBubbleBreakFreq(double freq)
{
	BUBBLE_BREAK_FREQ = freq;
}

void MACGrid::setBubbleRadius(double radiusMin, double radiusMax, int segment)
{
	bubbleRadiusList.clear();

	if(radiusMin >radiusMax)
	{
		double temp = radiusMin;
		radiusMin = radiusMax;
		radiusMax = temp;
	}

	//if(radiusMin == radiusMax)
	//	segment = 1;

	segment = 1;
	while(radiusMax / radiusMin >= 2 )
	{
		radiusMax = radiusMax / 2;
		segment++;
	}


	//if(segment <= 1)
	//{
	//	double radius = (radiusMax +radiusMin) / 2.0f;
	//	bubbleRadiusList.push_back(radius);
	//}
	//else
	//{
		for(int i = 0; i < segment ; i++)
		{

			//double radius = radiusMin * pow((radiusMax / radiusMin), (double)i / (segment - 1));
			double radius = radiusMin * pow(pow(2, 1.0f / 3.0f), i);
			bubbleRadiusList.push_back(radius);
		}
	//}
}

int MACGrid::getBubbleRadiusCount()
{
	return bubbleRadiusList.size();
}

double MACGrid::getBubbleRadius(int index)
{
	return bubbleRadiusList[index];
}


void MACGrid::advectVelocity(double dt)
{
    // TODO: compute new velocities and store in target


	// loops through every MacGrid face; defines i, j, k
	FOR_EACH_FACE {
		// mU = GridDataX, mV = GridDataY, mW = GridDataZ
		double velU = mU( i, j, k );
		double velV = mV( i, j, k );
		double velW = mW( i, j, k );

		// old, seemingly incorrect, method
		// compute velocity gradient - rate of change in each component direction
		//vec3 velocityGradient( ( mU( i+1, j, k ) - mU( i, j, k ) ) / theCellSize,
		//					   ( mV( i, j+1, k ) - mV( i, j, k ) ) / theCellSize,
		//					   ( mW( i, j, k+1 ) - mW( i, j, k ) ) / theCellSize );
		// solve for advection
		//velU = velU + ( dt * -1.0f * velU * Dot( UNIT_X, velocityGradient ) );
		//velV = velV + ( dt * -1.0f * velV * Dot( UNIT_Y, velocityGradient ) );
		//velW = velW + ( dt * -1.0f * velW * Dot( UNIT_Z, velocityGradient ) );

		// compute face center positions
		vec3 centerPosition = getCenter( i, j, k );
		vec3 grid_x_bottom_border_pos = centerPosition - vec3( theCellSize * 0.5f, 0.0f, 0.0f );
		vec3 grid_y_bottom_border_pos = centerPosition - vec3( 0.0f, theCellSize * 0.5f, 0.0f );
		vec3 grid_z_bottom_border_pos = centerPosition - vec3( 0.0f, 0.0f, theCellSize * 0.5f );


		////////////////////////////////////////////////////
		// compute velocities at face centers
		////////////////////////////////////////////////////

		vec3 grid_x_bottom_border_vel, grid_y_bottom_border_vel, grid_z_bottom_border_vel;

		// compute velocity at grid_x_bottom_border_pos
		if ( i == 0 ) {
			// low boundary cell
			grid_x_bottom_border_vel[VX] = velU;
			grid_x_bottom_border_vel[VY] = 0.5f * ( velV + mV( i, j+1, k ) );
			grid_x_bottom_border_vel[VZ] = 0.5f * ( velW + mW( i, j, k+1 ) );
		}
		else if ( i == theDim[MACGrid::X] ) {
			// high boundary cell
			grid_x_bottom_border_vel[VX] = velU;
			grid_x_bottom_border_vel[VY] = 0.5f * ( mV( i-1, j, k ) + mV( i-1, j+1, k ) );
			grid_x_bottom_border_vel[VZ] = 0.5f * ( mW( i-1, j, k ) + mW( i-1, j, k+1 ) );
		}
		else {
			// not boundary cell - cell is somehwere in middle of container
			grid_x_bottom_border_vel[VX] = velU;
			grid_x_bottom_border_vel[VY] = 0.25f * ( velV + mV( i-1, j, k ) + mV( i, j+1, k ) + mV( i-1, j+1, k ) );
			grid_x_bottom_border_vel[VZ] = 0.25f * ( velW + mW( i-1, j, k ) + mW( i, j, k+1 ) + mW( i-1, j, k+1 ) );
		}

		// compute velocity at grid_y_bottom_border_pos
		if ( j == 0 ) {
			// low boundary cell
			grid_y_bottom_border_vel[VX] = 0.5f * ( velU + mU( i+1, j, k ) );
			grid_y_bottom_border_vel[VY] = velV;
			grid_y_bottom_border_vel[VZ] = 0.5f * ( velW + mW( i, j, k+1 ) );
		}
		else if ( j == theDim[MACGrid::Y] ) {
			// high boundary cell
			grid_y_bottom_border_vel[VX] = 0.5f * (  mU( i, j-1, k ) + mU( i+1, j-1, k ) );
			grid_y_bottom_border_vel[VY] = velV;
			grid_y_bottom_border_vel[VZ] = 0.5f * ( mW( i, j-1, k ) + mW( i, j-1, k+1 ) );
		}
		else {
			// not boundary cell - cell is somehwere in middle of container
			grid_y_bottom_border_vel[VX] = 0.25f * ( velU + mU( i, j-1, k ) + mU( i+1, j, k ) + mU( i+1, j-1, k ) );
			grid_y_bottom_border_vel[VY] = velV;
			grid_y_bottom_border_vel[VZ] = 0.25f * ( velW + mW( i, j-1, k ) + mW( i, j, k+1 ) + mW( i, j-1, k+1 ) );
		}

		// compute velocity at grid_z_bottom_border_pos
		if ( k == 0 ) {
			// low boundary cell
			grid_z_bottom_border_vel[VX] = 0.5f * ( velU + mU( i+1, j, k ) );
			grid_z_bottom_border_vel[VY] = 0.5f * ( velV + mV( i, j+1, k ) );
			grid_z_bottom_border_vel[VZ] = velW;
		}
		else if ( k == theDim[MACGrid::Z] ) {
			// high boundary cell
			grid_z_bottom_border_vel[VX] = 0.5f * ( mU( i, j, k-1 ) + mU( i+1, j, k-1 ) );
			grid_z_bottom_border_vel[VY] = 0.5f * ( mV( i, j, k-1 ) + mV( i, j+1, k-1 ) );
			grid_z_bottom_border_vel[VZ] = velW;
		}
		else {
			// not boundary cell - cell is somehwere in middle of container
			grid_z_bottom_border_vel[VX] = 0.25f * ( velU + mU( i, j, k-1 ) + mU( i+1, j, k ) + mU( i+1, j, k-1 ) );
			grid_z_bottom_border_vel[VY] = 0.25f * ( velV + mV( i, j, k-1 ) + mV( i, j+1, k ) + mV( i, j+1, k-1 ) );
			grid_z_bottom_border_vel[VZ] = velW;
		}

		// solve for advection
		velU = getVelocityX( grid_x_bottom_border_pos - ( dt * grid_x_bottom_border_vel ) );
		velV = getVelocityY( grid_y_bottom_border_pos - ( dt * grid_y_bottom_border_vel ) );
		velW = getVelocityZ( grid_z_bottom_border_pos - ( dt * grid_z_bottom_border_vel ) );

		// store in target, checking boundaries
		if ( i > 0 && i < theDim[MACGrid::X] ) {
			target.mU( i, j, k ) = velU;
		}
		else {
			target.mU( i, j, k ) = 0.0f;
		}

		if ( j > 0 && j < theDim[MACGrid::Y] ) {
			target.mV( i, j, k ) = velV;
		}
		else {
			target.mV( i, j, k ) = 0.0f;
		}

		if ( k != 0 && k != theDim[MACGrid::Z] ) {
			target.mW( i, j, k ) = velW;
		}
		else {
			target.mW( i, j, k ) = 0.0f;
		}
	}

    // save result to object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature( double dt )
{
    // TODO: calculate new temp and store in target


	FOR_EACH_CELL {
		// velocity at cell center
		vec3 velocity( (mU( i, j, k ) + mU( i+1, j, k ))/2.0f, (mV( i, j, k ) + mV( i, j+1, k ))/2.0f, (mW( i, j, k )+mW( i, j, k+1 ))/2.0f );

		// trace back particle position using known velocity
		vec3 pos = getCenter( i, j, k );
		pos -= dt * velocity;

		// interpolate temperature for passed-in position, and store in target
		target.mT( i,j,k ) = getTemperature( pos );
	}

	// save result to object
	mT = target.mT;
}

void MACGrid::advectDensity( double dt )
{
    // TODO: calculate new densitities and store in target

	// use an identical trace back method to the one used in MACGrid::advectTemperature()
	FOR_EACH_CELL {
		vec3 velocity( (mU( i, j, k ) + mU( i+1, j, k ))/2.0f, (mV( i, j, k ) + mV( i, j+1, k ))/2.0f, (mW( i, j, k )+mW( i, j, k+1 ))/2.0f );
		vec3 pos = getCenter(i, j, k);
		pos -= dt * velocity;
		target.mD(i,j,k) = getDensity(pos);
	}
	mD = target.mD;
}

void MACGrid::computeBouyancy( double dt )
{
	// TODO: calculate bouyancy and store in target
	// TODO: tune alpha and beta parameters

	double alpha = 0.1f;
	double beta = 0.01f;
	double ambient_temp = 270.0f;

	FOR_EACH_CELL {
		target.mTemp( i, j, k ) = -1.0f * alpha * mD( i, j, k ) + beta * ( mT( i, j, k )  - ambient_temp );
	}

	FOR_EACH_CELL {
		if ( j != 0 ) {
			// get cell center position
			vec3 cell_center_pos = getCenter( i, j, k );

			// get y face center position
			vec3 y_face_center_pos = cell_center_pos - vec3( 0.0f, theCellSize * 0.5f, 0.0f );

			// interpolate buoyancy force from cell center to face center
			double buoyancy_force = target.mTemp.interpolate( y_face_center_pos );
		
			// acceleration = force on particle / mass of particle
			double acceleration = buoyancy_force / PARTICLE_MASS;

			// perform explicit Euler integration to update existing velocity with applied force
			// v' = v + at
			target.mV( i, j, k ) = mV( i, j, k ) + dt * acceleration;

			// this is incorrect
			//target.mV( i, j, k ) = mV( i, j, k ) + ( target.mTemp( i, j, k ) + target.mTemp( i, j-1, k ) ) / 2.0f;
		}
	}

	mV = target.mV;
}

void MACGrid::computeVorticityConfinement( double dt )
{
	// TODO: calculate vorticity confinement forces


	double epsilon = 1.0f;
	double two_cell_widths = 2.0f * theCellSize;
	double very_small = pow( 10.0f, -20.0f );

	FOR_EACH_CELL {

		// TODO: ask about boundary conditions here
		// index out of bounds when i == 0, 1, theDim[MACGrid::X]
		// index out of bounds when j == 0, 1, theDim[MACGrid::Y]
		// index out of bounds when k == 0, 1, theDim[MACGrid::Z]
		// velocities return 0 by default for these cases
		vec3 omegaGradient( ( getOmegaVector( i+1, j, k ).Length() - getOmegaVector( i-1, j, k ).Length() ) / two_cell_widths,
							( getOmegaVector( i, j+1, k ).Length() - getOmegaVector( i, j-1, k ).Length() ) / two_cell_widths,
							( getOmegaVector( i, j, k+1 ).Length() - getOmegaVector( i, j, k-1 ).Length() ) / two_cell_widths );
	
		// add very_small to prevent divide by zero
		vec3 normal = omegaGradient / ( omegaGradient.Length() + very_small );

		//vec3 temp = getOmegaVector( i, j, k );
		vec3 confinement = epsilon * theCellSize * normal.Cross( getOmegaVector( i, j, k ) );
		target.mConfForceX(i, j, k) = confinement[0];
		target.mConfForceY(i, j, k) = confinement[1];
		target.mConfForceZ(i, j, k) = confinement[2];
	}
	
	// TODO: ask how to apply computed forces to velocity field

	// take vorticity confinement forces computed at cell centers and approximate at faces for velocity field 
	FOR_EACH_CELL {

		// TODO: ask about boundary conditions here
		// currently, faces with 0 indices are being ignored

		// first, approximate vorticity force at face center with linear interpolation
		// next, acceleration = force on particle / mass of particle
		// finally, perform explicit Euler integration to update existing velocity with applied force
		// v' = v + at
		if ( i != 0 ) {
			double vorticity_force_x = ( target.mConfForceX( i, j, k ) + target.mConfForceX( i-1, j, k ) ) / 2.0f;
			double acceleration_x = vorticity_force_x / PARTICLE_MASS;
			target.mU( i, j, k ) = mU( i, j, k ) + dt * acceleration_x;
		}
		if ( j != 0 ) {
			double vorticity_force_y = ( target.mConfForceY( i, j, k ) + target.mConfForceY( i, j-1, k ) ) / 2.0f;
			double acceleration_y = vorticity_force_y / PARTICLE_MASS;
			target.mV( i, j, k ) = mV( i, j, k ) + dt * acceleration_y;
		}
		if ( k != 0 ) {
			double vorticity_force_z = ( target.mConfForceZ( i, j, k ) + target.mConfForceZ( i, j, k-1 ) ) / 2.0f;
			double acceleration_z = vorticity_force_z / PARTICLE_MASS;
			target.mW( i, j, k ) = mW( i, j, k ) + dt * acceleration_z;	
		}
	}

	// save result to object
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::computeViscosityForce(double dt)
{
	FOR_EACH_CELL 
	{
		//velocity difference between two grids

		double velDiff_IPlus1_I  = mU( i+1, j, k ) - mU( i, j, k );
		double velDiff_I_IMinus1 = mU( i, j, k ) - mU( i-1, j, k );
		double velDiff_JPlus1_J  = mU( i, j+1, k ) - mU( i, j, k );
		double velDiff_J_JMinus1 = mU( i, j, k ) - mU( i, j-1, k );
		double velDiff_KPlus1_K  = mU( i, j, k+1 ) - mU( i, j, k );
		double velDiff_K_KMinus1 = mU( i, j, k ) - mU( i, j, k-1 );

		if( i > 0 )
			velDiff_I_IMinus1 = mU( i, j, k ) - mU( i-1, j, k );
		else 
			velDiff_I_IMinus1 = 0.0f;

		if( i < theDim[MACGrid::X])
			velDiff_IPlus1_I  = mU( i+1, j, k ) - mU( i, j, k );
		else
			velDiff_IPlus1_I = 0.0f;

		if( j > 0 )
			velDiff_J_JMinus1  = mV( i, j, k ) - mV( i, j-1, k );
		else
			velDiff_J_JMinus1  = 0.0f;

		if( i < theDim[MACGrid::Y])
			velDiff_JPlus1_J  = mV( i, j+1, k ) - mV( i, j, k );
		else
			velDiff_JPlus1_J = 0.0f;

		if( k > 0 )
			velDiff_K_KMinus1 = mW( i, j, k ) - mW( i, j, k-1 );
		else
			velDiff_K_KMinus1 = 0.0f;

		if( k < theDim[MACGrid::Z] )
			velDiff_KPlus1_K = mW( i, j, k+1 ) - mW( i, j, k );
		else
			velDiff_KPlus1_K = 0.0f;

		// acceleration = force on particle / mass of particle
		// finally, perform explicit Euler integration to update existing velocity with applied force
		// v' = v + at
		if ( i != 0 ) {
			double viscosityForceX = FLUID_VISCOSITY * (velDiff_IPlus1_I - velDiff_I_IMinus1) / theCellSize / theCellSize;
			double acceleration_x = viscosityForceX / PARTICLE_MASS;
			target.mU( i, j, k ) = mU( i, j, k ) + dt * acceleration_x;
		}
		if ( j != 0 ) {
			double viscosityForceY = FLUID_VISCOSITY * (velDiff_JPlus1_J - velDiff_J_JMinus1) / theCellSize / theCellSize;
			double acceleration_y = viscosityForceY / PARTICLE_MASS;
			target.mV( i, j, k ) = mV( i, j, k ) + dt * acceleration_y;
		}
		if ( k != 0 ) {
			double viscosityForceZ = FLUID_VISCOSITY * (velDiff_KPlus1_K - velDiff_K_KMinus1) / theCellSize / theCellSize;
			double acceleration_z = viscosityForceZ / PARTICLE_MASS;
			target.mW( i, j, k ) = mW( i, j, k ) + dt * acceleration_z;	
		}

	}

	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

vec3 MACGrid::getOmegaVector( int i, int j, int k )
{
	double two_cell_widths = 2.0f * theCellSize;

	// return omega vector
	return vec3( ( mW( i, j+1, k ) - mW( i, j-1, k ) ) / two_cell_widths - ( mV( i, j, k+1 ) - mV( i, j, k-1 ) ) / two_cell_widths,
				 ( mU( i, j, k+1 ) - mU( i, j, k-1 ) ) / two_cell_widths - ( mW( i+1, j, k ) - mW( i-1, j, k ) ) / two_cell_widths,
				 ( mV( i+1, j, k ) - mV( i-1, j, k ) ) / two_cell_widths - ( mU( i, j+1, k ) - mU( i, j-1, k ) ) / two_cell_widths );
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
   //computeViscosityForce(dt);
}

void MACGrid::project( double dt )
{
	// TODO: solve Ap = d for pressure

	// TODO: IMPLEMENT assert( checkDivergence() ) AS A SANITY CHECK


	////////////////////////////////////////////////////
	// solve for pressure
	////////////////////////////////////////////////////

	double constant = -1.0f * theCellSize * theCellSize / dt;

	// pressure coefficient matrix
	GridDataMatrix A;

	// d is divergence vector; p is pressure vector
	GridData d, p;
	d.initialize( 0.0f );
	p.initialize( 0.0f );

	FOR_EACH_CELL {
		// fill divergence vector - ( constant * density * velocity gradient )
		d(i, j, k) = constant * FLUID_DENSITY * ( ( mU( i+1, j, k ) - mU( i, j, k ) ) / theCellSize + 
												  ( mV( i, j+1, k ) - mV( i, j, k ) ) / theCellSize +
												  ( mW( i, j, k+1 ) - mW( i, j, k ) ) / theCellSize );

		int num_neighbors = 6;
		bool has_neighbor_plus_x = true, has_neighbor_plus_y = true, has_neighbor_plus_z = true;

		// count neighbors, and determine whether neighbors exist along positive directions
		if ( i <= 0 ) {
			--num_neighbors;
		}
		if ( j <= 0 ) {
			--num_neighbors;
		}
		if ( k <= 0 ) {
			--num_neighbors;
		}
		if ( i+1 >= theDim[MACGrid::X] ) {
			--num_neighbors;
			has_neighbor_plus_x = false;
		}
		if ( j+1 >= theDim[MACGrid::Y] ) {
			--num_neighbors;
			has_neighbor_plus_y = false;
		}
		if ( k+1 >= theDim[MACGrid::Z] ) {
			--num_neighbors;
			has_neighbor_plus_z = false;
		}

		// set A.diag - number of neighbors the current cell has
		A.diag( i, j, k ) = num_neighbors;

		// set A.plusI - neighbor cell in positive x direction
		if ( has_neighbor_plus_x ) {
			A.plusI( i, j, k ) = -1.0f;
		}
		else {
			A.plusI( i, j, k ) = 0.0f;
		}

		// set A.plusJ - neighbor cell in positive y direction
		if ( has_neighbor_plus_y ) {
			A.plusJ( i, j, k ) = -1.0f;
		}
		else {
			A.plusJ( i, j, k ) = 0.0f;
		}

		// set A.plusK - neighbor cell in positive z direction
		if ( has_neighbor_plus_z ) {
			A.plusK( i, j, k ) = -1.0f;
		}
		else {
			A.plusK( i, j, k ) = 0.0f;
		}
	}

	// solve pressure vector by solving system of linear equations
	int max_num_iterations = 100;
	double tolerance = 0.00001f;
	conjugateGradient( A, p, d, max_num_iterations, tolerance );

	// store in target
	target.mP = p;
	

	////////////////////////////////////////////////////
	// apply computed pressures to velocity field
	////////////////////////////////////////////////////

	FOR_EACH_FACE {
		double velU = mU( i, j, k );
		double velV = mV( i, j, k );
		double velW = mW( i, j, k );

		double current_pressure = target.mP( i, j, k );

		// to check boundary conditions
		double pressure_i_minus_1, pressure_j_minus_1, pressure_k_minus_1;

		// set pressure_i_minus_1
		if ( i-1 < 0 ) {
			pressure_i_minus_1 = 0.0f;
		}
		else {
			pressure_i_minus_1 = target.mP( i-1, j, k );
		}

		// set pressure_j_minus_1
		if ( j-1 < 0 ) {
			pressure_j_minus_1 = 0.0f;
		}
		else {
			pressure_j_minus_1 = target.mP( i, j-1, k );
		}

		// set pressure_k_minus_1
		if ( k-1 < 0 ) {
			pressure_k_minus_1 = 0.0f;
		}
		else {
			pressure_k_minus_1 = target.mP( i, j, k-1 );
		}

		// apply computed pressures to velocity field
		// if face is an outside border of container, then velocity at that face must be 0
		if ( i == 0 || i == theDim[MACGrid::X] ) {
			velU = 0.0f;
		}
		else {
			velU -= dt * ( 1.0f / FLUID_DENSITY ) * ( ( target.mP( i, j, k ) - pressure_i_minus_1 ) / theCellSize );
		}
		
		if ( j == 0 || j == theDim[MACGrid::Y] ) {
			velV = 0.0f;
		}
		else {		
			velV -= dt * ( 1.0f / FLUID_DENSITY ) * ( ( target.mP( i, j, k ) - pressure_j_minus_1 ) / theCellSize );
		}
		
		if ( k == 0 || k == theDim[MACGrid::Z] ) {
			velW = 0.0f;
		}
		else {
			velW -= dt * ( 1.0f / FLUID_DENSITY ) * ( ( target.mP( i, j, k ) - pressure_k_minus_1 ) / theCellSize );
		}

		// store in target with additional container boundary checking
		// ex: for a 2x2 matrix, the j value for mU can be only 0 or 1, but j can have values of 0, 1, or 2 for a 2x2 matrix
		if ( j < theDim[MACGrid::Y] && k < theDim[MACGrid::Z] ) {
			target.mU( i, j, k ) = velU;
		}
		if ( i < theDim[MACGrid::X] && k < theDim[MACGrid::Z] ) {
			target.mV( i, j, k ) = velV;
		}
		if ( i < theDim[MACGrid::X] && j < theDim[MACGrid::Y] ) {
			target.mW( i, j, k ) = velW;
		}
	}
	
	// save result to object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;

	// debug - test if system is divergence free
	FOR_EACH_CELL {
		double sum = (mU(i+1,j,k) - mU(i,j,k)) + 
					 (mV(i,j+1,k) - mV(i,j,k)) +
					 (mW(i,j,k+1) - mW(i,j,k));


		if ( abs( sum ) > 0.01 ) {
			bool non_divergence_free = true;
		}
	}
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	//FOR_EACH_CELL {

	//	int numFluidNeighbors = 0;
	//	if (i-1 >= 0) {
	//		AMatrix.plusI(i-1,j,k) = -1;
	//		numFluidNeighbors++;
	//	}
	//	if (i+1 < theDim[MACGrid::X]) {
	//		AMatrix.plusI(i,j,k) = -1;
	//		numFluidNeighbors++;
	//	}
	//	if (j-1 >= 0) {
	//		AMatrix.plusJ(i,j-1,k) = -1;
	//		numFluidNeighbors++;
	//	}
	//	if (j+1 < theDim[MACGrid::Y]) {
	//		AMatrix.plusJ(i,j,k) = -1;
	//		numFluidNeighbors++;
	//	}
	//	if (k-1 >= 0) {
	//		AMatrix.plusK(i,j,k-1) = -1;
	//		numFluidNeighbors++;
	//	}
	//	if (k+1 < theDim[MACGrid::Z]) {
	//		AMatrix.plusK(i,j,k) = -1;
	//		numFluidNeighbors++;
	//	}
	//	// Set the diagonal:
	//	AMatrix.diag(i,j,k) = numFluidNeighbors;
	//}
}



/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}



