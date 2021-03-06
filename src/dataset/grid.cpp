
/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Grid: Irregular, Curvilinear, Cartesian grid
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cfloat>
#include "grid.h"
#include "core/interpolator.h"

using namespace std;


namespace edda {


//////////////////////////////////////////////////////////////////////////
//
//	definition of Cartesian Grid Class
//
//////////////////////////////////////////////////////////////////////////

CartesianGrid::CartesianGrid(int xdim, int ydim, int zdim)
{
  m_nDimension[0] = xdim;  //the grid dimensions in C space 
  m_nDimension[1] = ydim;
  m_nDimension[2] = zdim;
}

CartesianGrid::CartesianGrid()
{
  Reset();
}

CartesianGrid::~CartesianGrid()
{
}

void CartesianGrid::Reset()
{
  m_nDimension[0] = m_nDimension[1] = m_nDimension[2] = 0;
}

int *CartesianGrid::getDimension()
{
  return m_nDimension;
}

ReturnStatus CartesianGrid::getIndex(int i, int j, int k, size_t &idx)
{
  if (i>=0 && j>=0 && k>=0 && i<m_nDimension[0] && j<m_nDimension[1] && k<m_nDimension[2]) {
    idx = i+m_nDimension[0]*(j+m_nDimension[1]*(size_t)k);
    return SUCCESS;
  }
  return OUT_OF_BOUND;
}

//////////////////////////////////////////////////////////////////////////
//
//
//	definition of Regular Cartesian Grid Class
//
//
//////////////////////////////////////////////////////////////////////////

// constructor and deconstructor
RegularCartesianGrid::RegularCartesianGrid()
    :CartesianGrid(0,0,0)
{
  Reset();
}

RegularCartesianGrid::RegularCartesianGrid(int xdim, int ydim, int zdim)
    :CartesianGrid(xdim, ydim, zdim)
{
  VECTOR3 a(0,0,0);   // the default is from 0 to xdim-1, etc.
  VECTOR3 b(xdim-1, ydim-1, zdim-1);
  setBoundary(a, b);
}

RegularCartesianGrid::~RegularCartesianGrid()
{
}

void RegularCartesianGrid::reset(void)
{
  mappingFactorX = mappingFactorY = mappingFactorZ = 0.0;
  oneOvermappingFactorX = oneOvermappingFactorY = oneOvermappingFactorZ = 0.0;
  gridSpacing = 1.0;
}

void RegularCartesianGrid::boundary(VECTOR3& minB, VECTOR3& maxB)
{
  minB = m_vMinBound;
  maxB = m_vMaxBound;
}

//////////////////////////////////////////////////////////////////////////
// set bounding box
//////////////////////////////////////////////////////////////////////////
void RegularCartesianGrid::setBoundary(VECTOR3& minB, VECTOR3& maxB)
{
  m_vMinBound = minB;
  m_vMaxBound = maxB;
  m_vMinRealBound.set(m_vMinBound[0], m_vMinBound[1], m_vMinBound[2], 0.0);
  m_vMaxRealBound.set(m_vMaxBound[0], m_vMaxBound[1], m_vMaxBound[2], 0.0);
  mappingFactorX = (float)(xdim()-1)/(m_vMaxBound[0] - m_vMinBound[0]);
  mappingFactorY = (float)(ydim()-1)/(m_vMaxBound[1] - m_vMinBound[1]);
  mappingFactorZ = (float)(zdim()-1)/(m_vMaxBound[2] - m_vMinBound[2]);
  oneOvermappingFactorX = (m_vMaxBound[0] - m_vMinBound[0])/(float)(xdim()-1);
  oneOvermappingFactorY = (m_vMaxBound[1] - m_vMinBound[1])/(float)(ydim()-1);
  oneOvermappingFactorZ = (m_vMaxBound[2] - m_vMinBound[2])/(float)(zdim()-1);
 
  // grid spacing
  gridSpacing = min(min(oneOvermappingFactorX, oneOvermappingFactorY), 
		    oneOvermappingFactorZ);
}

//////////////////////////////////////////////////////////////////////////
// set real bounding box (does not include ghost cells)
//////////////////////////////////////////////////////////////////////////
void RegularCartesianGrid::setRealBoundary(VECTOR4& minB, VECTOR4& maxB)
{
  m_vMinRealBound = minB;
  m_vMaxRealBound = maxB;
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point pos is in bounding box (counting ghost cells)
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::isInBBox(VECTOR3& pos)
{
  if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
      (pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
      (pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]))
    return true;
  else
    return false;
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point pos is in bounding box (not counting ghost cells)
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::isInRealBBox(VECTOR3& pos)
{
  if( (pos[0] >= m_vMinRealBound[0]) && (pos[0] <= m_vMaxRealBound[0]) &&
      (pos[1] >= m_vMinRealBound[1]) && (pos[1] <= m_vMaxRealBound[1]) &&
      (pos[2] >= m_vMinRealBound[2]) && (pos[2] <= m_vMaxRealBound[2]))
    return true;
  else
    return false;
}

bool RegularCartesianGrid::isInRealBBox(VECTOR3& pos, float t)
{
  if( (pos[0] >= m_vMinRealBound[0]) && (pos[0] <= m_vMaxRealBound[0]) &&
      (pos[1] >= m_vMinRealBound[1]) && (pos[1] <= m_vMaxRealBound[1]) &&
      (pos[2] >= m_vMinRealBound[2]) && (pos[2] <= m_vMaxRealBound[2]) &&
	  (     t >= m_vMinRealBound[3]) && (     t <= m_vMaxRealBound[3]) )
    return true;
  else
    return false;
}

// compute a default boundary 
void RegularCartesianGrid::computeBBox(void)
{
  VECTOR3 minB(0, 0, 0);  // default is from zero to xdim-1, etc.
  VECTOR3 maxB((float)(xdim()-1), (float)(ydim()-1), (float)(zdim()-1));

  setBoundary(minB, maxB);
}

//////////////////////////////////////////////////////////////////////////
// for Cartesian grid, this funcion means whether the physical point is 
// in the boundary
//////////////////////////////////////////////////////////////////////////
#if 0
ReturnStatus at_phys(VECTOR3& pos, AbstractDataArray *array, boost::any output)
{
	// whether in the bounding box
    if (!isInBBox(pos))
        return OUT_OF_BOUND;
}
#endif

//////////////////////////////////////////////////////////////////////////
// get vertex list of a cell
// input
//		cellId:		cell Id
//		cellType:	cell type
// output
//		vVertices: the vertex lis of the cell
//////////////////////////////////////////////////////////////////////////
ReturnStatus  RegularCartesianGrid::getCellVertices(int cellId,
                                          std::vector<size_t> &vVertices)
{
  int totalCell = xcelldim() * ycelldim() * zcelldim();
  size_t xidx, yidx, zidx, index;

  if((cellId < 0) || (cellId >= totalCell))
    return FAIL;

  vVertices.clear();
  zidx = cellId / (xcelldim() * ycelldim());
  yidx = (cellId % (xcelldim() * ycelldim())) / xcelldim();
  xidx = cellId - zidx * xcelldim() * ycelldim() - yidx * xcelldim();

  for(int kFor = 0; kFor < 2; kFor++)
    for(int jFor = 0; jFor < 2; jFor++)
      for(int iFor = 0; iFor < 2; iFor++)
	{
	  index = (zidx+kFor) * ydim() * xdim() + (yidx + jFor) * xdim() + (xidx + iFor);
	  vVertices.push_back(index);
	}
  return SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// get the physical coordinate of the vertex
//
// input:
// verIdx: index of vertex
// output:
// pos: physical coordinate of vertex
//////////////////////////////////////////////////////////////////////////
ReturnStatus RegularCartesianGrid::at_vertex(int verIdx, VECTOR3& pos)
{
  int xidx, yidx, zidx;
  int w=xdim(), h=ydim(), d=zdim();
  int totalVer = w*h*d;
  if((verIdx < 0) || (verIdx >= totalVer))
    return FAIL;

  zidx = verIdx / (w*h);
  yidx = (verIdx % (w*h)) / w; // bug fixed 07/04/2014
  xidx = verIdx % w;

  float xpos = m_vMinBound[0] + xidx*oneOvermappingFactorX; 
  float ypos = m_vMinBound[1] + yidx*oneOvermappingFactorY; 
  float zpos = m_vMinBound[2] + zidx*oneOvermappingFactorZ; 

  // pos.set((float)xidx, (float)yidx, (float)zidx);
  pos.set((float)xpos, (float)ypos, (float)zpos);
  return SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// whether the point in the physical position is in the cell
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// output:
// pInfo.interpolant: interpolation coefficient
// return:		returns 1 if in cell
//////////////////////////////////////////////////////////////////////////
bool RegularCartesianGrid::isInCell(PointInfo& pInfo, const int cellId)
{
  if(!isInBBox(pInfo.phyCoord))
    return false;

  float cx, cy, cz; // computatnoial space x, y, and z 
  cx = (pInfo.phyCoord[0]-m_vMinBound[0])/oneOvermappingFactorX; 
  cy = (pInfo.phyCoord[1]-m_vMinBound[1])/oneOvermappingFactorY; 
  cz = (pInfo.phyCoord[2]-m_vMinBound[2])/oneOvermappingFactorZ; 

  int xidx, yidx, zidx;
  xidx = (int)floor(cx); 
  yidx = (int)floor(cy); 
  zidx = (int)floor(cz); 

  int inCell = zidx * ycelldim() * xcelldim() + yidx * xcelldim() + xidx;
  if(cellId == inCell)
    {
      pInfo.interpolant.set(cx - (float)xidx, cy - (float)yidx, cz - (float)zidx);
      return true;
    }
  else
    return true;
}

//////////////////////////////////////////////////////////////////////////
// get the cell id and also interpolating coefficients for the given
// physical position
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// fromCell:	if -1, initial point; otherwise this point is advected from others
// output:
// pInfo.inCell: in which cell
// pInfo.interpolant: interpolation coefficients
// return:	returns 1 if successful; otherwise returns -1
//////////////////////////////////////////////////////////////////////////
ReturnStatus RegularCartesianGrid::phys_to_cell(PointInfo& pInfo)
{
  if(!isInBBox(pInfo.phyCoord))
    return FAIL;

  float cx, cy, cz; // computatnoial space x, y, and z 
  cx = (pInfo.phyCoord[0]-m_vMinBound[0])/oneOvermappingFactorX; 
  cy = (pInfo.phyCoord[1]-m_vMinBound[1])/oneOvermappingFactorY; 
  cz = (pInfo.phyCoord[2]-m_vMinBound[2])/oneOvermappingFactorZ; 

  int xidx, yidx, zidx;
  xidx = (int)floor(cx); 
  yidx = (int)floor(cy);
  zidx = (int)floor(cz); 

  // special treatment for the very boundary position -- fixed by Jimmy 8/30/2014
  if (xidx == xcelldim()) xidx--;
  if (yidx == ycelldim()) yidx--;
  if (zidx == zcelldim()) zidx--;

  pInfo.inCell = zidx * ycelldim() * xcelldim() + yidx * xcelldim() + xidx;
  pInfo.interpolant.set(cx - (float)xidx, cy - (float)yidx, cz - (float)zidx);

  //std::cout << pInfo.inCell << ',' << pInfo.interpolant << endl;
  return SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// trilinear interpolation
// input
// nodeData:	8 corners of cube cell
// coeff:	bilinear interpolation coefficents
// output
// vData:	output
//////////////////////////////////////////////////////////////////////////
#if 0
void RegularCartesianGrid::interpolate(vector<boost::any>& vData,
                        VECTOR3 coeff,
                        boost::any output
                        )
{
  if (vData.size() == 0) printf(" vData panic.\n");
  output = triLerp(vData[0], vData[1], vData[2], vData[3],
                   vData[4], vData[5], vData[6], vData[7],
                   coeff.getData());
}
#endif

//////////////////////////////////////////////////////////////////////////
// get tetra volume
// input
// cellId:	which cell
// return the volume of this cell
//////////////////////////////////////////////////////////////////////////
double RegularCartesianGrid::cellVolume(int cellId)
{
  double volume = (double)oneOvermappingFactorX*oneOvermappingFactorY*oneOvermappingFactorZ;
  return volume;
}

void RegularCartesianGrid::boundaryIntersection(VECTOR3& intersectP,
                        VECTOR3& startP,
                        VECTOR3& endP,
						float* stepSize, 
						float oldStepSize)
{
    VECTOR3 hitPoint;
	float thit;
	float xMin, xMax, yMin, yMax, zMin, zMax;

    intersectP.set(-1, -1, -1);

	xMax = (float)(m_nDimension[0] - 1);
	yMax = (float)(m_nDimension[1] - 1);
	zMax = (float)(m_nDimension[2] - 1);
	xMin = yMin = zMin = 0;

    VECTOR3 planeN[6];
    planeN[0].set(1, 0, 0);			//right
    planeN[1].set(-1, 0, 0);		//left
    planeN[2].set(0, 1, 0);			//top
    planeN[3].set(0,-1, 0);			//bottom
    planeN[4].set(0, 0, 1);			//front
    planeN[5].set(0, 0,-1);			//back

    VECTOR3 rayNormal;
    rayNormal = endP - startP;
    rayNormal.normalize();

    /*float vertex[8][3] = {	{xMin, yMin, zMin},		// rear, left, bottom
				{xMin, yMax, zMin},		// rear, left, top
				{xMax, yMin, zMin},		// rear, right, bottom
				{xMax, yMax, zMin},		// rear, right, top
				{xMin, yMin, zMax},		// front, left, bottom
				{xMin, yMax, zMax},		// front, left, top
				{xMax, yMin, zMax},		// front, right, bottom
				{xMax, yMax, zMax}};	// front, right, toop
                */
	
    VECTOR3 A;
	A = startP;

	// a point on the plane
    VECTOR3 B[6];
    B[0].set(xMax, (yMin + yMax)/2, (zMin + zMax)/2);			// right
    B[1].set(xMin, (yMin + yMax)/2, (zMin + zMax)/2);			// left
    B[2].set((xMax + xMin)/2, yMax, (zMin + zMax)/2);			// top
    B[3].set((xMax + xMin)/2, yMin, (zMin + zMax)/2);			// bottom
    B[4].set((xMax + xMin)/2, (yMin + yMax)/2, zMax);			// front
    B[5].set((xMax + xMin)/2, (yMin + yMax)/2, zMin);			// rear

    VECTOR3 bMinusA[6];
    bMinusA[0] = B[0] - A;
    bMinusA[1] = B[1] - A;
    bMinusA[2] = B[2] - A;
    bMinusA[3] = B[3] - A;
    bMinusA[4] = B[4] - A;
    bMinusA[5] = B[5] - A;

	float nDotBMinusA[6];
	nDotBMinusA[0] = dot(planeN[0], bMinusA[0]);
	nDotBMinusA[1] = dot(planeN[1], bMinusA[1]);
	nDotBMinusA[2] = dot(planeN[2], bMinusA[2]);
	nDotBMinusA[3] = dot(planeN[3], bMinusA[3]);
	nDotBMinusA[4] = dot(planeN[4], bMinusA[4]);
	nDotBMinusA[5] = dot(planeN[5], bMinusA[5]);

	float nDotc[6];
	nDotc[0] = dot(rayNormal, planeN[0]);
	nDotc[1] = dot(rayNormal, planeN[1]);
	nDotc[2] = dot(rayNormal, planeN[2]);
	nDotc[3] = dot(rayNormal, planeN[3]);
	nDotc[4] = dot(rayNormal, planeN[4]);
	nDotc[5] = dot(rayNormal, planeN[5]);

	// intersect with the right plan?
	if(nDotc[0] != 0)							// ray does not parallel with the right plane
	{
		// get hit point
		thit = nDotBMinusA[0]/nDotc[0];
		if(thit >= 0)
		{
            hitPoint.set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if((hitPoint[0] >= (xMax-EPS)) && (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS))&& (hitPoint[1] <= (yMax+EPS))		&&
				(hitPoint[2] >= (zMin-EPS))&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[0] = xMax;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the left plan?
	if(nDotc[1] != 0)							// ray does not parallel with the left plane
	{
		// get hit point
		thit = nDotBMinusA[1]/nDotc[1];
		if(thit >= 0)
		{
            hitPoint.set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if((hitPoint[0] >= (xMin-EPS)) && (hitPoint[0] <= (xMin+EPS)) &&
				(hitPoint[1] >= (yMin-EPS))&& (hitPoint[1] <= (yMax+EPS))		&&
				(hitPoint[2] >= (zMin-EPS))&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[0] = xMin;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the top plan?
	if(nDotc[2] != 0)							// ray does not parallel with the top plane
	{
		// get hit point
		thit = nDotBMinusA[2]/nDotc[2];
		if(thit >= 0)
		{
            hitPoint.set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if((hitPoint[0] >= (xMin-EPS))	&& (hitPoint[0] <= (xMax+EPS))		&&
				(hitPoint[1] >= (yMax-EPS))	&& (hitPoint[1] <= (yMax+EPS)) &&
				(hitPoint[2] >= (zMin-EPS))	&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[1] = yMax;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the bottom plan?
	if(nDotc[3] != 0)							// ray does not parallel with the bottom plane
	{
		// get hit point
		thit = nDotBMinusA[3]/nDotc[3];
		if(thit >= 0)
		{
            hitPoint.set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if( (hitPoint[0] >= (xMin-EPS))&& (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS))&& (hitPoint[1] <= (yMin+EPS)) &&
				(hitPoint[2] >= (zMin-EPS))&& (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[1] = yMin;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the front plan?
	if(nDotc[4] != 0)							// ray does not parallel with the front plane
	{
		// get hit point
		thit = nDotBMinusA[4]/nDotc[4];
		if(thit >= 0)
		{
            hitPoint.set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if( (hitPoint[0] >= (xMin-EPS)) && (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS)) && (hitPoint[1] <= (yMax+EPS)) &&
				(hitPoint[2] >= (zMax-EPS)) && (hitPoint[2] <= (zMax+EPS)))
			{
				// in plane
				hitPoint[2] = zMax;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}

	// intersect with the rear plan?
	if(nDotc[5] != 0)							// ray does not parallel with the rear plane
	{
		// get hit point
		thit = nDotBMinusA[5]/nDotc[5];
		if(thit >= 0)
		{
            hitPoint.set(A[0] + rayNormal[0]*thit, A[1] + rayNormal[1]*thit, A[2] + rayNormal[2]*thit);

			// check whether this point is in the plane range
			if( (hitPoint[0] >= (xMin-EPS)) && (hitPoint[0] <= (xMax+EPS)) &&
				(hitPoint[1] >= (yMin-EPS)) && (hitPoint[1] <= (yMax+EPS)) &&
				(hitPoint[2] >= (zMin-EPS)) && (hitPoint[2] <= (zMin+EPS)))
			{
				// in plane
				hitPoint[2] = zMin;
				intersectP = hitPoint;
				*stepSize = getStepSize(hitPoint, startP, endP, oldStepSize);
				return;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//																		//
//	definition of Irregular CartesianGrid Class							//
//																		//
//////////////////////////////////////////////////////////////////////////

/*
// constructor and deconstructor
IrregularCartesianGrid::IrregularCartesianGrid():CartesianGrid()
{
    Reset();
}

IrregularCartesianGrid::IrregularCartesianGrid(int xdim, int ydim, int zdim):CartesianGrid(xdim, ydim, zdim)
{
    Reset();
}

IrregularCartesianGrid::~IrregularCartesianGrid()
{
	delete[] m_pXSpacing;
	delete[] m_pYSpacing;
	delete[] m_pZSpacing;
    Reset();
}

void IrregularCartesianGrid::Reset(void)
{
	m_pXSpacing = m_pYSpacing = m_pZSpacing = NULL;
}
*/

#if 0
//////////////////////////////////////////////////////////////////////////
//																		//
//	definition of Irregular Grid Class									//
//																		//
//////////////////////////////////////////////////////////////////////////
// constructor and deconstructor
IrregularGrid::IrregularGrid()
{
    Reset();
}

IrregularGrid::IrregularGrid(int nodeNum, int tetraNum, 
							 CVertex* pVertexGeom, 
							 CTetra* pTetra, 
							 TVertex* pVertexTopo)
{
    Reset();
	m_nNodeNum = nodeNum;
	m_nTetraNum = tetraNum;
	m_pVertexGeom = pVertexGeom;
	m_pTetra = pTetra;
	m_pVertexTopo = pVertexTopo;
}

IrregularGrid::~IrregularGrid()
{
	if(m_pVertexGeom != NULL)
		delete[] m_pVertexGeom;

	if(m_pTetra != NULL)
		delete[] m_pTetra;

	if(m_pTetraInfo != NULL)
		delete[] m_pTetraInfo;

	if(m_pVertexTopo != NULL)
		delete[] m_pVertexTopo;
}

void IrregularGrid::Reset(void)
{
	m_nNodeNum = m_nTetraNum = 0;
	m_bTetraInfoInit = false;
	m_pVertexGeom = NULL;
	m_pTetra = NULL;
	m_pTetraInfo = NULL;
	m_pVertexTopo = NULL;
}

void IrregularGrid::setTetraInfoInit(bool bInit)
{
	m_bTetraInfoInit = bInit;
}

bool IrregularGrid::GetTetraInfoInit(void)
{
	return m_bTetraInfoInit;
}

//////////////////////////////////////////////////////////////////////////
// set bounding box
//////////////////////////////////////////////////////////////////////////
void IrregularGrid::setBoundary(VECTOR3& minB, VECTOR3& maxB)
{
	m_vMinBound = minB;
	m_vMaxBound = maxB;
}

void IrregularGrid::Boundary(VECTOR3& minB, VECTOR3& maxB)
{
	minB = m_vMinBound;
	maxB = m_vMaxBound;
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point pos is in bounding box
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::isInBBox(VECTOR3& pos)
{
	if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
		(pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
		(pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]))
		return true;
	else
		return false;
}

bool IrregularGrid::isInRealBBox(VECTOR3& pos)
{
  if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
      (pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
      (pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]))
    return true;
  else {
    return false;
  }
}

bool IrregularGrid::isInRealBBox(VECTOR3& pos, float t)
{
  if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
      (pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
      (pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]) &&
	  (     t >= m_vMinBound[3]) && (     t <= m_vMaxBound[3]) )
    return true;
  else {
    return false;
  }
}

void IrregularGrid::ComputeBBox(void)
{
	VECTOR3 minB, maxB;

    minB.set(FLT_MAX, FLT_MAX, FLT_MAX);
    maxB.set(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	for(int iFor = 0; iFor < m_nNodeNum; iFor++)
	{
		// x
		if(m_pVertexGeom[iFor].position[0] < minB[0])
			minB[0] = m_pVertexGeom[iFor].position[0];
		if(m_pVertexGeom[iFor].position[0] > maxB[0])
			maxB[0] = m_pVertexGeom[iFor].position[0];

		// y
		if(m_pVertexGeom[iFor].position[1] < minB[1])
			minB[1] = m_pVertexGeom[iFor].position[1];
		if(m_pVertexGeom[iFor].position[1] > maxB[1])
			maxB[1] = m_pVertexGeom[iFor].position[1];

		// z
		if(m_pVertexGeom[iFor].position[2] < minB[2])
			minB[2] = m_pVertexGeom[iFor].position[2];
		if(m_pVertexGeom[iFor].position[2] > maxB[2])
			maxB[2] = m_pVertexGeom[iFor].position[2];
	}

    setBoundary(minB, maxB);
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point is in the boundary
//////////////////////////////////////////////////////////////////////////
#if 0
bool IrregularGrid::at_phys(VECTOR3& pos)
{
	PointInfo pInfo;

	// whether in the bounding box
	if(!isInBBox(pos))
		return false;

    pInfo.set(pos, pInfo.interpolant, -1, -1);
	if(phys_to_cell(pInfo) == 1)
		return true;

	return false;
}
#endif

//////////////////////////////////////////////////////////////////////////
// get vertex list of a cell
//////////////////////////////////////////////////////////////////////////
ReturnStatus IrregularGrid::getCellVertices(int cellId,
								   vector<int>& vVertices)
{
	if((cellId < 0) || (cellId >= m_nTetraNum))
        return FAIL;

	vVertices.clear();
	for(int iFor = 0; iFor < 4; iFor++)
	{
		vVertices.push_back(m_pTetra[cellId].ver[iFor]);
	}

    return SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// from physical coordinate to natural coordinate
//
// input:
// pCoord:	physical position
// output:
// nCoord:	interpolation coefficients
// return:	return false if determinant is 0
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::Physical2NaturalCoord(VECTOR3& nCoord, VECTOR3& pCoord, int cellId)
{
	VECTOR3 B, v, AmulB;
	int idx, i;
	float determinant, oneOverDeter;

	try
	{
		i = 0;
		idx = m_pTetra[cellId].ver[i];
        v.set(m_pVertexGeom[idx].position[0], m_pVertexGeom[idx].position[1], m_pVertexGeom[idx].position[2]);

        B.set(pCoord[0] - v[0], pCoord[1] - v[1], pCoord[2] - v[2]);

		if(GetTetraInfoInit() == false)
		{
            setTetraInfoInit(true);
			if(m_pTetraInfo == NULL)
				m_pTetraInfo = new TetraInfo[m_nTetraNum];

			for(i = 0; i < m_nTetraNum; i++)
			{
				PreGetP2NMatrix(m_pTetraInfo[i].p2nMatrix, i);
				m_pTetraInfo[i].volume = cellVolume(i);
			}
		}

		AmulB = m_pTetraInfo[cellId].p2nMatrix * B;

		determinant = m_pTetraInfo[cellId].volume * (float)6.0;
		if(determinant == 0) return false;

		oneOverDeter = (float)1.0 / determinant;
		nCoord[0] = AmulB[0] * oneOverDeter;
		nCoord[1] = AmulB[1] * oneOverDeter;
		nCoord[2] = AmulB[2] * oneOverDeter;
		return true;
	}
	catch(...)
	{
		printf("Exceptions happen in function IrregularGrid::Physical2NaturalCoord()!\n");
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////
// get the physical coordinate of the vertex
//
// input:
// verIdx: index of vertex
// output:
// pos: physical coordinate of vertex
//////////////////////////////////////////////////////////////////////////
ReturnStatus IrregularGrid::at_vertex(int verIdx, VECTOR3& pos)
{
	if((verIdx < 0) || (verIdx >= m_nNodeNum))
        return FAIL;

    pos.set(m_pVertexGeom[verIdx].position[0],
			m_pVertexGeom[verIdx].position[1],
			m_pVertexGeom[verIdx].position[2]);
    return SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// whether the point in the physical position is in the cell
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// output:
// pInfo.interpolant: natural coordinate
// return:		returns 1 if in cell
//////////////////////////////////////////////////////////////////////////
bool IrregularGrid::isInCell(PointInfo& pInfo, const int cellId)
{
	VECTOR3 nCoord;

	if(!Physical2NaturalCoord(nCoord, pInfo.phyCoord, cellId))			// determinant is 0
	{
		return false;
	}

    pInfo.interpolant.set(nCoord[0], nCoord[1], nCoord[2]);
	// four conditions to be in a tetra
	// all three natural coordinates >= 0, and the sum of them <= 1
	if( (nCoord[0] >= 0.0) && (nCoord[1] >= 0.0) && (nCoord[2] >= 0.0) && ((1.0 - nCoord[0] - nCoord[1] - nCoord[2]) >= 0))	
		return true;
	else															// not in tetra
		return false;
}

//////////////////////////////////////////////////////////////////////////
//	find the next tetra, from which the point goes out
//
//	input: 
//	tetraId:	the tetra in check
//	pInfo:		the physical and natural coordinates of point
//
//	output:		the next tetra to search
//////////////////////////////////////////////////////////////////////////
//changed by lijie, find the tetra with the highest absolute value
int IrregularGrid::nextTetra(PointInfo& pInfo, int tetraId)
{
	
	int nextT;
	VECTOR3 nCoord;
/*
    nCoord.set(pInfo.interpolant[0], pInfo.interpolant[1], pInfo.interpolant[2]);

	if(nCoord[0] < 0.0)
		nextT = m_pTetra[tetraId].tetra[1];

	if(nCoord[1] < 0.0)
		nextT = m_pTetra[tetraId].tetra[2];

	if(nCoord[2] < 0.0)
		nextT = m_pTetra[tetraId].tetra[3];

	if((1.0 - nCoord[0] - nCoord[1] - nCoord[2]) < 0.0)
		nextT = m_pTetra[tetraId].tetra[0];

	//printf("nextT=%d\n",tetraId);
	return nextT;
	*/
	
	 nextT=-1;
//	VECTOR3 nCoord;

    nCoord.set(pInfo.interpolant[0], pInfo.interpolant[1], pInfo.interpolant[2]);

	float tmp_ncord[4];
	tmp_ncord[0]=1.0 - nCoord[0] - nCoord[1] - nCoord[2];
	tmp_ncord[1]=nCoord[0];
	tmp_ncord[2]=nCoord[1];
	tmp_ncord[3]=nCoord[2];

	bool found=false;
	int next_id;
	float maxi=0;
	for(int i=0; i<4; i++)
	{
		if(tmp_ncord[i]<0)
		{
			if(fabs(tmp_ncord[i])>=maxi)
			{
				maxi=fabs(tmp_ncord[i]);
				next_id=i;
				found=true;
				break;
			}
		}
	}
	
	if(found==true)
	nextT = m_pTetra[tetraId].tetra[next_id];
	
	//printf("ncord=%f %f %f %f\n",tmp_ncord[0],tmp_ncord[1],tmp_ncord[2],tmp_ncord[3]);
	//printf("nextT=%d\n",next_id);
	return nextT;

}

//////////////////////////////////////////////////////////////////////////
// get the cell id and also interpolating coefficients for the given
// physical position
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// fromCell:	if -1, initial point; otherwise this point is advected from others
// output:
// pInfo.inCell: in which cell
// pInfo.interpolant: natural coordinate
// return:	returns 1 if successful; otherwise returns -1
//////////////////////////////////////////////////////////////////////////
ReturnStatus IrregularGrid::phys_to_cell(PointInfo& pInfo)
{
	int iFor, tetraInCheck, nextT;
	bool bFind;
	bool* pNonCheckedCells;

    if(pInfo.fromCell == FAIL)		// no previous information
	{
		bFind = false;
		for(iFor = 0; iFor < m_nTetraNum; iFor++)
		{
			if(isInCell(pInfo, iFor))			// in this cell
			{
				pInfo.inCell = iFor;
				bFind = true;
                return SUCCESS;
			}
			if((bFind == false) && (iFor == m_nTetraNum))
                return FAIL;
		}
	}
	else							// generated from previous point
	{
		bFind = false;
		tetraInCheck = pInfo.fromCell;

		// in order to prevent circular check
		pNonCheckedCells = new bool[m_nTetraNum];
		for(iFor = 0; iFor < m_nTetraNum; iFor++)
			pNonCheckedCells[iFor] = false;

                while(!bFind && tetraInCheck>=0)
		{
			if(isInCell(pInfo, tetraInCheck))		// whether in this tetra
			{
				//pInfo.inCell = iFor;//wrong
				pInfo.inCell = 	tetraInCheck;		//changed by lijie
				delete pNonCheckedCells;
                return SUCCESS;
			}
			else									// not in tetraInCheck
			{
				// find which face to go out
				pNonCheckedCells[tetraInCheck] = true;
				nextT = nextTetra(pInfo, tetraInCheck);
				tetraInCheck = nextT;
                                if(tetraInCheck>=0 && pNonCheckedCells[tetraInCheck] == true)	// circular check
				{
					delete pNonCheckedCells;
                    return FAIL;
				}
			}
		}
	}

    return FAIL;
}

//////////////////////////////////////////////////////////////////////////
// barycentric interpolation
// input
// nodeData:	four corners of tetra
// coeff:		barycentric coordinate
// output
// vData:		output
//////////////////////////////////////////////////////////////////////////
void IrregularGrid::interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, VECTOR3 coeff)
{
	float fTemp[4];
    double fCoeff[3];

	fCoeff[0] = coeff[0];
	fCoeff[1] = coeff[1];
	fCoeff[2] = coeff[2];

	for(int iFor = 0; iFor < 3; iFor++)
	{
		fTemp[0] = vData[0][iFor];
		fTemp[1] = vData[1][iFor];
		fTemp[2] = vData[2][iFor];
		fTemp[3] = vData[3][iFor];

		nodeData[iFor] = BaryInterp(fTemp, fCoeff);
	}
}

//////////////////////////////////////////////////////////////////////////
// get tetra volume
// input
// cellId:	which cell
// return the volume of this cell
//////////////////////////////////////////////////////////////////////////
double IrregularGrid::cellVolume(int cellId)
{
	double volume;
    VECTOR3 v[4];
	int idx, i;
	
	for(i = 0; i < 4; i++)
	{
		idx = m_pTetra[cellId].ver[i];
        v[i].set(m_pVertexGeom[idx].position[0], m_pVertexGeom[idx].position[1], m_pVertexGeom[idx].position[2]);
	}

	volume = (v[1][0] - v[0][0])*((v[2][1] - v[0][1])*(v[3][2] - v[0][2]) - (v[2][2] - v[0][2])*(v[3][1] - v[0][1])) + 
			 (v[2][0] - v[0][0])*((v[0][1] - v[1][1])*(v[3][2] - v[0][2]) - (v[0][2] - v[1][2])*(v[3][1] - v[0][1])) +
			 (v[3][0] - v[0][0])*((v[1][1] - v[0][1])*(v[2][2] - v[0][2]) - (v[1][2] - v[0][2])*(v[2][1] - v[0][1]));

	volume /= 6.0;
	return volume;
}

//////////////////////////////////////////////////////////////////////////
// used to compute the natural coordinate of points
//////////////////////////////////////////////////////////////////////////


//new version, modified by lijie
void IrregularGrid::PreGetP2NMatrix(MATRIX3& m, int cellId)
{
	VECTOR3 v[5];
	int idx, i;
	
	//printf("cell %d vertices=\n",cellId);

	for(i = 1; i < 5; i++)
	{
		idx = m_pTetra[cellId].ver[i-1];
        v[i].set(m_pVertexGeom[idx].position[0], m_pVertexGeom[idx].position[1], m_pVertexGeom[idx].position[2]);
	//	printf("%f %f %f \n",v[i].x(),v[i].y(),v[i].z());
	}

	/*
	m[0][0] = (v[4].z()-v[1].z()) * (v[3].y() - v[4].y()) - (v[3].z()-v[4].z()) * (v[4].y()-v[1].y());
	m[0][1] = (v[4].z()-v[1].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[4].y()-v[1].y());
	m[0][2] = (v[2].z()-v[3].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[2].y()-v[3].y());

	m[1][0] = (v[4].x()-v[1].x()) * (v[3].z() - v[4].z()) - (v[3].x()-v[4].x()) * (v[4].z()-v[1].z());
	m[1][1] = (v[4].x()-v[1].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[4].z()-v[1].z());
	m[1][2] = (v[2].x()-v[3].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[2].z()-v[3].z());

	m[2][0] = (v[4].y()-v[1].y()) * (v[3].x() - v[4].x()) - (v[3].y()-v[4].y()) * (v[4].x()-v[1].x());
	m[2][1] = (v[4].y()-v[1].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[4].x()-v[1].x());
	m[2][2] = (v[2].y()-v[3].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[2].x()-v[3].x());
	*/

	m[0][0] = (v[4].z()-v[1].z()) * (v[3].y() - v[4].y()) - (v[3].z()-v[4].z()) * (v[4].y()-v[1].y());
	m[1][0] = (v[4].z()-v[1].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[4].y()-v[1].y());
	m[2][0] = (v[2].z()-v[3].z()) * (v[1].y() - v[2].y()) - (v[1].z()-v[2].z()) * (v[2].y()-v[3].y());

	m[0][1] = (v[4].x()-v[1].x()) * (v[3].z() - v[4].z()) - (v[3].x()-v[4].x()) * (v[4].z()-v[1].z());
	m[1][1] = (v[4].x()-v[1].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[4].z()-v[1].z());
	m[2][1] = (v[2].x()-v[3].x()) * (v[1].z() - v[2].z()) - (v[1].x()-v[2].x()) * (v[2].z()-v[3].z());

	m[0][2] = (v[4].y()-v[1].y()) * (v[3].x() - v[4].x()) - (v[3].y()-v[4].y()) * (v[4].x()-v[1].x());
	m[1][2] = (v[4].y()-v[1].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[4].x()-v[1].x());
	m[2][2] = (v[2].y()-v[3].y()) * (v[1].x() - v[2].x()) - (v[1].y()-v[2].y()) * (v[2].x()-v[3].x());

}

void IrregularGrid::GetDimension(int& xdim, int& ydim, int& zdim)
{
    xdim = m_nNodeNum; ydim = m_nTetraNum; zdim = 0;

}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////
// with the newly intersection point, get the interpolated step size for this new point
///////////////////////////////////////////////////////////////////////////////////////////////////////
float getStepSize(VECTOR3& p, VECTOR3& p1, VECTOR3& p2, float oldStepSize)
{
	float newStepSize;
	float p1p, p1p2;

	p1p = (float)sqrt((p[0] - p1[0]) * (p[0] - p1[0])+
					  (p[1] - p1[1]) * (p[1] - p1[1])+
					  (p[2] - p1[2]) * (p[2] - p1[2]));

	p1p2 = (float)sqrt((p2[0] - p1[0]) * (p2[0] - p1[0])+
					   (p2[1] - p1[1]) * (p2[1] - p1[1])+
					   (p2[2] - p1[2]) * (p2[2] - p1[2]));

	assert(p1p2 != 0);

	newStepSize = oldStepSize * (p1p / p1p2);

	return newStepSize;
}

} // namespace edda
