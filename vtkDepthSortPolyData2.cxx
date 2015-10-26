/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDepthSortPolyData2.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkDepthSortPolyData2.h"

#include "vtkCamera.h"
#include "vtkCellData.h"
#include "vtkGenericCell.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkProp3D.h"
#include "vtkTransform.h"
#include "vtkUnsignedIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCharArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkShortArray.h"
#include "vtkSignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkIdTypeArray.h"
#include "vtkDataArray.h"

#include <algorithm>
#include <limits>
#include <cstdlib>

vtkStandardNewMacro(vtkDepthSortPolyData2);

vtkCxxSetObjectMacro(vtkDepthSortPolyData2,Camera,vtkCamera);

vtkDepthSortPolyData2::vtkDepthSortPolyData2()
{
  this->Camera = NULL;
  this->Prop3D = NULL;
  this->Direction = VTK_DIRECTION_BACK_TO_FRONT;
  this->DepthSortMode = VTK_SORT_FIRST_POINT;
  this->Vector[0] = this->Vector[1] = 0.0;
  this->Vector[2] = 0.0;
  this->Origin[0] = this->Origin[1] = this->Origin[2] = 0.0;
  this->Transform = vtkTransform::New();
  this->SortScalars = 0;
}

vtkDepthSortPolyData2::~vtkDepthSortPolyData2()
{
  this->Transform->Delete();

  if ( this->Camera )
    {
    this->Camera->Delete();
    }

  //Note: vtkProp3D is not deleted to avoid reference count cycle
}

// Don't reference count to avoid nasty cycle
void vtkDepthSortPolyData2::SetProp3D(vtkProp3D *prop3d)
{
  if ( this->Prop3D != prop3d )
    {
    this->Prop3D = prop3d;
    this->Modified();
    }
}

vtkProp3D *vtkDepthSortPolyData2::GetProp3D()
{
  return this->Prop3D;
}

struct CellDepth
{
  CellDepth() : z(), id(0) {}
  CellDepth(double az, double aid) : z(az), id(aid) {}
  CellDepth(const CellDepth &o) : z(o.z), id(o.id) {}
  void operator=(const CellDepth &o) { z = o.z; id = o.id; }
  double z;
  vtkIdType id;
};

static inline
bool greater(const CellDepth &l, const CellDepth &r)
{ return l.z > r.z; }

static inline
bool less(const CellDepth &l, const CellDepth &r)
{ return l.z < r.z; }

template <typename T> struct ptsVTKTt {};
template <> struct ptsVTKTt<float> { typedef vtkFloatArray vtkType; };
template <> struct ptsVTKTt<double> { typedef vtkDoubleArray vtkType; };
template <> struct ptsVTKTt<int> { typedef vtkIntArray vtkType; };
template <> struct ptsVTKTt<unsigned int> { typedef vtkUnsignedIntArray vtkType; };
template <> struct ptsVTKTt<char> { typedef vtkCharArray vtkType; };
template <> struct ptsVTKTt<signed char> { typedef vtkSignedCharArray vtkType; };
template <> struct ptsVTKTt<unsigned char> { typedef vtkUnsignedCharArray vtkType; };
template <> struct ptsVTKTt<long> { typedef vtkLongArray vtkType; };
template <> struct ptsVTKTt<unsigned long> { typedef vtkUnsignedLongArray vtkType; };
template <> struct ptsVTKTt<long long> { typedef vtkLongLongArray vtkType; };
template <> struct ptsVTKTt<unsigned long long> { typedef vtkUnsignedLongLongArray vtkType; };
template <> struct ptsVTKTt<short> { typedef vtkShortArray vtkType; };
template <> struct ptsVTKTt<unsigned short> { typedef vtkUnsignedShortArray vtkType; };
//template <> struct ptsVTKTt<vtkIdType> { typedef vtkIdTypeArray vtkType; };

template <typename T>
T getMin(vtkIdType *pids, vtkIdType nPids, const T *px)
{
  T mn = std::numeric_limits<T>::max();
  for (vtkIdType i = 0; i < nPids; ++i)
    {
    vtkIdType ii = pids[i];
    mn = px[3*ii] < mn ? px[3*ii] : mn;
    }
  return mn;
}

template <typename T>
T getMax(vtkIdType *pids, vtkIdType nPids, const T *px)
{
  T mx = -std::numeric_limits<T>::max();
  for (vtkIdType i = 0; i < nPids; ++i)
    {
    vtkIdType ii = pids[i];
    mx = px[3*ii] > mx ? px[3*ii] : mx;
    }
  return mx;
}

template <typename T>
void getCellCenters(vtkPolyData *pds, vtkDataArray *gpts, vtkIdType nCells,
    T * __restrict__ cx, T * __restrict__ cy, T * __restrict__ cz)
{
  if (nCells < 1)
    {
    return;
    }

  // ok, the whole point of this is to work with points in
  // their native precision. here we go.
  typename ptsVTKTt<T>::vtkType *pts
    = static_cast<typename ptsVTKTt<T>::vtkType*>(gpts);

  T *ppts = pts->GetPointer(0);
  T *px = ppts;
  T *py = ppts + 1;
  T *pz = ppts + 2;

  // this call insures that BuildCells gets done if it's
  // needed and we can use the faster GetCellPoints api
  // that doesn't check
  pds->GetCellType(0);

  vtkIdType *pids = NULL;
  vtkIdType nPids = 0;
  for (vtkIdType cid = 0; cid < nCells; ++cid)
    {
    // get the cell points using the "fast" api
    pds->GetCellPoints(cid, nPids, pids);

    // compute the cell centroid.
    T xmn = getMin(pids, nPids, px);
    T xmx = getMax(pids, nPids, px);
    cx[cid] = (xmn + xmx)/T(2);

    T ymn = getMin(pids, nPids, py);
    T ymx = getMax(pids, nPids, py);
    cy[cid] = (ymn + ymx)/T(2);

    T zmn = getMin(pids, nPids, pz);
    T zmx = getMax(pids, nPids, pz);
    cz[cid] = (zmn + zmx)/T(2);
    }
}

template <typename T>
void getCellCenterDepth(vtkPolyData *pds, vtkDataArray *gpts, vtkIdType nCells,
    T x0, T y0, T z0, T vx, T vy, T vz, std::vector<CellDepth> &depth)
{
  if (nCells < 1)
    {
    return;
    }

  // ok, the whole point of this is to work with points in
  // their native precision. here we go.
  typename ptsVTKTt<T>::vtkType *pts
    = static_cast<typename ptsVTKTt<T>::vtkType*>(gpts);

  T *ppts = pts->GetPointer(0);

  // now get pointers to each of the x,y,z values
  T *px = ppts;
  T *py = ppts + 1;
  T *pz = ppts + 2;

  // this call insures that BuildCells gets done if it's
  // needed and we can use the faster GetCellPoints api
  // that doesn't check
  pds->GetCellType(0);

  vtkIdType *pids = NULL;
  vtkIdType nPids = 0;
  for (vtkIdType cid = 0; cid < nCells; ++cid)
    {
    // get the cell points using the "fast" api
    pds->GetCellPoints(cid, nPids, pids);

    // compute the cell bounds
    T xmn = getMin(pids, nPids, px);
    T xmx = getMax(pids, nPids, px);

    T ymn = getMin(pids, nPids, py);
    T ymx = getMax(pids, nPids, py);

    T zmn = getMin(pids, nPids, pz);
    T zmx = getMax(pids, nPids, pz);

    // distance to the centroid along v
    T d = (((xmn + xmx)/T(2)) - x0)*vx
      + (((ymn + ymx)/T(2)) - y0)*vy + (((zmn + zmx)/T(2)) - z0)*vz;

    depth[cid].z = static_cast<double>(d);
    depth[cid].id = cid;
    }
}

template <typename T>
void getCellPoint0(vtkPolyData *pds, vtkDataArray *gpts, vtkIdType nCells,
    T * __restrict__ cx, T * __restrict__ cy, T * __restrict__ cz)
{
  if (nCells < 1)
    {
    return;
    }

  // ok, the whole point of this is to work with points in
  // their native precision. here we go.
  typename ptsVTKTt<T>::vtkType *pts
    = static_cast<typename ptsVTKTt<T>::vtkType*>(gpts);

  T *ppts = pts->GetPointer(0);

  // now get pointers to x,y,z values
  T *px = ppts;
  T *py = ppts + 1;
  T *pz = ppts + 2;

  // the following seemingly pointless call is because there is no
  // vtkPolyData::NeedToBuildCells method. this call insures that
  // BuildCells gets done if it's needed and we can use the faster
  // GetCellPoints api on our performance critical loop below.
  pds->GetCellType(0);

  vtkIdType *pids = NULL;
  vtkIdType nPids = 0;
  for (vtkIdType cid = 0; cid < nCells; ++cid)
    {
    // get the cell points using the "fast" api
    pds->GetCellPoints(cid, nPids, pids);

    vtkIdType ii = pids[0];
    cx[cid] = px[3*ii];
    cy[cid] = py[3*ii];
    cz[cid] = pz[3*ii];
    }
}

int vtkDepthSortPolyData2::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType numCells=input->GetNumberOfCells();

  // Compute the sort vector
  double vector[3] = {0.0};
  double origin[3] = {0.0};
  if ( this->Direction == VTK_DIRECTION_SPECIFIED_VECTOR )
    {
    for (int i=0; i<3; i++)
      {
      vector[i] = this->Vector[i];
      origin[i] = this->Origin[i];
      }
    }
  else //compute view vector
    {
    if (this->Camera == NULL)
      {
      vtkErrorMacro(<<"Need a camera to sort");
      return 0;
      }
    this->ComputeProjectionVector(vector, origin);
    }

  // Create temporary input, almost anything we
  // do would modify it.
  vtkPolyData *tmpInput = vtkPolyData::New();
  tmpInput->CopyStructure(input);
  vtkGenericCell *cell = vtkGenericCell::New();

  // Compute the depth value
  std::vector<CellDepth> depth(numCells);
  if (numCells)
    {
    if ((this->DepthSortMode == VTK_SORT_FIRST_POINT)
      || (this->DepthSortMode == VTK_SORT_BOUNDS_CENTER))
      {
      vtkDataArray *pts = tmpInput->GetPoints()->GetData();
      switch (pts->GetDataType())
        {
        vtkTemplateMacro(
          size_t nn = sizeof(VTK_TT)*numCells;
          // compute cell centers
          VTK_TT *cx = static_cast<VTK_TT*>(malloc(nn));
          VTK_TT *cy = static_cast<VTK_TT*>(malloc(nn));
          VTK_TT *cz = static_cast<VTK_TT*>(malloc(nn));
          if (this->DepthSortMode == VTK_SORT_FIRST_POINT)
            {
            getCellPoint0(tmpInput, pts, numCells, cx, cy, cz);
            }
         else
            {
            getCellCenters(tmpInput, pts, numCells, cx, cy, cz);
            }
          // compute distance to cell center
          VTK_TT x0 = static_cast<VTK_TT>(origin[0]);
          VTK_TT y0 = static_cast<VTK_TT>(origin[1]);
          VTK_TT z0 = static_cast<VTK_TT>(origin[2]);
          VTK_TT vx = static_cast<VTK_TT>(vector[0]);
          VTK_TT vy = static_cast<VTK_TT>(vector[1]);
          VTK_TT vz = static_cast<VTK_TT>(vector[2]);
          VTK_TT *d = static_cast<VTK_TT*>(malloc(nn));
          for (vtkIdType cid = 0; cid < numCells; ++cid)
            {
            d[cid] = (cx[cid] - x0)*vx + (cy[cid] - y0)*vy + (cz[cid] - z0)*vz;
            }
          free(cx);
          free(cy);
          free(cz);
          // finally store the result. this is done here in
          // a separate loop because the auto-vectorizer fails
          // when it's in the above loop.
          for (vtkIdType cid = 0; cid < numCells; ++cid)
            {
            CellDepth &cdep = depth[cid];
            cdep.z = static_cast<double>(d[cid]);
            cdep.id = cid;
            }
          free(d);
          );
        }
      }
    else // VTK_SORT_PARAMETRIC_CENTER
      {
      double *w = new double [input->GetMaxCellSize()];
      double x[3] = {0.0};
      double p[3] = {0.0};
      for (vtkIdType cid = 0; cid < numCells; ++cid)
        {
        tmpInput->GetCell(cid, cell);
        int subId = cell->GetParametricCenter(p);
        cell->EvaluateLocation(subId, p, x, w);
  
        double z = (x[0] - origin[0])*vector[0]
            + (x[1] - origin[1])*vector[1] + (x[2] - origin[2])*vector[2];
  
        CellDepth &d = depth[cid];
        d.z = z;
        d.id = cid;
        }
      delete [] w;
      }
  
    this->UpdateProgress(0.20);
  
    // Sort the depths
    if (this->Direction == VTK_DIRECTION_FRONT_TO_BACK)
      {
      std::sort(depth.begin(), depth.end(), less);
      }
    else
      {
      std::sort(depth.begin(), depth.end(), greater);
      }
    }

  this->UpdateProgress(0.60);

  // Generate sorted output
  vtkIdTypeArray *sortScalars = NULL;
  vtkIdType *scalars = NULL;
  if (this->SortScalars)
    {
    sortScalars = vtkIdTypeArray::New();
    sortScalars->SetNumberOfTuples(numCells);
    scalars = sortScalars->GetPointer(0);
    }

  vtkCellData *inCD = input->GetCellData();
  vtkCellData *outCD = output->GetCellData();
  outCD->CopyAllocate(inCD);
  output->Allocate(tmpInput, numCells);

  vtkIdType *pids = NULL;
  vtkIdType nPids = 0;
  for (vtkIdType i = 0; i < numCells; ++i)
    {
    // get the cell points using the "fast" api
    vtkIdType cid = depth[i].id;
    tmpInput->GetCellPoints(cid, nPids, pids);
    int ctype = tmpInput->GetCellType(cid);

    // add to the ouput
    vtkIdType newId = output->InsertNextCell(ctype, nPids, pids);

    // copy over data
    outCD->CopyData(inCD, cid, newId);

    if (this->SortScalars)
      {
      scalars[newId] = newId;
      }
    }

  this->UpdateProgress(0.90);

  // Points are left alone
  output->SetPoints(input->GetPoints());
  output->GetPointData()->PassData(input->GetPointData());

  // add the sort indices
  if (this->SortScalars)
    {
    int idx = output->GetCellData()->AddArray(sortScalars);
    output->GetCellData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    sortScalars->Delete();
    }

  // Clean up and get out
  tmpInput->Delete();
  cell->Delete();

  // we explicitly pre-allocated the number of cells we
  // need. assuming Allocate does something sane then
  // squeezing is probably uneccessary.
  //output->Squeeze();

  return 1;
}

void vtkDepthSortPolyData2::ComputeProjectionVector(double vector[3],
                                                   double origin[3])
{
  double *focalPoint = this->Camera->GetFocalPoint();
  double *position = this->Camera->GetPosition();

  // If a camera is present, use it
  if ( !this->Prop3D )
    {
    for(int i=0; i<3; i++)
      {
      vector[i] = focalPoint[i] - position[i];
      origin[i] = position[i];
      }
    }

  else  //Otherwise, use Prop3D
    {
    double focalPt[4], pos[4];
    int i;

    this->Transform->SetMatrix(this->Prop3D->GetMatrix());
    this->Transform->Push();
    this->Transform->Inverse();

    for(int i=0; i<4; ++i)
      {
      focalPt[i] = focalPoint[i];
      pos[i] = position[i];
      }

    this->Transform->TransformPoint(focalPt,focalPt);
    this->Transform->TransformPoint(pos,pos);

    for (i=0; i<3; i++)
      {
      vector[i] = focalPt[i] - pos[i];
      origin[i] = pos[i];
      }
    this->Transform->Pop();
  }
}

unsigned long int vtkDepthSortPolyData2::GetMTime()
{
  unsigned long mTime=this->Superclass::GetMTime();

  if ( this->Direction != VTK_DIRECTION_SPECIFIED_VECTOR )
    {
    unsigned long time;
    if ( this->Camera != NULL )
      {
      time = this->Camera->GetMTime();
      mTime = ( time > mTime ? time : mTime );
      }

    if ( this->Prop3D != NULL )
      {
      time = this->Prop3D->GetMTime();
      mTime = ( time > mTime ? time : mTime );
      }
    }

  return mTime;
}

void vtkDepthSortPolyData2::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->Camera )
    {
    os << indent << "Camera:\n";
    this->Camera->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Camera: (none)\n";
    }

  if ( this->Prop3D )
    {
    os << indent << "Prop3D:\n";
    this->Prop3D->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Prop3D: (none)\n";
    }

  os << indent << "Direction: ";
  if ( this->Direction == VTK_DIRECTION_BACK_TO_FRONT )
    {
    os << "Back To Front" << endl;
    }
  else if ( this->Direction == VTK_DIRECTION_FRONT_TO_BACK )
    {
    os << "Front To Back";
    }
  else
    {
    os << "Specified Direction: ";
    os << "(" << this->Vector[0] << ", " << this->Vector[1] << ", "
       << this->Vector[2] << ")\n";
    os << "Specified Origin: ";
    os << "(" << this->Origin[0] << ", " << this->Origin[1] << ", "
       << this->Origin[2] << ")\n";
    }

  os << indent << "Depth Sort Mode: ";
  if ( this->DepthSortMode == VTK_SORT_FIRST_POINT )
    {
    os << "First Point" << endl;
    }
  else if ( this->DepthSortMode == VTK_SORT_BOUNDS_CENTER )
    {
    os << "Bounding Box Center" << endl;
    }
  else
    {
    os << "Paramteric Center" << endl;
    }

  os << indent << "Sort Scalars: " << (this->SortScalars ? "On\n" : "Off\n");
}
