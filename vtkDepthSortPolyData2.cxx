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

namespace {
template <typename T>
bool greater(const T *l, const T *r)
{ return l[0] > r[0]; }

template <typename T>
bool less(const T *l, const T *r)
{ return l[0] < r[0]; }

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

template <typename T>
T getCellMin(vtkIdType *pids, vtkIdType nPids, const T *px)
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
T getCellMax(vtkIdType *pids, vtkIdType nPids, const T *px)
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
void getCellCenterDepth(vtkPolyData *pds, vtkDataArray *gpts,
    vtkIdType nCells, T x0, T y0, T z0, T vx, T vy, T vz,
    T *&d, T **&dd)
{
  if (nCells < 1)
    {
    return;
    }

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

  d = static_cast<T*>(malloc(sizeof(T)*nCells));
  dd = static_cast<T**>(malloc(sizeof(T*)*nCells));
  for (vtkIdType cid = 0; cid < nCells; ++cid)
    {
    // get the cell points using the fast api
    vtkIdType *pids = NULL;
    vtkIdType nPids = 0;
    pds->GetCellPoints(cid, nPids, pids);

    // compute the cell bounds
    T xmn = getCellMin(pids, nPids, px);
    T xmx = getCellMax(pids, nPids, px);

    T ymn = getCellMin(pids, nPids, py);
    T ymx = getCellMax(pids, nPids, py);

    T zmn = getCellMin(pids, nPids, pz);
    T zmx = getCellMax(pids, nPids, pz);

    // compute the distance to the cell center
    d[cid] = (((xmn + xmx)/T(2)) - x0)*vx
      + (((ymn + ymx)/T(2)) - y0)*vy + (((zmn + zmx)/T(2)) - z0)*vz;

    // set up the sort structure
    dd[cid] = d + cid;
    }
}

template <typename T>
void getCellPoint0Depth(vtkPolyData *pds, vtkDataArray *gpts,
    vtkIdType nCells, T x0, T y0, T z0, T vx, T vy, T vz,
    T *&d, T **&dd)
{
  if (nCells < 1)
    {
    return;
    }

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

  d = static_cast<T*>(malloc(sizeof(T)*nCells));
  dd = static_cast<T**>(malloc(sizeof(T*)*nCells));
  for (vtkIdType cid = 0; cid < nCells; ++cid)
    {
    // get the cell points using the fast api
    vtkIdType *pids = NULL;
    vtkIdType nPids = 0;
    pds->GetCellPoints(cid, nPids, pids);

    // compute the distance to the cell's first point
    vtkIdType ii = pids[0];
    d[cid] = (px[3*ii] - x0)*vx + (py[3*ii] - y0)*vy + (pz[3*ii] - z0)*vz;

    // set up the sort structure
    dd[cid] = d + cid;
    }
}
};

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


  // Compute the sort vector
  double vector[3] = {0.0};
  double origin[3] = {0.0};
  if (this->Direction == VTK_DIRECTION_SPECIFIED_VECTOR)
    {
    memcpy(vector, this->Vector, 3*sizeof(double));
    memcpy(origin, this->Origin, 3*sizeof(double));
    }
  else //compute view vector
    {
    if (!this->Camera)
      {
      vtkErrorMacro("Need a camera to sort");
      return 0;
      }
    this->ComputeProjectionVector(vector, origin);
    }

  // create temporary input
  vtkPolyData *tmpInput = vtkPolyData::New();
  tmpInput->CopyStructure(input);

  vtkIdType nCells=input->GetNumberOfCells();

  vtkIdType *orderedIds
    = static_cast<vtkIdType*>(malloc(nCells*sizeof(vtkIdType)));

  if (nCells)
    {
    if ((this->DepthSortMode == VTK_SORT_FIRST_POINT)
      || (this->DepthSortMode == VTK_SORT_BOUNDS_CENTER))
      {
      vtkDataArray *pts = tmpInput->GetPoints()->GetData();
      switch (pts->GetDataType())
        {
        vtkTemplateMacro(

          // work in native precision
          VTK_TT x0 = static_cast<VTK_TT>(origin[0]);
          VTK_TT y0 = static_cast<VTK_TT>(origin[1]);
          VTK_TT z0 = static_cast<VTK_TT>(origin[2]);
          VTK_TT vx = static_cast<VTK_TT>(vector[0]);
          VTK_TT vy = static_cast<VTK_TT>(vector[1]);
          VTK_TT vz = static_cast<VTK_TT>(vector[2]);

          // compute the cell's depth
          VTK_TT *d = NULL; // depths
          VTK_TT **dd = NULL; // sort structure
          if (this->DepthSortMode == VTK_SORT_FIRST_POINT)
            {
            ::getCellPoint0Depth(
              tmpInput, pts, nCells, x0, y0, z0, vx, vy, vz, d, dd);
            }
          else
            {
            ::getCellCenterDepth(
              tmpInput, pts, nCells, x0, y0, z0, vx, vy, vz, d, dd);
            }

          // sort
          if (this->Direction == VTK_DIRECTION_FRONT_TO_BACK)
            {
            std::sort(dd, dd + nCells, ::less<VTK_TT>);
            }
          else
            {
            std::sort(dd, dd + nCells, ::greater<VTK_TT>);
            }

          // get the cell id order
          for (vtkIdType i = 0; i < nCells; ++i)
            {
            orderedIds[i] = dd[i] - d;
            }

          free(d);
          free(dd);
          );
        }
      }
    else // VTK_SORT_PARAMETRIC_CENTER
      {
      vtkGenericCell *cell = vtkGenericCell::New();

      double x[3] = {0.0};
      double p[3] = {0.0};

      size_t maxCellSize = input->GetMaxCellSize();
      double *w = static_cast<double*>(malloc(maxCellSize*sizeof(double)));
      double *d = static_cast<double*>(malloc(nCells*sizeof(double)));
      double **dd = static_cast<double**>(malloc(nCells*sizeof(double*)));
      for (vtkIdType cid = 0; cid < nCells; ++cid)
        {
        tmpInput->GetCell(cid, cell);
        int subId = cell->GetParametricCenter(p);
        cell->EvaluateLocation(subId, p, x, w);

        // compute the distance
        d[cid] = (x[0] - origin[0])*vector[0]
            + (x[1] - origin[1])*vector[1] + (x[2] - origin[2])*vector[2];

        // set up the sort structure
        dd[cid] = d + cid;
        }

      // sort
      if (this->Direction == VTK_DIRECTION_FRONT_TO_BACK)
        {
        std::sort(dd, dd + nCells, ::less<double>);
        }
      else
        {
        std::sort(dd, dd + nCells, ::greater<double>);
        }

      // get the cell id order
      for (vtkIdType i = 0; i < nCells; ++i)
        {
        orderedIds[i] = dd[i] - d;
        }

      free(w);
      free(d);
      free(dd);
      cell->Delete();
      }
    }

  // construct the output
  vtkCellData *inCD = input->GetCellData();
  vtkCellData *outCD = output->GetCellData();
  outCD->CopyAllocate(inCD);
  output->Allocate(tmpInput, nCells);
  output->SetPoints(input->GetPoints());
  output->GetPointData()->PassData(input->GetPointData());

  for (vtkIdType i = 0; i < nCells; ++i)
    {
    // get the cell points using the "fast" api
    vtkIdType *pids = NULL;
    vtkIdType nPids = 0;
    vtkIdType cid = orderedIds[i];
    tmpInput->GetCellPoints(cid, nPids, pids);
    int ctype = tmpInput->GetCellType(cid);

    // add to the ouput
    vtkIdType newId = output->InsertNextCell(ctype, nPids, pids);

    // copy over data
    outCD->CopyData(inCD, cid, newId);
    }

  if (this->SortScalars)
    {
    // add the sort indices
    vtkIdTypeArray *sortScalars = vtkIdTypeArray::New();
    sortScalars->SetArray(orderedIds, nCells, 0, 0);
    output->GetCellData()->SetScalars(sortScalars);
    sortScalars->Delete();
    }
  else
    {
    free(orderedIds);
    }

  tmpInput->Delete();

  return 1;
}

void vtkDepthSortPolyData2::ComputeProjectionVector(double vector[3],
                                                   double origin[3])
{
  double *focalPoint = this->Camera->GetFocalPoint();
  double *position = this->Camera->GetPosition();

  // If a camera is present, use it
  if (!this->Prop3D)
    {
    memcpy(origin, position, 3*sizeof(double));
    for(int i = 0; i < 3; ++i)
      {
      vector[i] = focalPoint[i] - position[i];
      }
    }
  else  //Otherwise, use Prop3D
    {
    this->Transform->SetMatrix(this->Prop3D->GetMatrix());
    this->Transform->Push();
    this->Transform->Inverse();

    double focalPt[4];
    memcpy(focalPt, focalPoint, 3*sizeof(double));
    focalPt[3] = 1.0;
    this->Transform->TransformPoint(focalPt, focalPt);

    double pos[4];
    memcpy(pos, position, 3*sizeof(double));
    pos[3] = 1.0;
    this->Transform->TransformPoint(pos, pos);

    memcpy(origin, pos, 3*sizeof(double));

    for (int i = 0; i < 3; ++i)
      {
      vector[i] = focalPt[i] - pos[i];
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
