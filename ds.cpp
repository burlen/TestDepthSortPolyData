#include <iostream>
using namespace std;

#include <vtkCamera.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDepthSortPolyData.h>
#include "vtkDepthSortPolyData2.h"
#include <sys/time.h>

static double elapsed(timeval &tv0, timeval &tv1)
{ return tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec)/1e6; }

int getModeNumber(char code)
{
    switch (code)
    {
        case '1': return VTK_SORT_FIRST_POINT; break;
        case '2': return VTK_SORT_BOUNDS_CENTER; break;
        case '3': return VTK_SORT_PARAMETRIC_CENTER; break;
    }
    return -1;
}

const char *getModeString(char code)
{
    switch (code)
    {
        case '1': return "first point"; break;
        case '2': return "axes aligned bounds"; break;
        case '3': return "parametric center"; break;
    }
    return "invalid";
}

const char *getCodeString(char code)
{
    switch (code)
    {
        case '1': return "old"; break;
        case '2': return "new"; break;
    }
    return "invalid";
}

int main(int ac, char **av)
{
    if (ac != 4)
    {
        cerr << "ERROR: must provide 1) mode 2) input file 3) output file." << endl << endl
            << "where: mode is a two char string. The first char can be either \"1\" for the old class or "
            "\"2\" for the new class. the second char can be either \"1\" to sort by cells "
            "first point, \"2\" to sort by the cells axis aligned bounding box centroid, or "
            "\"3\" to sort by the cells paramteric centroid. If the output file is \"\" output "
            "is skipped." << endl;
        return -1;
    }

    cerr << "Reading ... ";
    vtkPolyDataReader *r = vtkPolyDataReader::New();
    r->SetFileName(av[2]);
    r->Update();
    cerr << "ok!" << endl;


    vtkCamera *cam = vtkCamera::New();
    cam->SetPosition(128, 128, 1024);
    cam->SetFocalPoint(128, 128, 128);

    vtkPolyDataAlgorithm *a = NULL;
    if (av[1][0] == '1')
    {
        vtkDepthSortPolyData *ds = vtkDepthSortPolyData::New();
        ds->SetCamera(cam);
        ds->SortScalarsOn();
        ds->SetDepthSortMode(getModeNumber(av[1][1]));
        a = static_cast<vtkDepthSortPolyData*>(ds);
    }
    else
    {
        vtkDepthSortPolyData2 *ds = vtkDepthSortPolyData2::New();
        ds->SetCamera(cam);
        ds->SortScalarsOn();
        ds->SetDepthSortMode(getModeNumber(av[1][1]));
        a = static_cast<vtkDepthSortPolyData2*>(ds);
    }

    a->SetInputConnection(r->GetOutputPort());
    r->Delete();

    cerr << "Depth sort ... ";

    struct timeval tv0,tv1;
    gettimeofday(&tv0, 0);

    a->Update();

    gettimeofday(&tv1, 0);
    cerr << "ok!" << endl
        << getCodeString(av[1][0]) << " code sorting by " << getModeString(av[1][1])
        << " = " << elapsed(tv0, tv1) << " seconds" << endl;

    if (strlen(av[3]))
    {
        cerr << "Writing ... ";
        vtkPolyDataWriter *w = vtkPolyDataWriter::New();
        w->SetFileName(av[3]);
        w->SetInputData(a->GetOutput());
        w->Write();
        w->Delete();
        cerr << "ok!" << endl;
    }

    a->Delete();
    cam->Delete();

    return 0;
}
