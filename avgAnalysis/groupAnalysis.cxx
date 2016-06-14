/*
* DiagnosticIndex research code
* Beatriz Paniagua UNC
*
*/
#include "PCAModelBuilder.h"
#include "StatisticalModel.h"
#include "DataManager.h"
#include "vtkStandardMeshRepresenter.h"
#include <boost/scoped_ptr.hpp>

#include "vtkDirectory.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "CommonTypes.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include <vtkDelimitedTextReader.h> //CSVloader
#include <vtkTable.h> //CSVloader

using namespace statismo;


typedef vtkStandardMeshRepresenter RepresenterType;
typedef DataManager<vtkPolyData> DataManagerType;
typedef PCAModelBuilder<vtkPolyData> ModelBuilderType;
typedef StatisticalModel<vtkPolyData> StatisticalModelType;
typedef std::vector<std::string> StringVectorType;

int getdir (std::string dir, std::vector<std::string> &files, const std::string& extension=".*") {
vtkDirectory *directory = vtkDirectory::New();
directory->Open(dir.c_str());
for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
const char* filename = directory->GetFile(i);
if (extension == ".*" || std::string(filename).find(extension) != std::string::npos)
files.push_back(filename);
}
directory->Delete();
return 0;
}

vtkPolyData* ScaleVTK (vtkPolyData* reference)
{
    // Sum up Original Points
    double sum = 0;

    for (int i = 0; i < reference->GetNumberOfPoints(); i++)
    {
    double curPoint[3];
    reference->GetPoint(i, curPoint);
    for( unsigned int dim = 0; dim < 3; dim++ )
        {
          sum += curPoint[dim]*curPoint[dim];
        }

    }

    //Calculate Scale
    double scale = sqrt(sum) / (reference->GetNumberOfPoints()*3);

    // Create a New Point Set "scaledpoint" with the scaled Values
    vtkPoints * scaledpoints = vtkPoints::New();
    scaledpoints = reference->GetPoints();

    for( int pointID = 0; pointID < reference->GetNumberOfPoints(); pointID++ )
      {
        double curPoint[3];
        reference->GetPoint(pointID, curPoint);
        double scaledPoint[3];
      for( unsigned int dim = 0; dim < 3; dim++ )
        {
        scaledPoint[dim] = curPoint[dim] / scale;
        }
      scaledpoints->SetPoint(pointID, scaledPoint);
      }
    // Set the scaled Points back to original mesh
    vtkPolyData* output = vtkPolyData::New();
    output->DeepCopy(reference);
    output->SetPoints(scaledpoints);

    return output;
}

vtkPolyData* loadVTKPolyData(const std::string& filename)
{
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkPolyData* pd = vtkPolyData::New();
    pd->ShallowCopy(reader->GetOutput());

    vtkPolyData* pds = ScaleVTK(pd);

    return pds;
}

void saveSample(const vtkPolyData* pd, const std::string& resdir, const std::string& basename) {
    std::string filename = resdir +std::string("/") + basename;

    vtkPolyDataWriter* w = vtkPolyDataWriter::New();
#if (VTK_MAJOR_VERSION == 5 )
    w->SetInput(const_cast<vtkPolyData*>(pd));
#else
    w->SetInputData(const_cast<vtkPolyData*>(pd));
#endif
    w->SetFileName(filename.c_str());
    w->Update();
}

vtkPolyData* averageVTK (std::string datadir)
{
    StringVectorType filenames;
    getdir(datadir, filenames, ".vtk");
    if (filenames.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadir << " exiting.";
    exit(-1);
    }

    vtkPolyData *polydata0 = loadVTKPolyData(datadir+"/" + filenames[0]);
    vtkPoints *avgPoints = vtkPoints::New();
    avgPoints=polydata0->GetPoints();

    int numMeshes = filenames.size();
    for( int index = 1; index < numMeshes; index++ )
    {
      vtkPolyData* polydata = vtkPolyData::New();
      polydata=loadVTKPolyData(datadir+"/" + filenames[index]);

      vtkPoints *meshPoints = vtkPoints::New();
      meshPoints=polydata->GetPoints();
      vtkPoints *tmpPoints = vtkPoints::New();
      for( unsigned int pointID = 0; pointID < meshPoints->GetNumberOfPoints(); pointID++ )
      {
        double curPoint[3];
        double avgPoint[3];
        double tmpPoint[3];
        meshPoints->GetPoint(pointID, curPoint);
        avgPoints->GetPoint(pointID, avgPoint);
        for( unsigned int dim = 0; dim < 3; dim++ )
        {
            tmpPoint[dim] = curPoint[dim] + avgPoint[dim];
        }
        tmpPoints->InsertPoint(pointID, tmpPoint);
      }

      avgPoints = tmpPoints;
    }

    vtkPoints *tmpPoints = vtkPoints::New();
    for( unsigned int pointID = 0; pointID < avgPoints->GetNumberOfPoints(); pointID++ )
    {
          double avgPoint[3];
          double tmpPoint[3];
              avgPoints->GetPoint(pointID, avgPoint);
          for( unsigned int dim = 0; dim < 3; dim++ )
          {
            tmpPoint[dim] = avgPoint[dim] / (numMeshes + 1);
          }
          tmpPoints->InsertPoint(pointID, tmpPoint);
    }
    avgPoints = tmpPoints;

    polydata0->SetPoints(avgPoints);

    return polydata0;
}

vtkPoints* subtractMesh (vtkPolyData* mesh1, vtkPolyData* mesh2)
{
    vtkPoints* subtractedPoints = vtkPoints::New();

    vtkPoints* pointsMesh1 = vtkPoints::New();
    pointsMesh1 = mesh1->GetPoints();
    vtkPoints* pointsMesh2 = vtkPoints::New();
    pointsMesh2 = mesh2->GetPoints();

    for( unsigned int pointID = 0; pointID < mesh1->GetNumberOfPoints(); pointID++ )
    {
        double mesh1point[3];
        double mesh2point[3];
        double tmpPoint[3];

        pointsMesh1->GetPoint(pointID, mesh1point);
        pointsMesh2->GetPoint(pointID, mesh2point);

        for( unsigned int dim = 0; dim < 3; dim++ )
        {
            tmpPoint[dim] = mesh1point[dim] - mesh2point[dim];
        }
        subtractedPoints->InsertPoint(pointID, tmpPoint);
    }

    return subtractedPoints;
}

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName.c_str());
    return infile.good();
}

double ComputeOAindex (vtkPoints* diff)
{

    double OAindex;
    double sum = 0;

    for( unsigned int pointID = 0; pointID < diff->GetNumberOfPoints(); pointID++ )
    {
        double point[3];

        diff->GetPoint(pointID, point);

        for( unsigned int dim = 0; dim < 3; dim++ )
        {
            sum = pow ( point[dim],2) + sum ;
        }

    }

    OAindex = sqrt ( sum / (diff->GetNumberOfPoints()*3) );
    return OAindex;
}

int main (int argc, char ** argv)
{

    // Variables that can be changed by the user This will be implemented in the GUI
    int NumOfGroups = 8;
    std::string token = "Data-cog";

    std::string GroupLUT("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/MasterLookUp.csv");

    std::string datadirHC("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/ControlGroup");
    std::string datadirOA("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/OA");
    std::string datadirBoth("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Both");

    std::string datadirG01("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G01");
    std::string datadirG02("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G02");
    std::string datadirG03("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G03");
    std::string datadirG04("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G04");
    std::string datadirG05("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G05");
    std::string datadirG06("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G06");
    std::string datadirG07("/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/" + token + "/Grouping/G07");

    //Read .CSV look up table with group information and VTK
    vtkSmartPointer<vtkDelimitedTextReader> CSVreader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    CSVreader->SetFieldDelimiterCharacters(",");
    CSVreader->SetFileName(GroupLUT.c_str());
    CSVreader->SetHaveHeaders(false);
    CSVreader->Update();
    vtkSmartPointer<vtkTable> table = CSVreader->GetOutput();


    StringVectorType filenamesToClassify;
    getdir(datadirBoth, filenamesToClassify, ".vtk");
    if (filenamesToClassify.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirBoth << " exiting.";
    exit(-1);
    }

    std::cout << "Number of samples to classify " << filenamesToClassify.size() << std::endl;

    VectorType totalGroupVariances (8);
    double totalVarPooled=0;

    Eigen::MatrixXd OAindex_all (NumOfGroups,filenamesToClassify.size());
    VectorType OAindex (filenamesToClassify.size());
    VectorType OAindex_groupAssignment (filenamesToClassify.size());
    VectorType OAindex_groupReal (filenamesToClassify.size());

    // Load all avgs from each group
    vtkPolyData* VTKavgG01 = averageVTK(datadirG01);
    vtkPolyData* VTKavgG02 = averageVTK(datadirG02);
    vtkPolyData* VTKavgG03 = averageVTK(datadirG03);
    vtkPolyData* VTKavgG04 = averageVTK(datadirG04);
    vtkPolyData* VTKavgG05 = averageVTK(datadirG05);
    vtkPolyData* VTKavgG06 = averageVTK(datadirG06);
    vtkPolyData* VTKavgG07 = averageVTK(datadirG07);
    vtkPolyData* VTKavgHC = averageVTK(datadirHC);

    for (unsigned sample = 0 ; sample < filenamesToClassify.size() ; sample++)
    {

        int group = 0;

        for ( int r = 0; r < table->GetNumberOfRows() ; r++ )
        {

            if ( filenamesToClassify[sample] == table->GetValue(r,0).ToString())
            {
                group = atoi(table->GetValue(r,1).ToString().c_str());
            }

         }

        std::cout << sample <<  ": Group assignment for " << filenamesToClassify[sample] << " is " << group << std::endl;
        OAindex_groupReal[sample] = group;


        // ###################################
        //Computing OAindex for all samples
        // ###################################

        // Read ShapeOA in vtk format
        std::string ShapeOAFilename = datadirBoth + "/" + filenamesToClassify[sample];
        vtkPolyData* VTKShapeOA = loadVTKPolyData(ShapeOAFilename);

        vtkPoints* diffG01 = subtractMesh(VTKShapeOA,VTKavgG01);
        vtkPoints* diffG02 = subtractMesh(VTKShapeOA,VTKavgG02);
        vtkPoints* diffG03 = subtractMesh(VTKShapeOA,VTKavgG03);
        vtkPoints* diffG04 = subtractMesh(VTKShapeOA,VTKavgG04);
        vtkPoints* diffG05 = subtractMesh(VTKShapeOA,VTKavgG05);
        vtkPoints* diffG06 = subtractMesh(VTKShapeOA,VTKavgG06);
        vtkPoints* diffG07 = subtractMesh(VTKShapeOA,VTKavgG07);
        vtkPoints* diffHC = subtractMesh(VTKShapeOA,VTKavgHC);

        double minOAindex = 99999999999999;
        double maxOAindex = -1;
        double groupAssignment = 0;

        OAindex_all(0,sample) = ComputeOAindex(diffG01);
        OAindex_all(1,sample) = ComputeOAindex(diffG02);
        OAindex_all(2,sample) = ComputeOAindex(diffG03);
        OAindex_all(3,sample) = ComputeOAindex(diffG04);
        OAindex_all(4,sample) = ComputeOAindex(diffG05);
        OAindex_all(5,sample) = ComputeOAindex(diffG06);
        OAindex_all(6,sample) = ComputeOAindex(diffG07);
        OAindex_all(7,sample) = ComputeOAindex(diffHC);

        for (int i = 0 ; i < NumOfGroups ; i ++ )
        {
            if (OAindex_all(i,sample) < minOAindex)
            {
                minOAindex = OAindex_all(i,sample);
                groupAssignment = i+1;
            }

        }

        OAindex(sample) = minOAindex;
        OAindex_groupAssignment(sample) = groupAssignment;
        std::cout << "--------------------" << std::endl;

    }

    std::cout << "----- Writing OA index results -----" << std::endl;
//
    std::ofstream outputOAindexAll;
    outputOAindexAll.open( "/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/Data/avgAnalysis/OAindexAll_groups_" + token + "scale.csv" );
    outputOAindexAll << "GroupID,";
    for (int ij = 0; ij < filenamesToClassify.size()-1; ij ++)
        outputOAindexAll << filenamesToClassify[ij] << "," ;
    outputOAindexAll << filenamesToClassify[filenamesToClassify.size()-1] << std::endl;

    for (int ij = 0; ij < OAindex_all.rows(); ij ++)
    {
        outputOAindexAll << ij+1 << "," ;
        for (int kl = 0; kl < OAindex_all.cols()-1; kl ++)
        {
            outputOAindexAll << OAindex_all(ij,kl) << "," ;
        }
        outputOAindexAll << OAindex_all(ij,(OAindex_all.cols()-1)) << std::endl;
    }
    outputOAindexAll.close();

    std::ofstream outputOAindex;
    outputOAindex.open( "/Users/bpaniagua/Work/Projects/CMF/TMJR01/OAIndex/Code/Data/avgAnalysis/OAindex_groups_" + token + "scale.csv" );
    outputOAindex << "subjectID,OAindex,GroupAssignment,GroupReal,Classification" << std::endl;

    int missclassifications = 0;
    for (int ij = 0; ij < OAindex.size(); ij ++)
    {
        int missclassified = 0;
        if (OAindex_groupAssignment[ij] != OAindex_groupReal[ij])
        {
            missclassified = 1;
            missclassifications++;
        }
        outputOAindex << filenamesToClassify[ij] << "," << OAindex[ij] << "," << OAindex_groupAssignment[ij] << "," << OAindex_groupReal[ij] << "," << missclassified << std::endl;
    }
    outputOAindex << "Number of missclassifications " << missclassifications << "/" << (missclassifications * 100) / filenamesToClassify.size() << "%" << std::endl;
    outputOAindex.close();


}
