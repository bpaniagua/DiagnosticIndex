/*
* This file is part of the statismo library.
*
* Author: Marcel Luethi (marcel.luethi@unibas.ch)
*
* Copyright (c) 2011 University of Basel
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* Redistributions of source code must retain the above copyright notice,
* this list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
*
* Neither the name of the project's author nor the names of its
* contributors may be used to endorse or promote products derived from
* this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrixtype;

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

vtkPolyData* loadVTKPolyData(const std::string& filename) {
vtkPolyDataReader* reader = vtkPolyDataReader::New();
reader->SetFileName(filename.c_str());
reader->Update();
vtkPolyData* pd = vtkPolyData::New();
pd->ShallowCopy(reader->GetOutput());
return pd;
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

vtkPolyData* VectorToVTK (const VectorType v, vtkPolyData* reference)
{
    unsigned cont_v = 0;

    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
    output->DeepCopy(reference);

    vtkPoints* points = vtkPoints::New();
    points = output->GetPoints();

    for (unsigned i = 0 ; i < output->GetNumberOfPoints() ; i++)
    {
        double p[3];
        p[0]=v(cont_v);
        p[1]=v(cont_v+1);
        p[2]=v(cont_v+2);
        cont_v =  cont_v+3;

        points->InsertPoint(i,p);
    }

    output->SetPoints(points);

    return output;
}

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName.c_str());
    return infile.good();
}

StatisticalModelType* buildSSM (std::string datadir, std::string filenamereference)
{
    StringVectorType filenames;
    getdir(datadir, filenames, ".vtk");
    if (filenames.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadir << " exiting.";
    exit(-1);
    }

    std::cout << "Building SSM of data in " << datadir << std::endl;
    vtkPolyData* reference = loadVTKPolyData(filenamereference);
    boost::scoped_ptr<RepresenterType> representer(RepresenterType::Create(reference));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManager(DataManagerType::Create(representer.get()));

    for (unsigned i = 0; i < filenames.size() ; i++)
    {
        vtkPolyData* dataset = loadVTKPolyData(datadir + "/" + filenames[i]);
        // We provde the filename as a second argument.
        // It will be written as metadata, and allows us to more easily figure out what we did later.
        dataManager->AddDataset(dataset, filenames[i]);
        // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
        dataset->Delete();
    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
    StatisticalModelType* model = modelBuilder->BuildNewModel(dataManager->GetData(), 0.01);

    std::cout << "Total variance " << model->GetPCAVarianceVector().sum() << " number of Eigenmodes " << model->GetPCAVarianceVector().size() << std::endl;
    return model;
}

StatisticalModelType* buildSSMCrossvalidation (std::string datadir, std::string filenamereference, std::string filename)
{
    StringVectorType filenames;
    getdir(datadir, filenames, ".vtk");
    if (filenames.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadir << " exiting.";
    exit(-1);
    }

  //  std::cout << "Building SSM of data in " << datadir << std::endl;
    vtkPolyData* reference = loadVTKPolyData(filenamereference);
    boost::scoped_ptr<RepresenterType> representer(RepresenterType::Create(reference));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManager(DataManagerType::Create(representer.get()));

    for (unsigned i = 0; i < filenames.size() ; i++)
    {
        if (filename.compare(filenames[i]) != 0)
        {
            vtkPolyData* dataset = loadVTKPolyData(datadir + "/" + filenames[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManager->AddDataset(dataset, filenames[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
        }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
    StatisticalModelType* model = modelBuilder->BuildNewModel(dataManager->GetData(), 0.01);

  //  std::cout << "Total variance " << model->GetPCAVarianceVector().sum() << " number of Eigenmodes " << model->GetPCAVarianceVector().size() << std::endl;
    return model;
}

void analyzeVariance (StatisticalModelType* model)
{
    VectorType PCAVariance = model->GetPCAVarianceVector();
    double totalVar = model->GetPCAVarianceVector().sum();
    std::cout << "### Model variance analysis , total variance " << totalVar << "###"<< std::endl;

    int sum =0;
    double per = 0;
    for (int i = 0 ; i < PCAVariance.size() ; i++)
    {
        sum = sum + PCAVariance[i];
        double percentage = (PCAVariance[i] * 100) / totalVar;
        per = per + percentage;
        std::cout << "Eigenvalue " << i << " variance " << PCAVariance[i] << " (" << percentage <<" / " << per << ")" << std::endl;
    }
}

double ComputeOAindex (StatisticalModelType* model, vtkPolyData* sample, double totalVarPooled, int numOfEigenmodes)
{
    VectorType ShapeLoads = model->ComputeCoefficientsForDataset(sample);
    VectorType eigenVectors = model->GetPCAVarianceVector();

    double OAindex;
    double sum = 0;

    double finalnumberofeigenmodes;

    if (numOfEigenmodes == -1)
        finalnumberofeigenmodes = ShapeLoads.size();
    else
        finalnumberofeigenmodes = numOfEigenmodes;

    for (unsigned l = 0 ; l < finalnumberofeigenmodes ; l++)
    {
        sum = sum + ( pow(ShapeLoads(l),2) );
    }

    OAindex = sqrt( sum )/finalnumberofeigenmodes;

    return OAindex;
}

double ComputeOAindex2 (StatisticalModelType* model, vtkPolyData* sample, double totalVarPooled, int numOfEigenmodes)
{

   VectorType ShapeLoads = model->ComputeCoefficientsForDataset(sample);
   VectorType eigenVectors = model->GetPCAVarianceVector();

   double OAindex;
   double sum = 0;

   double finalnumberofeigenmodes;

   if (numOfEigenmodes == -1)
       finalnumberofeigenmodes = ShapeLoads.size();
   else
       finalnumberofeigenmodes = numOfEigenmodes;

   for (unsigned l = 0 ; l < finalnumberofeigenmodes ; l++)
   {
       sum = sum + ( pow( (ShapeLoads(l)* (eigenVectors(l))/totalVarPooled ) ,2) );
   }

   OAindex = sqrt( sum )/finalnumberofeigenmodes;

   return OAindex;

}

int main (int argc, char ** argv)
{

    // Variables that can be changed by the user This will be implemented in the GUI
    int NumOfGroups = 8;
    int NumOfEigenmodes = atoi(argv[1]);

    if( argc < 3 )
    {
        std::cout << "Usage: diagnosticIndex <numberOfEigenModes> <filename>" << std::endl;
        return 0;
    }

    std::string GroupLUT("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/LookUpGroup3.csv");
    std::string datadirModels("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Data/Grouping/HD5Groups");

    std::string datadirHC("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/ControlGroup");
    std::string datadirOA("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/OA");
    std::string datadirBoth("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Both");

    std::string datadirG01("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G01");
    std::string datadirG02("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G02");
    std::string datadirG03("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G03");
    std::string datadirG04("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G04");
    std::string datadirG05("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G05");
    std::string datadirG06("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G06");
    std::string datadirG07("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Grouping/G07");

    std::string datadirRepresenters("/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/DataCO-cog/Representers");

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

    //###########################################
    // Build all SSM for all groups
    //##########################################
    VectorType totalGroupVariances (8);
    double totalVarPooled=0;

    std::string representerG01 = datadirRepresenters + "/avgG01.vtk";
    StatisticalModelType* modelG01 = buildSSM(datadirG01,representerG01 );

    std::string representerG02 = datadirRepresenters + "/avgG02.vtk";
    StatisticalModelType* modelG02 = buildSSM(datadirG02,representerG02);

    std::string representerG03 = datadirRepresenters + "/avgG03.vtk";
    StatisticalModelType* modelG03 = buildSSM(datadirG03,representerG03);

    std::string representerG04 = datadirRepresenters + "/avgG04.vtk";
    StatisticalModelType* modelG04 = buildSSM(datadirG04,representerG04);

    std::string representerG05 = datadirRepresenters + "/avgG05.vtk";
    StatisticalModelType* modelG05 = buildSSM(datadirG05,representerG05);

    std::string representerG06 = datadirRepresenters + "/avgG06.vtk";
    StatisticalModelType* modelG06 = buildSSM(datadirG06,representerG06);

    std::string representerG07 = datadirRepresenters + "/avgG07.vtk";
    StatisticalModelType* modelG07 = buildSSM(datadirG07,representerG07);

    std::string representerHC = datadirRepresenters + "/avgHC.vtk";
    StatisticalModelType* modelHC = buildSSM(datadirHC,representerHC);

   /* std::string representerBoth = datadirRepresenters + "/avgBoth.vtk";
    StatisticalModelType* modelBoth = buildSSM(datadirBoth, representerBoth);
    double totalVarBoth = modelBoth->GetPCAVarianceVector().sum();*/

    //################################################

    MatrixType OAindex_all (NumOfGroups,filenamesToClassify.size());
    MatrixType OAindex_all_idx (NumOfGroups,filenamesToClassify.size());
    VectorType OAindex (filenamesToClassify.size());
    VectorType OAindex_groupAssignment (filenamesToClassify.size());
    VectorType OAindex_groupReal (filenamesToClassify.size());

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


        // ################################################
        //Rebuilding SSM without a sample (leave one out)
        // ################################################
        if ( group == 1 )
        {
            modelG01 = buildSSMCrossvalidation(datadirG01, representerG01, filenamesToClassify[sample]);
        }
        else
        {
            if ( group == 2 )
            {
                modelG02 = buildSSMCrossvalidation(datadirG02, representerG02, filenamesToClassify[sample]);
            }
            else
            {
                if ( group == 3 )
                {
                    modelG03 = buildSSMCrossvalidation(datadirG03, representerG03, filenamesToClassify[sample]);
                }
                else
                {
                    if ( group == 4 )
                    {
                        modelG04 = buildSSMCrossvalidation(datadirG04, representerG04, filenamesToClassify[sample]);
                    }
                    else
                    {
                        if ( group == 5 )
                        {
                            modelG05 = buildSSMCrossvalidation(datadirG05, representerG05, filenamesToClassify[sample]);
                        }
                        else
                        {
                            if ( group == 6 )
                            {
                                modelG06 = buildSSMCrossvalidation(datadirG06, representerG06, filenamesToClassify[sample]);
                            }
                            else
                            {
                                if ( group == 7 )
                                {
                                    modelG07 = buildSSMCrossvalidation(datadirG07, representerG07, filenamesToClassify[sample]);
                                }
                                else
                                {
                                    if ( group == 8 )
                                    {
                                       modelHC = buildSSMCrossvalidation(datadirHC, representerHC, filenamesToClassify[sample]);
                                    }
                                    else
                                    {
                                        std::cerr << "No group information found for " << filenamesToClassify[sample] << std::endl;
                                        exit(-1);
                                    }
                                }

                            }
                        }
                    }
                }
            }
        }

        //Computing variances per group and pooled
        totalGroupVariances(0) = modelG01->GetPCAVarianceVector().sum();
        totalGroupVariances(1) = modelG02->GetPCAVarianceVector().sum();
        totalGroupVariances(2) = modelG03->GetPCAVarianceVector().sum();
        totalGroupVariances(3) = modelG04->GetPCAVarianceVector().sum();
        totalGroupVariances(4) = modelG05->GetPCAVarianceVector().sum();
        totalGroupVariances(5) = modelG06->GetPCAVarianceVector().sum();
        totalGroupVariances(6) = modelG07->GetPCAVarianceVector().sum();
        totalGroupVariances(7) = modelHC->GetPCAVarianceVector().sum();

        double sum = 0;
        for (int i = 0; i <totalGroupVariances.size(); i++)
        {
            sum = sum + totalGroupVariances(i);
        }
        totalVarPooled = sum/totalGroupVariances.size();



        // ###################################
        //Computing OAindex for all samples
        // ###################################

        // Read ShapeOA in vtk format
        std::string ShapeOAFilename = datadirBoth + "/" + filenamesToClassify[sample];
        vtkPolyData* VTKShapeOA = loadVTKPolyData(ShapeOAFilename);

        double minOAindex = 99999999999999;
        double maxOAindex = -1;
        double groupAssignment = 0;

        OAindex_all(0,sample) = ComputeOAindex(modelG01,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(1,sample) = ComputeOAindex(modelG02,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(2,sample) = ComputeOAindex(modelG03,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(3,sample) = ComputeOAindex(modelG04,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(4,sample) = ComputeOAindex(modelG05,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(5,sample) = ComputeOAindex(modelG06,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(6,sample) = ComputeOAindex(modelG07,VTKShapeOA,totalVarPooled,NumOfEigenmodes);
        OAindex_all(7,sample) = ComputeOAindex(modelHC,VTKShapeOA,totalVarPooled,NumOfEigenmodes);

        VectorType col = OAindex_all.col(sample);
        std::sort(col.data(),col.data()+col.size());

        for (int i = 0 ; i < col.size() ; i ++ )
        {
            for (int j = 0 ; j < col.size() ; j ++ )
            {
                if (OAindex_all(i,sample) == col(j))
                {
                    OAindex_all_idx(j,sample)=i+1;
                }
            }
        }

        OAindex(sample) = col(0);
        OAindex_groupAssignment(sample) = OAindex_all_idx(0,sample);
        std::cout << "--------------------" << std::endl;

    }

    std::cout << "----- Writing OA index results -----" << std::endl;
//
    std::ofstream outputOAindexAll;
    outputOAindexAll.open( "/NIRAL/projects5/CMF/TMJR01/OAIndex/Code/Results/OAindexAll_groups_cog.csv" );
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
    outputOAindex.open( argv[2] );
    outputOAindex << "subjectID,OAindex,GroupReal,GroupAssignment,Classification,2ndAssigGroup,3rdAssigGroup,4thAssigGroup,5thAssigGroup" << std::endl;

    int missclassifications = 0;
    for (int ij = 0; ij < OAindex.size(); ij ++)
    {
        int missclassified = 0;
        if (OAindex_groupAssignment[ij] != OAindex_groupReal[ij])
        {
            missclassified = 1;
            missclassifications++;
        }
        outputOAindex << filenamesToClassify[ij] << "," << OAindex[ij] << "," << OAindex_groupReal[ij] <<"," << OAindex_groupAssignment[ij] <<
                         "," << missclassified << "," << OAindex_all_idx(1,ij) << "," << OAindex_all_idx(2,ij) << "," << OAindex_all_idx(3,ij) << "," << OAindex_all_idx(4,ij)<< std::endl;
    }
    outputOAindex << "Number of missclassifications " << missclassifications << "/" << (missclassifications * 100) / filenamesToClassify.size() << "%" << std::endl;
    outputOAindex.close();


}
