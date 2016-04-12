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

int main (int argc, char ** argv)
{
    std::string GroupLUT("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/LookUpGroup.csv");
    std::string GroupLUTtest("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/LookUpGroup-testorient.csv");

    std::string datadirHC("/projects5/CMF/TMJR01/OAIndex/JP/Data/ControlGroup");
    std::string datadirOA("/projects5/CMF/TMJR01/OAIndex/JP/Data/OA");
    std::string datadirBoth("/projects5/CMF/TMJR01/OAIndex/JP/Data/Both");

    std::string datadirG01("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G01");
    std::string datadirG02("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G02");
    std::string datadirG03("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G03");
    std::string datadirG04("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G04");
    std::string datadirG05("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G05");
    std::string datadirG06("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G06");
    std::string datadirG07("/NIRAL/projects5/CMF/TMJR01/OAIndex/JP/Grouping/G07");

    // All the statismo classes have to be parameterized with the RepresenterType.
    // For building a shape model with vtk, we use the vtkPolyDataRepresenter.

    StringVectorType filenamesHC;
    getdir(datadirHC, filenamesHC, ".vtk");
    if (filenamesHC.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirHC << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesOA;
    getdir(datadirOA, filenamesOA, ".vtk");
    if (filenamesOA.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirOA << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesBoth;
    getdir(datadirBoth, filenamesBoth, ".vtk");
    if (filenamesBoth.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirBoth << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG01;
    getdir(datadirG01, filenamesG01, ".vtk");
    if (filenamesG01.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG01 << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG02;
    getdir(datadirG02, filenamesG02, ".vtk");
    if (filenamesG02.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG02 << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG03;
    getdir(datadirG03, filenamesG03, ".vtk");
    if (filenamesG01.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG01 << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG04;
    getdir(datadirG04, filenamesG04, ".vtk");
    if (filenamesG04.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG04 << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG05;
    getdir(datadirG05, filenamesG05, ".vtk");
    if (filenamesG05.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG05 << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG06;
    getdir(datadirG06, filenamesG06, ".vtk");
    if (filenamesG06.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG06 << " exiting.";
    exit(-1);
    }

    StringVectorType filenamesG07;
    getdir(datadirG07, filenamesG07, ".vtk");
    if (filenamesG07.size() == 0) {
        std::cerr << "did not find any vtk files in directory " << datadirG07 << " exiting.";
    exit(-1);
    }

    //Read .CSV look up table with group information and VTK
    vtkSmartPointer<vtkDelimitedTextReader> CSVreader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    CSVreader->SetFieldDelimiterCharacters(",");
    CSVreader->SetFileName(GroupLUTtest.c_str());
    CSVreader->SetHaveHeaders(true);
    CSVreader->Update();
    vtkSmartPointer<vtkTable> table = CSVreader->GetOutput();


   // std::cout << "Number of samples to classify " << filenamesBoth.size() << " -- OA = " << filenamesOA.size() << "  -- HC = "
     //         << filenamesHC.size() << std::endl;

    std::cout << "Number of samples to classify " << filenamesOA.size() << std::endl;

    //###########################################
    // Build all SSM for all groups
    //##########################################
    std::cout << "Building SSM of G01" << std::endl;
    vtkPolyData* referenceG01 = loadVTKPolyData(datadirG01 + "/" + filenamesG01[0]);
    boost::scoped_ptr<RepresenterType> representerG01(RepresenterType::Create(referenceG01));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG01(DataManagerType::Create(representerG01.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG01 (referenceG01->GetNumberOfPoints()*3, filenamesG01.size()-1);
    for (unsigned i = 0; i < filenamesG01.size() ; i++)
    {
        //if (i != sample)
        //{
            vtkPolyData* dataset = loadVTKPolyData(datadirG01 + "/" + filenamesG01[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG01->AddDataset(dataset, filenamesG01[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
        //}

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
    StatisticalModelType* modelG01 = modelBuilder->BuildNewModel(dataManagerG01->GetData(), 0);

    //################################################
    std::cout << "Building SSM of G02" << std::endl;
    vtkPolyData* referenceG02 = loadVTKPolyData(datadirG02 + "/" + filenamesG02[0]);
    boost::scoped_ptr<RepresenterType> representerG02(RepresenterType::Create(referenceG02));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG02(DataManagerType::Create(representerG02.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG02 (referenceG02->GetNumberOfPoints()*3, filenamesG02.size()-1);
    for (unsigned i = 0; i < filenamesG02.size() ; i++)
    {
//                        if (i != sample)
//                      {
            vtkPolyData* dataset = loadVTKPolyData(datadirG02 + "/" + filenamesG02[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG02->AddDataset(dataset, filenamesG02[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
//                    }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelG02 = modelBuilder->BuildNewModel(dataManagerG02->GetData(), 0);

    //################################################

    std::cout << "Building SSM of G03" << std::endl;
    vtkPolyData* referenceG03 = loadVTKPolyData(datadirG03 + "/" + filenamesG03[0]);
    boost::scoped_ptr<RepresenterType> representerG03(RepresenterType::Create(referenceG03));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG03(DataManagerType::Create(representerG03.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG03 (referenceG03->GetNumberOfPoints()*3, filenamesG03.size()-1);
    for (unsigned i = 0; i < filenamesG03.size() ; i++)
    {
//                        if (i != sample)
//                      {
            vtkPolyData* dataset = loadVTKPolyData(datadirG03 + "/" + filenamesG03[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG03->AddDataset(dataset, filenamesG03[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
//                    }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelG03 = modelBuilder->BuildNewModel(dataManagerG03->GetData(), 0);

    //################################################

    std::cout << "Building SSM of G04" << std::endl;
    vtkPolyData* referenceG04 = loadVTKPolyData(datadirG04 + "/" + filenamesG04[0]);
    boost::scoped_ptr<RepresenterType> representerG04(RepresenterType::Create(referenceG04));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG04(DataManagerType::Create(representerG04.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG04 (referenceG04->GetNumberOfPoints()*3, filenamesG04.size()-1);
    for (unsigned i = 0; i < filenamesG04.size() ; i++)
    {
//                        if (i != sample)
//                      {
            vtkPolyData* dataset = loadVTKPolyData(datadirG04 + "/" + filenamesG04[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG04->AddDataset(dataset, filenamesG04[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
//                    }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelG04 = modelBuilder->BuildNewModel(dataManagerG04->GetData(), 0);

    //################################################

    std::cout << "Building SSM of G05" << std::endl;
    vtkPolyData* referenceG05 = loadVTKPolyData(datadirG05 + "/" + filenamesG05[0]);
    boost::scoped_ptr<RepresenterType> representerG05(RepresenterType::Create(referenceG05));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG05(DataManagerType::Create(representerG05.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG05 (referenceG05->GetNumberOfPoints()*3, filenamesG05.size()-1);
    for (unsigned i = 0; i < filenamesG05.size() ; i++)
    {
//                        if (i != sample)
//                      {
            vtkPolyData* dataset = loadVTKPolyData(datadirG05 + "/" + filenamesG05[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG05->AddDataset(dataset, filenamesG05[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
//                    }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelG05 = modelBuilder->BuildNewModel(dataManagerG05->GetData(), 0);

    //################################################

    std::cout << "Building SSM of G06" << std::endl;
    vtkPolyData* referenceG06 = loadVTKPolyData(datadirG06 + "/" + filenamesG06[0]);
    boost::scoped_ptr<RepresenterType> representerG06(RepresenterType::Create(referenceG06));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG06(DataManagerType::Create(representerG06.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG06 (referenceG06->GetNumberOfPoints()*3, filenamesG06.size()-1);
    for (unsigned i = 0; i < filenamesG06.size() ; i++)
    {
//                        if (i != sample)
//                      {
            vtkPolyData* dataset = loadVTKPolyData(datadirG06 + "/" + filenamesG06[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG06->AddDataset(dataset, filenamesG06[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
//                    }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelG06 = modelBuilder->BuildNewModel(dataManagerG06->GetData(), 0);

    //################################################

    std::cout << "Building SSM of G07" << std::endl;
    vtkPolyData* referenceG07 = loadVTKPolyData(datadirG07 + "/" + filenamesG07[0]);
    boost::scoped_ptr<RepresenterType> representerG07(RepresenterType::Create(referenceG07));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerG07(DataManagerType::Create(representerG07.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixG07 (referenceG07->GetNumberOfPoints()*3, filenamesG07.size()-1);
    for (unsigned i = 0; i < filenamesG07.size() ; i++)
    {
//                        if (i != sample)
//                      {
            vtkPolyData* dataset = loadVTKPolyData(datadirG07 + "/" + filenamesG07[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerG07->AddDataset(dataset, filenamesG07[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
//                    }

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelG07 = modelBuilder->BuildNewModel(dataManagerG07->GetData(), 0);

    //###########################################
    // Build all SSM for healthy
    //##########################################
    std::cout << "Building SSM of Healthy" << std::endl;
    vtkPolyData* referenceHC = loadVTKPolyData(datadirHC + "/" + filenamesHC[0]);
    boost::scoped_ptr<RepresenterType> representerHC(RepresenterType::Create(referenceHC));
    // We create a datamanager and provide it with a pointer to the representer
    boost::scoped_ptr<DataManagerType> dataManagerHC(DataManagerType::Create(representerHC.get()));
    // Now we add our data to the data manager
    // load the data and add it to the data manager.
    MatrixType dataMatrixHC (referenceHC->GetNumberOfPoints()*3, filenamesHC.size()-1);
    for (unsigned i = 0; i < filenamesHC.size() ; i++)
    {
        //if (i != sample)
        //{
            vtkPolyData* dataset = loadVTKPolyData(datadirHC + "/" + filenamesHC[i]);
            // We provde the filename as a second argument.
            // It will be written as metadata, and allows us to more easily figure out what we did later.
            dataManagerHC->AddDataset(dataset, filenamesHC[i]);
            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
            dataset->Delete();
        //}

    }
    // To actually build a model, we need to create a model builder object.
    // Calling the build model with a list of samples from the data manager, returns a new model.
    // The second parameter to BuildNewModel is the variance of the noise on our data
    StatisticalModelType* modelHC = modelBuilder->BuildNewModel(dataManagerHC->GetData(), 0);

    //################################################

    Eigen::MatrixXd OAindex_all (8,filenamesOA.size());
    VectorType OAindex (filenamesOA.size());
    VectorType OAindex_groupAssignment (filenamesOA.size());
    VectorType OAindex_groupReal (filenamesOA.size());

    for (unsigned sample = 0 ; sample < filenamesOA.size() ; sample++)
    {

        int group = 0;

        for ( int r = 0; r < table->GetNumberOfRows() ; r++ )
        {

            if ( filenamesOA[sample] == table->GetValue(r,0).ToString())
            {
                group = atoi(table->GetValue(r,1).ToString().c_str());
            }

         }

        std::cout << sample <<  ": Group assignment for " << filenamesOA[sample] << " is " << group << std::endl;
        OAindex_groupReal[sample] = group;


        // ################################################
        //Rebuilding SSM without a sample (leave one out)
        // ################################################
        if ( group == 1 )
        {
            std::cout << "Re-Building SSM of G01 without sample " << filenamesOA[sample] << std::endl;
            if ( filenamesG01[0] == filenamesOA[sample] )
            {
                vtkPolyData* referenceG01 = loadVTKPolyData(datadirG01 + "/" + filenamesG01[1]);
            }
            else
            {
                vtkPolyData* referenceG01 = loadVTKPolyData(datadirG01 + "/" + filenamesG01[0]);
            }
            boost::scoped_ptr<RepresenterType> representerG01(RepresenterType::Create(referenceG01));
            // We create a datamanager and provide it with a pointer to the representer
            boost::scoped_ptr<DataManagerType> dataManagerG01(DataManagerType::Create(representerG01.get()));
            // Now we add our data to the data manager
            // load the data and add it to the data manager.
            MatrixType dataMatrixG01 (referenceG01->GetNumberOfPoints()*3, filenamesG01.size()-1);
            for (unsigned i = 0; i < filenamesG01.size() ; i++)
            {
                if (i != sample)
                {
                    vtkPolyData* dataset = loadVTKPolyData(datadirG01 + "/" + filenamesG01[i]);
                    // We provde the filename as a second argument.
                    // It will be written as metadata, and allows us to more easily figure out what we did later.
                    dataManagerG01->AddDataset(dataset, filenamesG01[i]);
                    // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                    dataset->Delete();
                }

            }
            // To actually build a model, we need to create a model builder object.
            // Calling the build model with a list of samples from the data manager, returns a new model.
            // The second parameter to BuildNewModel is the variance of the noise on our data
            ModelBuilderType* modelBuilder = ModelBuilderType::Create();
            modelG01 = modelBuilder->BuildNewModel(dataManagerG01->GetData(), 0);

        }
        else
        {
            if ( group == 2 )
            {
                std::cout << "Re-Building SSM of G02 without sample " << filenamesOA[sample] << std::endl;
                if ( filenamesG02[0] == filenamesOA[sample] )
                {
                    vtkPolyData* referenceG02 = loadVTKPolyData(datadirG02 + "/" + filenamesG02[1]);
                }
                else
                {
                    vtkPolyData* referenceG02 = loadVTKPolyData(datadirG02 + "/" + filenamesG02[0]);
                }
                boost::scoped_ptr<RepresenterType> representerG02(RepresenterType::Create(referenceG02));
                // We create a datamanager and provide it with a pointer to the representer
                boost::scoped_ptr<DataManagerType> dataManagerG02(DataManagerType::Create(representerG02.get()));
                // Now we add our data to the data manager
                // load the data and add it to the data manager.
                MatrixType dataMatrixG02 (referenceG02->GetNumberOfPoints()*3, filenamesG02.size()-1);
                for (unsigned i = 0; i < filenamesG02.size() ; i++)
                {
                    if (i != sample)
                    {
                        vtkPolyData* dataset = loadVTKPolyData(datadirG02 + "/" + filenamesG02[i]);
                        // We provde the filename as a second argument.
                        // It will be written as metadata, and allows us to more easily figure out what we did later.
                        dataManagerG02->AddDataset(dataset, filenamesG02[i]);
                        // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                        dataset->Delete();
                    }

                }
                // To actually build a model, we need to create a model builder object.
                // Calling the build model with a list of samples from the data manager, returns a new model.
                // The second parameter to BuildNewModel is the variance of the noise on our data
                ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                modelG02 = modelBuilder->BuildNewModel(dataManagerG02->GetData(), 0);
            }
            else
            {
                if ( group == 3 )
                {
                    std::cout << "Re-Building SSM of G03 without sample " << filenamesOA[sample] << std::endl;
                    if ( filenamesG03[0] == filenamesOA[sample] )
                    {
                        vtkPolyData* referenceG03 = loadVTKPolyData(datadirG03 + "/" + filenamesG03[1]);
                    }
                    else
                    {
                        vtkPolyData* referenceG03 = loadVTKPolyData(datadirG03 + "/" + filenamesG03[0]);
                    }
                    boost::scoped_ptr<RepresenterType> representerG03(RepresenterType::Create(referenceG03));
                    // We create a datamanager and provide it with a pointer to the representer
                    boost::scoped_ptr<DataManagerType> dataManagerG03(DataManagerType::Create(representerG03.get()));
                    // Now we add our data to the data manager
                    // load the data and add it to the data manager.
                    MatrixType dataMatrixG03 (referenceG03->GetNumberOfPoints()*3, filenamesG03.size()-1);
                    for (unsigned i = 0; i < filenamesG03.size() ; i++)
                    {
                        if (i != sample)
                        {
                            vtkPolyData* dataset = loadVTKPolyData(datadirG03 + "/" + filenamesG03[i]);
                            // We provde the filename as a second argument.
                            // It will be written as metadata, and allows us to more easily figure out what we did later.
                            dataManagerG03->AddDataset(dataset, filenamesG03[i]);
                            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                            dataset->Delete();
                        }

                    }
                    // To actually build a model, we need to create a model builder object.
                    // Calling the build model with a list of samples from the data manager, returns a new model.
                    // The second parameter to BuildNewModel is the variance of the noise on our data
                    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                    modelG03 = modelBuilder->BuildNewModel(dataManagerG03->GetData(), 0);
                }
                else
                {
                    if ( group == 4 )
                    {
                        std::cout << "Re-Building SSM of G04 without sample " << filenamesOA[sample] << std::endl;
                        if ( filenamesG04[0] == filenamesOA[sample] )
                        {
                            vtkPolyData* referenceG04 = loadVTKPolyData(datadirG04 + "/" + filenamesG04[1]);
                        }
                        else
                        {
                            vtkPolyData* referenceG04 = loadVTKPolyData(datadirG04 + "/" + filenamesG04[0]);
                        }
                        boost::scoped_ptr<RepresenterType> representerG04(RepresenterType::Create(referenceG04));
                        // We create a datamanager and provide it with a pointer to the representer
                        boost::scoped_ptr<DataManagerType> dataManagerG04(DataManagerType::Create(representerG04.get()));
                        // Now we add our data to the data manager
                        // load the data and add it to the data manager.
                        MatrixType dataMatrixG04 (referenceG04->GetNumberOfPoints()*3, filenamesG04.size()-1);
                        for (unsigned i = 0; i < filenamesG04.size() ; i++)
                        {
                            if (i != sample)
                            {
                                vtkPolyData* dataset = loadVTKPolyData(datadirG04 + "/" + filenamesG04[i]);
                                // We provde the filename as a second argument.
                                // It will be written as metadata, and allows us to more easily figure out what we did later.
                                dataManagerG04->AddDataset(dataset, filenamesG04[i]);
                                // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                                dataset->Delete();
                            }

                        }
                        // To actually build a model, we need to create a model builder object.
                        // Calling the build model with a list of samples from the data manager, returns a new model.
                        // The second parameter to BuildNewModel is the variance of the noise on our data
                        ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                        modelG04 = modelBuilder->BuildNewModel(dataManagerG04->GetData(), 0);
                    }
                    else
                    {
                        if ( group == 5 )
                        {
                            std::cout << "Re-Building SSM of G05 without sample " << filenamesOA[sample] << std::endl;
                            if ( filenamesG05[0] == filenamesOA[sample] )
                            {
                                vtkPolyData* referenceG05 = loadVTKPolyData(datadirG05 + "/" + filenamesG05[1]);
                            }
                            else
                            {
                                vtkPolyData* referenceG05 = loadVTKPolyData(datadirG05 + "/" + filenamesG05[0]);
                            }
                            boost::scoped_ptr<RepresenterType> representerG05(RepresenterType::Create(referenceG05));
                            // We create a datamanager and provide it with a pointer to the representer
                            boost::scoped_ptr<DataManagerType> dataManagerG05(DataManagerType::Create(representerG05.get()));
                            // Now we add our data to the data manager
                            // load the data and add it to the data manager.
                            MatrixType dataMatrixG05 (referenceG05->GetNumberOfPoints()*3, filenamesG05.size()-1);
                            for (unsigned i = 0; i < filenamesG05.size() ; i++)
                            {
                                if (i != sample)
                                {
                                    vtkPolyData* dataset = loadVTKPolyData(datadirG05 + "/" + filenamesG05[i]);
                                    // We provde the filename as a second argument.
                                    // It will be written as metadata, and allows us to more easily figure out what we did later.
                                    dataManagerG05->AddDataset(dataset, filenamesG05[i]);
                                    // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                                    dataset->Delete();
                                }

                            }
                            // To actually build a model, we need to create a model builder object.
                            // Calling the build model with a list of samples from the data manager, returns a new model.
                            // The second parameter to BuildNewModel is the variance of the noise on our data
                            ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                            modelG05 = modelBuilder->BuildNewModel(dataManagerG05->GetData(), 0);
                        }
                        else
                        {
                            if ( group == 6 )
                            {
                                std::cout << "Re-Building SSM of G06 without sample " << filenamesOA[sample] << std::endl;
                                if ( filenamesG06[0] == filenamesOA[sample] )
                                {
                                    vtkPolyData* referenceG06 = loadVTKPolyData(datadirG06 + "/" + filenamesG06[1]);
                                }
                                else
                                {
                                    vtkPolyData* referenceG06 = loadVTKPolyData(datadirG06 + "/" + filenamesG06[0]);
                                }
                                boost::scoped_ptr<RepresenterType> representerG06(RepresenterType::Create(referenceG06));
                                // We create a datamanager and provide it with a pointer to the representer
                                boost::scoped_ptr<DataManagerType> dataManagerG06(DataManagerType::Create(representerG06.get()));
                                // Now we add our data to the data manager
                                // load the data and add it to the data manager.
                                MatrixType dataMatrixG06 (referenceG06->GetNumberOfPoints()*3, filenamesG06.size()-1);
                                for (unsigned i = 0; i < filenamesG06.size() ; i++)
                                {
                                    if (i != sample)
                                    {
                                        vtkPolyData* dataset = loadVTKPolyData(datadirG06 + "/" + filenamesG06[i]);
                                        // We provde the filename as a second argument.
                                        // It will be written as metadata, and allows us to more easily figure out what we did later.
                                        dataManagerG06->AddDataset(dataset, filenamesG06[i]);
                                        // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                                        dataset->Delete();
                                    }

                                }
                                // To actually build a model, we need to create a model builder object.
                                // Calling the build model with a list of samples from the data manager, returns a new model.
                                // The second parameter to BuildNewModel is the variance of the noise on our data
                                ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                                modelG06 = modelBuilder->BuildNewModel(dataManagerG06->GetData(), 0);
                            }
                            else
                            {
                                if ( group == 7 )
                                {
                                    std::cout << "Re-Building SSM of G07 without sample " << filenamesOA[sample] << std::endl;
                                    if ( filenamesG07[0] == filenamesOA[sample] )
                                    {
                                        vtkPolyData* referenceG07 = loadVTKPolyData(datadirG07 + "/" + filenamesG07[1]);
                                    }
                                    else
                                    {
                                        vtkPolyData* referenceG07 = loadVTKPolyData(datadirG07 + "/" + filenamesG07[0]);
                                    }
                                    boost::scoped_ptr<RepresenterType> representerG07(RepresenterType::Create(referenceG07));
                                    // We create a datamanager and provide it with a pointer to the representer
                                    boost::scoped_ptr<DataManagerType> dataManagerG07(DataManagerType::Create(representerG07.get()));
                                    // Now we add our data to the data manager
                                    // load the data and add it to the data manager.
                                    MatrixType dataMatrixG07 (referenceG07->GetNumberOfPoints()*3, filenamesG07.size()-1);
                                    for (unsigned i = 0; i < filenamesG07.size() ; i++)
                                    {
                                        if (i != sample)
                                        {
                                            vtkPolyData* dataset = loadVTKPolyData(datadirG07 + "/" + filenamesG07[i]);
                                            // We provde the filename as a second argument.
                                            // It will be written as metadata, and allows us to more easily figure out what we did later.
                                            dataManagerG07->AddDataset(dataset, filenamesG07[i]);
                                            // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                                            dataset->Delete();
                                        }

                                    }
                                    // To actually build a model, we need to create a model builder object.
                                    // Calling the build model with a list of samples from the data manager, returns a new model.
                                    // The second parameter to BuildNewModel is the variance of the noise on our data
                                    ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                                    modelG07 = modelBuilder->BuildNewModel(dataManagerG07->GetData(), 0);
                                }
                                else
                                {
                                    if ( group == 8 )
                                    {
                                        std::cout << "Re-Building SSM of HC without sample " << filenamesOA[sample] << std::endl;
                                        if ( filenamesHC[0] == filenamesOA[sample] )
                                        {
                                            vtkPolyData* referenceHC = loadVTKPolyData(datadirHC + "/" + filenamesHC[1]);
                                        }
                                        else
                                        {
                                            vtkPolyData* referenceHC = loadVTKPolyData(datadirHC + "/" + filenamesHC[0]);
                                        }
                                        boost::scoped_ptr<RepresenterType> representerHC(RepresenterType::Create(referenceHC));
                                        // We create a datamanager and provide it with a pointer to the representer
                                        boost::scoped_ptr<DataManagerType> dataManagerHC(DataManagerType::Create(representerHC.get()));
                                        // Now we add our data to the data manager
                                        // load the data and add it to the data manager.
                                        MatrixType dataMatrixHC (referenceHC->GetNumberOfPoints()*3, filenamesHC.size()-1);
                                        for (unsigned i = 0; i < filenamesHC.size() ; i++)
                                        {
                                            if (i != sample)
                                            {
                                                vtkPolyData* dataset = loadVTKPolyData(datadirHC + "/" + filenamesHC[i]);
                                                // We provde the filename as a second argument.
                                                // It will be written as metadata, and allows us to more easily figure out what we did later.
                                                dataManagerHC->AddDataset(dataset, filenamesHC[i]);
                                                // it is save to delete the dataset after it was added, as the datamanager direclty copies it.
                                                dataset->Delete();
                                            }

                                        }
                                        // To actually build a model, we need to create a model builder object.
                                        // Calling the build model with a list of samples from the data manager, returns a new model.
                                        // The second parameter to BuildNewModel is the variance of the noise on our data
                                        ModelBuilderType* modelBuilder = ModelBuilderType::Create();
                                        modelHC = modelBuilder->BuildNewModel(dataManagerHC->GetData(), 0);
                                    }
                                    else
                                    {
                                        std::cerr << "No group information found for " << filenamesOA[sample] << std::endl;
                                        exit(-1);
                                    }
                                }

                            }
                        }
                    }
                }
            }
        }

        // ###################################
        //Computing OAindex for all samples
        // ###################################

        // Read ShapeOA in vtk format
        vtkPolyDataReader* reader = vtkPolyDataReader::New();
        std::string ShapeOAFilename = datadirBoth + "/" + filenamesOA[sample];
        reader->SetFileName(ShapeOAFilename.c_str());
        reader->Update();
        vtkPolyData* VTKShapeOA = vtkPolyData::New();
        VTKShapeOA->ShallowCopy(reader->GetOutput());

        VectorType ShapeOAVectorLoadsG01 = modelG01->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsG02 = modelG02->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsG03 = modelG03->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsG04 = modelG04->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsG05 = modelG05->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsG06 = modelG06->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsG07 = modelG07->ComputeCoefficientsForDataset(VTKShapeOA);
        VectorType ShapeOAVectorLoadsHC = modelHC->ComputeCoefficientsForDataset(VTKShapeOA);

        double minOAindex = 99999999999999;
        double groupAssignment = 0;
        double sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG01.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG01(l),2);
        }

        OAindex_all(0,sample) = sqrt(sum);
        if (OAindex_all(0,sample) < minOAindex)
        {
            minOAindex = OAindex_all(0,sample);
            groupAssignment = 1;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG02.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG02(l),2);
        }

        OAindex_all(1,sample) = sqrt(sum);
        if (OAindex_all(1,sample) < minOAindex)
        {
            minOAindex = OAindex_all(1,sample);
            groupAssignment = 2;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG03.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG03(l),2);
        }

        OAindex_all(2,sample) = sqrt(sum);
        if (OAindex_all(2,sample) < minOAindex)
        {
            minOAindex = OAindex_all(2,sample);
            groupAssignment = 3;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG04.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG04(l),2);
        }

        OAindex_all(3,sample) = sqrt(sum);
        if (OAindex_all(3,sample) < minOAindex)
        {
            minOAindex = OAindex_all(3,sample);
            groupAssignment = 4;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG05.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG05(l),2);
        }

        OAindex_all(4,sample) = sqrt(sum);
        if (OAindex_all(4,sample) < minOAindex)
        {
            minOAindex = OAindex_all(4,sample);
            groupAssignment = 5;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG06.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG06(l),2);
        }

        OAindex_all(5,sample) = sqrt(sum);
        if (OAindex_all(5,sample) < minOAindex)
        {
            minOAindex = OAindex_all(5,sample);
            groupAssignment = 6;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsG07.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsG07(l),2);
        }

        OAindex_all(6,sample) = sqrt(sum);
        if (OAindex_all(6,sample) < minOAindex)
        {
            minOAindex = OAindex_all(6,sample);
            groupAssignment = 7;
        }

        sum = 0;
        for (unsigned l = 0 ; l < ShapeOAVectorLoadsHC.size() ; l++)
        {
            sum = sum + pow(ShapeOAVectorLoadsHC(l),2);
        }

        OAindex_all(7,sample) = sqrt(sum);
        if (OAindex_all(7,sample) < minOAindex)
        {
            minOAindex = OAindex_all(7,sample);
            groupAssignment = 8;
        }


        OAindex(sample) = minOAindex;
        OAindex_groupAssignment(sample) = groupAssignment;
        std::cout << "--------------------" << std::endl;
    }

    std::cout << "----- Writing OA index results -----" << std::endl;

    std::ofstream outputOAindexAll;
    outputOAindexAll.open( "/projects5/CMF/TMJR01/OAIndex/JP/Results/OAindexAll_groups.csv" );
    outputOAindexAll << "GroupID,";
    for (int ij = 0; ij < filenamesOA.size()-1; ij ++)
        outputOAindexAll << filenamesOA[ij] << "," ;
    outputOAindexAll << filenamesOA[filenamesOA.size()-1] << std::endl;

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
    outputOAindex.open( "/projects5/CMF/TMJR01/OAIndex/JP/Results/OAindex_groups.csv" );
    outputOAindex << "subjectID,OAindex,GroupAssignment,GroupReal,Classification" << std::endl;

    for (int ij = 0; ij < OAindex.size(); ij ++)
    {
        int missclassified = 0;
        if (OAindex_groupAssignment[ij] != OAindex_groupReal[ij])
        {
            missclassified = 1;
        }
        outputOAindex << filenamesOA[ij] << "," << OAindex[ij] << "," << OAindex_groupAssignment[ij] << "," << OAindex_groupReal[ij] << "," << missclassified << std::endl;
    }
    outputOAindex.close();


}
