/* Copyright (C) 2004-2006 by Gerard Gorman
   Copyright (C) 2006- Imperial College London and others.
   
   Please see the AUTHORS file in the main source directory for a full
   list of copyright holders.
   
   Dr Gerard J Gorman
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London
   
   g.gorman@imperial.ac.uk
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
*/

/* Copyright (C) 2004-2006 by Gerard Gorman
   Copyright (C) 2006- Imperial College London and others.
   
   Please see the AUTHORS file in the main source directory for a full
   list of copyright holders.
   
   Dr Gerard J Gorman
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London
   
   g.gorman@imperial.ac.uk
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
*/

#include "confdefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <vtkMPIController.h>
#else
#include <vtkDummyController.h>
#endif

#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <vector>
#include <string>
#include <cassert>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include "tinyxml.h"

using namespace std;

static vtkUnstructuredGrid *dataSet;
static string fl_vtkFileName;
static unsigned ncnt;
static unsigned ecnt;

int pvtu_search_and_replace(TiXmlElement *pElement, const char *dir){
  if (!pElement) return 0;
  
  TiXmlAttribute* pAttrib=pElement->FirstAttribute();
  
  while (pAttrib){
    if(std::string(pAttrib->Name())=="Source"){
      pAttrib->SetValue(std::string(dir)+pAttrib->Value());
    }
    pAttrib=pAttrib->Next();
  }

  return 0;	
}

void pvtu_search_and_replace(TiXmlNode* pParent, const char *dir){
  if (!pParent) return;
  
  TiXmlNode* pChild;
  int t = pParent->Type();
  
  switch (t){
  case TiXmlNode::ELEMENT:
    if(std::string(pParent->Value())=="Piece"){
      pvtu_search_and_replace(pParent->ToElement(), dir);
    }
    break;
  default:
    break;
  }
  
  for(pChild = pParent->FirstChild(); pChild!=0; pChild=pChild->NextSibling()){
    pvtu_search_and_replace(pChild, dir);
  }
}

void pvtu_fix_path(const char* pFilename, const char *dir){
  TiXmlDocument doc(pFilename);
  bool loadOkay = doc.LoadFile();
  if(loadOkay){
    pvtu_search_and_replace(&doc, dir);
  }else{
    cerr<<"ERROR: Failed to load file "<<pFilename<<endl;
  }
  doc.SaveFile();
}

extern "C" {
  void vtkopen(char *outName, int *len1, char *vtkTitle, int *len2){
    static bool initialized = false;
    if(initialized) 
      assert(dataSet==NULL);
    initialized = true;
    
    dataSet = vtkUnstructuredGrid::New();
    string title(vtkTitle, *len2);
    fl_vtkFileName = string(outName, *len1);
    
    return;
  }

  void vtkwritemesh(int *NNodes, int *NElems, 
                   float *x, float *y, float *z,
                   int *enlist, int *elementTypes, int *elementSizes){
    ncnt = *NNodes;
    ecnt = *NElems;

    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToFloat();
    for(unsigned i=0; i<ncnt; i++){
      float xyz[3];
      xyz[0] = x[i]; 
      xyz[1] = y[i]; 
      xyz[2] = z[i]; 
    
      newPts->InsertPoint(i, xyz);
    }
    dataSet->SetPoints(newPts);

    vtkIdType cell[20];
    int *elem = enlist;
    dataSet->Allocate(ecnt);
    for(unsigned i=0; i<ecnt; i++){
      if(elementTypes[i] == 9){
        cell[0] = elem[0]-1;
        cell[1] = elem[1]-1;
        cell[2] = elem[3]-1; 
        cell[3] = elem[2]-1;
      }else if(elementTypes[i] == 12){
        cell[0] = elem[0]-1;
        cell[1] = elem[1]-1;
        cell[2] = elem[3]-1;
        cell[3] = elem[2]-1;
        cell[4] = elem[4]-1;
        cell[5] = elem[5]-1;
        cell[6] = elem[7]-1;
        cell[7] = elem[6]-1;
      }else{
        for(int j=0; j<elementSizes[i]; j++)
          cell[j] = elem[j]-1;
      }
    
      dataSet->InsertNextCell(elementTypes[i], elementSizes[i], cell);    
      elem+=elementSizes[i];
    }

    newPts->Delete();  
    return;
  }

  void vtkwritemeshd(int *NNodes, int *NElems,
                    double *x, double *y, double *z,
                    int *enlist, int *elementTypes, int *elementSizes){
    ncnt = *NNodes;
    ecnt = *NElems;

    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    for(unsigned i=0; i<ncnt; i++){
      double xyz[3];
      xyz[0] = x[i]; 
      xyz[1] = y[i]; 
      xyz[2] = z[i]; 
    
      newPts->InsertPoint(i, xyz);
    }
    dataSet->SetPoints(newPts);

    vtkIdType cell[20];
    int *elem = enlist;
    dataSet->Allocate(ecnt);
    for(unsigned i=0; i<ecnt; i++){
      if(elementTypes[i] == 9){
        cell[0] = elem[0]-1;
        cell[1] = elem[1]-1;
        cell[2] = elem[3]-1; 
        cell[3] = elem[2]-1;
      }else if(elementTypes[i] == 12){
        cell[0] = elem[0]-1;
        cell[1] = elem[1]-1;
        cell[2] = elem[3]-1;
        cell[3] = elem[2]-1;
        cell[4] = elem[4]-1;
        cell[5] = elem[5]-1;
        cell[6] = elem[7]-1;
        cell[7] = elem[6]-1;
      }else{
        for(int j=0; j<elementSizes[i]; j++)
          cell[j] = elem[j]-1;
      }
    
      dataSet->InsertNextCell(elementTypes[i], elementSizes[i], cell);    
      elem+=elementSizes[i];
    }

    newPts->Delete();  
    return;
  }
  
  void vtkstartn(){}

  void vtkwriteisn(int *vect, char *name, int *len){
    string tag(name, *len);
    vtkIntArray *newScalars = vtkIntArray::New();
    newScalars->SetName( tag.c_str() );
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(ncnt);

    for(unsigned i=0; i<ncnt; i++)
      newScalars->InsertValue(i, vect[i]);
  
    dataSet->GetPointData()->AddArray(newScalars);
    newScalars->Delete();
    return;
  }

  void vtkwritefsn(float *vect, char *name, int *len){
    string tag(name, *len);
    vtkFloatArray *newScalars = vtkFloatArray::New();
    newScalars->SetName( tag.c_str() );
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(ncnt);
  
    for(unsigned i=0; i<ncnt; i++)
      newScalars->InsertValue(i, vect[i]);
  
    dataSet->GetPointData()->AddArray(newScalars);
    newScalars->Delete();
    return;
  }

  void vtkwritedsn(double *vect, char *name, int *len){ 
    string tag(name, *len);
    vtkDoubleArray *newScalars = vtkDoubleArray::New();
    newScalars->SetName( tag.c_str() );
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(ncnt);
  
    for(unsigned i=0; i<ncnt; i++)
      newScalars->InsertValue(i, vect[i]);
  
    dataSet->GetPointData()->AddArray(newScalars);
    newScalars->Delete();
    return;
  }

  void vtkwritefvn(float *vx, float *vy, float *vz,
                      char *name, int *len){
    string tag(name, *len);
    vtkFloatArray *newVectors = vtkFloatArray::New();  
    newVectors->SetName( tag.c_str() );
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(ncnt);

    for(unsigned i=0; i<ncnt; i++){
      newVectors->SetTuple3(i, vx[i], vy[i], vz[i]);
    }

    dataSet->GetPointData()->AddArray(newVectors);
    newVectors->Delete();
    return;
  }

  void vtkwritedvn(double *vx, double *vy, double *vz,
                      char *name, int *len){
    string tag(name, *len);
    vtkDoubleArray *newVectors = vtkDoubleArray::New();  
    newVectors->SetName( tag.c_str() );
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(ncnt);

    for(unsigned i=0; i<ncnt; i++){
      newVectors->SetTuple3(i, vx[i], vy[i], vz[i]);
    }

    dataSet->GetPointData()->AddArray(newVectors);
    newVectors->Delete();
    return;
  }

  void vtkwriteftn(float *v1, float *v2, float *v3, 
                      float *v4, float *v5, float *v6, 
                      float *v7, float *v8, float *v9, 
                      char *name, int *len){
    string tag(name, *len);
    vtkFloatArray *newTensors = vtkFloatArray::New();  
    newTensors->SetName( tag.c_str() );
    newTensors->SetNumberOfComponents(9);
    newTensors->SetNumberOfTuples(ncnt);
  
    for(unsigned i=0; i<ncnt; i++){
      newTensors->SetTuple9(i, 
          v1[i], v2[i], v3[i],
          v4[i], v5[i], v6[i],
          v7[i], v8[i], v9[i]);
    }
  
    dataSet->GetPointData()->AddArray(newTensors);
    newTensors->Delete();
    return;
  }

  void vtkwritedtn(double *v1, double *v2, double *v3, 
                      double *v4, double *v5, double *v6, 
                      double *v7, double *v8, double *v9, 
                      char *name, int *len){
    string tag(name, *len);
    vtkDoubleArray *newTensors = vtkDoubleArray::New();  
    newTensors->SetName( tag.c_str() );
    newTensors->SetNumberOfComponents(9);
    newTensors->SetNumberOfTuples(ncnt);
  
    for(unsigned i=0; i<ncnt; i++){
      newTensors->SetTuple9(i, 
          v1[i], v2[i], v3[i],
          v4[i], v5[i], v6[i],
          v7[i], v8[i], v9[i]);
    }
  
    dataSet->GetPointData()->AddArray(newTensors);
    newTensors->Delete();
    return;
  }

  void vtkstartc(){}

  void vtkwriteisc(int *vect, char *name, int *len){
    string tag(name, *len);
    vtkIntArray *newScalars = vtkIntArray::New();
    newScalars->SetName( tag.c_str() );
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(ecnt);

    for(unsigned i=0; i<ecnt; i++)
      newScalars->InsertValue(i, vect[i]);
  
    dataSet->GetCellData()->AddArray(newScalars);
    newScalars->Delete();
    return;
  }

    void vtkwriteghostlevels(int *ghost_levels) {
        vtkUnsignedCharArray *ghosts = vtkUnsignedCharArray::New();
        ghosts->SetName(vtkDataSetAttributes::GhostArrayName());
        ghosts->SetNumberOfTuples(ecnt);
        
        for(unsigned i=0; i<ecnt; i++) {
            ghosts->SetValue(i, ghost_levels[i] ? vtkDataSetAttributes::HIDDENCELL : 0);
        }
        
        dataSet->GetCellData()->AddArray(ghosts);
        ghosts->Delete();
    }


  void vtkwritefsc(float *vect, char *name, int *len){
    string tag(name, *len);
    vtkFloatArray *newScalars = vtkFloatArray::New();
    newScalars->SetName( tag.c_str() );
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(ecnt);
  
    for(unsigned i=0; i<ecnt; i++)
      newScalars->InsertValue(i, vect[i]);
  
    dataSet->GetCellData()->AddArray(newScalars);
    newScalars->Delete();
    return;
  }

  void vtkwritedsc(double *vect, char *name, int *len){
    string tag(name, *len);
    vtkDoubleArray *newScalars = vtkDoubleArray::New();
    newScalars->SetName( tag.c_str() );
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(ecnt);
  
    for(unsigned i=0; i<ecnt; i++)
      newScalars->InsertValue(i, vect[i]);
  
    dataSet->GetCellData()->AddArray(newScalars);
    newScalars->Delete();
    return;
  }

  void vtkwritefvc(float *vx, float *vy, float *vz,
                      char *name, int *len){
    string tag(name, *len);
    vtkFloatArray *newVectors = vtkFloatArray::New();  
    newVectors->SetName( tag.c_str() );
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(ecnt);

    for(unsigned i=0; i<ecnt; i++){
      newVectors->SetTuple3(i, vx[i], vy[i], vz[i]);
    }

    dataSet->GetCellData()->AddArray(newVectors);
    newVectors->Delete();
    return;
  }

  void vtkwritedvc(double *vx, double *vy, double *vz,
                      char *name, int *len){
    string tag(name, *len);
    vtkDoubleArray *newVectors = vtkDoubleArray::New();  
    newVectors->SetName( tag.c_str() );
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(ecnt);

    for(unsigned i=0; i<ecnt; i++){
      newVectors->SetTuple3(i, vx[i], vy[i], vz[i]);
    }

    dataSet->GetCellData()->AddArray(newVectors);
    newVectors->Delete();
    return;
  }

  void vtkwriteftc(float *v1, float *v2, float *v3, 
                      float *v4, float *v5, float *v6, 
                      float *v7, float *v8, float *v9, 
                      char *name, int *len){
    string tag(name, *len);
    vtkFloatArray *newTensors = vtkFloatArray::New();  
    newTensors->SetName( tag.c_str() );
    newTensors->SetNumberOfComponents(9);
    newTensors->SetNumberOfTuples(ecnt);
  
    for(unsigned i=0; i<ecnt; i++){
      newTensors->SetTuple9(i, 
          v1[i], v2[i], v3[i],
          v4[i], v5[i], v6[i],
          v7[i], v8[i], v9[i]);
    }
  
    dataSet->GetCellData()->AddArray(newTensors);
    newTensors->Delete();
    return;
  }

  void vtkwritedtc(double *v1, double *v2, double *v3, 
                      double *v4, double *v5, double *v6, 
                      double *v7, double *v8, double *v9, 
                      char *name, int *len){
    string tag(name, *len);
    vtkDoubleArray *newTensors = vtkDoubleArray::New();  
    newTensors->SetName( tag.c_str() );
    newTensors->SetNumberOfComponents(9);
    newTensors->SetNumberOfTuples(ecnt);
  
    for(unsigned i=0; i<ecnt; i++){
      newTensors->SetTuple9(i, 
          v1[i], v2[i], v3[i],
          v4[i], v5[i], v6[i],
          v7[i], v8[i], v9[i]);
    }
  
    dataSet->GetCellData()->AddArray(newTensors);
    newTensors->Delete();
    return;
  }

  void vtkclose(){
    vtkXMLUnstructuredGridWriter *writer= vtkXMLUnstructuredGridWriter::New();
    vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
    
    writer->SetFileName( fl_vtkFileName.c_str() );
    // writer->SetDataModeToAppended();
    writer->SetDataModeToBinary();
    writer->SetEncodeAppendedData(1);       // Base64-encode appended payload
    writer->SetCompressorTypeToNone();  
    writer->EncodeAppendedDataOff();
    writer->SetInputData(dataSet);
    writer->SetCompressor(compressor);
    compressor->Delete();
  
    writer->Write();
    writer->Delete();
    dataSet->Delete();
    dataSet = NULL;
    return;
  }

  void _vtkpclose_nointerleave(const int *rank, const int *npartitions){
    if((*npartitions)<2){
      vtkclose();
      return;
    }
    
    vtkXMLPUnstructuredGridWriter *writer= vtkXMLPUnstructuredGridWriter::New();
    vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
    
    string filename = fl_vtkFileName;
    bool is_pvtu = false;
    string basename;
    if(fl_vtkFileName.size()>4){
      is_pvtu = string(fl_vtkFileName, fl_vtkFileName.size()-4, 4)=="pvtu";
      
      basename = string(fl_vtkFileName, 0, fl_vtkFileName.size()-5)+"/";
      mkdir(basename.c_str(), 0777);
      filename = basename+filename;
    }

    writer->SetFileName(filename.c_str());
    writer->SetNumberOfPieces(*npartitions);
    writer->SetGhostLevel(1);
    writer->SetStartPiece(*rank);
    writer->SetEndPiece(*rank);
    writer->SetInputData(dataSet);
    writer->SetCompressor(compressor);
    compressor->Delete();
    // writer->SetDataModeToAppended();
    writer->SetDataModeToBinary();
    writer->SetEncodeAppendedData(1);       // Base64-encode appended payload
    writer->SetCompressorTypeToNone();  
    writer->EncodeAppendedDataOff();

#ifdef HAVE_MPI
    if (!writer->GetController()) {
      vtkMPIController *cont = vtkMPIController::New();
      cont->Initialize();
      writer->SetController(cont);
      cont->Delete();
    }
    writer->SetWriteSummaryFile(1);
#else
    writer->SetWriteSummaryFile((*rank)==0);
#endif
    
    writer->Write();
    writer->Delete();
    dataSet->Delete();
    dataSet = NULL;
  
    if(is_pvtu && (*rank)==0){
      rename(filename.c_str(), fl_vtkFileName.c_str());
      pvtu_fix_path(fl_vtkFileName.c_str(), basename.c_str());
    }
    return;
  }

  void vtkpclose(int *rank, int *npartitions){
#ifdef HAVE_MPI
#define INTERLEAVE_IO_TRESHOLD 64000
    if(*npartitions>INTERLEAVE_IO_TRESHOLD){
      int nwrites = (int)(sqrt(*npartitions)+0.5);
      
      for(int lrank=0; lrank<nwrites; lrank++){
        if((*rank)%nwrites==lrank){
          _vtkpclose_nointerleave(rank, npartitions);
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }else{
#endif
      _vtkpclose_nointerleave(rank, npartitions);
#ifdef HAVE_MPI
    }
#endif
    return;
  }
  
  void vtksetactivescalars(char* name, int *len){
    string tag(name, *len);
    dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::SCALARS);
    return;
  }

  void vtksetactivevectors(char* name, int *len){
    string tag(name, *len);
    dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::VECTORS);
    return;
  }

  void vtksetactivetensors(char* name, int *len){
    string tag(name, *len);
    dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::TENSORS);
    return;
  }
}
 
 
