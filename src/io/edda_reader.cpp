#include <memory>
#include <vtkDataSet.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkFieldData.h>
#include <vtkStringArray.h>
#include <vtkStdString.h>
#include <vtkPoints.h>

#include "io/path.h"

#include "io/gmm_vtk_data_array.h"
#include "io/edda_reader.h"
#include "dataset/distr_array.h"
#include "distributions/histogram.h"

using namespace std;
using namespace edda;
using namespace dist;

namespace edda {

	DistrArray * readGMMArray(ifstream & myfile, int n)
	{
		int nc;
		myfile.read((char*)(&nc), sizeof (int));
		int isNumGaussianUniform;
		myfile.read((char*)(&isNumGaussianUniform), sizeof (int));
		if (isNumGaussianUniform == 1){
			int GMs;
			myfile.read((char*)(&GMs), sizeof (int));
			float* gmData = new float[GMs * 3 * nc*n];
			myfile.read((char*)(gmData), sizeof (float)*GMs * 3 * nc*n);

			if (nc == 1){
				shared_ary<GaussianMixture<5>> s_gmAry(new GaussianMixture<5>[n], n);
				//TODO:: based on GMs may use other than GMM5

				for (int c = 0; c < nc; c++){
					for (int j = 0; j < n; j++){
						std::vector<GMMTuple> models;
						for (int i = 0; i < GMs * 3; i += 3){
							int index = c*n*GMs * 3 + j*GMs * 3 + i;
							float p[3] = { gmData[index], gmData[index + 1], gmData[index+2] };
							GMMTuple curG;
							curG.m = p[0];
							curG.v = p[1];
							curG.w = p[2];
							curG.p[0] = p[0];
							curG.p[1] = p[1];
							curG.p[2] = p[2];
							models.push_back(curG);
						}
						GaussianMixture<5> curGM(models);
						s_gmAry[j] = curGM;
					}
				}
				return (DistrArray *)(new ScalarDistrArray<GaussianMixture<5>>(s_gmAry));
			}
			else if (nc == 3){
				shared_ary<Vector<GaussianMixture<5>, 3> > s_gmAry(new Vector<GaussianMixture<5>, 3>[n], n);
				//TODO:: based on GMs may use other than GMM5
				//TODO:: saving order of vector data?
				for (int j = 0; j < n; j++){
					Vector<GaussianMixture<5>, 3> curVectorDistr;

					for (int c = 0; c < nc; c++){
						std::vector<GMMTuple> models;
						for (int i = 0; i < GMs * 3; i += 3){
							int index = c*n*GMs * 3 + j*GMs * 3 + i;
							float p[3] = { gmData[index], gmData[index + 1], gmData[index + 2] };
							GMMTuple curG;
							curG.m = p[0];
							curG.v = p[1];
							curG.w = p[2];
							curG.p[0] = p[0];
							curG.p[1] = p[1];
							curG.p[2] = p[2];
							models.push_back(curG);
						}
						GaussianMixture<5> curGM(models);
						curVectorDistr[c] = curGM;				
					}

					s_gmAry[j] = curVectorDistr;
				}

				return (DistrArray *)(new VectorDistrArray<GaussianMixture<5>, 3>(s_gmAry));
			}
			else{
				delete gmData;
				myfile.close();
				throw NotImplementedException();
			}
		}
		else{
			myfile.close();
			throw NotImplementedException();
		}
	}


	DistrArray *readHistoArray(ifstream & myfile, int n)//vtkPointData *vtk_point_data)
	{
		int nc;
		myfile.read((char*)(&nc), sizeof (int));
		
		int isNumBinsUniform;
		myfile.read((char*)(&isNumBinsUniform), sizeof (int));
		int nbins;
		if (isNumBinsUniform == 1){
			myfile.read((char*)(&nbins), sizeof (int));
		}

		int isMinMaxValueUniform;
		myfile.read((char*)(&isMinMaxValueUniform), sizeof (int));
		float minv, maxv;
		if (isMinMaxValueUniform == 1){
			myfile.read((char*)(&minv), sizeof (float));
			myfile.read((char*)(&maxv), sizeof (float));
		}

		if (nc == 1){
			if (isNumBinsUniform == 1 && isMinMaxValueUniform == 1){
				//not tested yet
				float* histData = new float[nbins + 2];
				histData[0] = minv;
				histData[1] = maxv;

				shared_ary<Histogram> histAry(new Histogram[n], n);
				for (int nn = 0; nn < n; nn++){				
					myfile.read((char*)(histData + 2), sizeof (float)*nbins);
					histAry[nn] = Histogram(histData, nbins);
				}

				free(histData);

				return new ScalarDistrArray<Histogram>(histAry);
			}
			else if(isNumBinsUniform == 1 && isMinMaxValueUniform != 1){
				myfile.close();
				throw NotImplementedException();
			}
			else if (isNumBinsUniform != 1 && isMinMaxValueUniform == 1){
				myfile.close();
				throw NotImplementedException();
			}
			else{ //(isNumBinsUniform == 0 && isMinMaxValueUniform == 0)

				int* headerBinary_nbins = (int*)malloc(sizeof(int)*n);
				float* headerBinary_minMaxV = (float*)malloc(sizeof(float)* 2 * n);
				myfile.read((char*)(headerBinary_nbins), sizeof (int)*n);
				myfile.read((char*)(headerBinary_minMaxV), sizeof (float)* 2 * n);

				shared_ary<Histogram> histAry(new Histogram[n], n);

				for (int j = 0; j < n; j++){
					int n_bins = headerBinary_nbins[j];

					float* histData = new float[n_bins + 2];
					histData[0] = headerBinary_minMaxV[2*j];
					histData[1] = headerBinary_minMaxV[2*j+1];

					myfile.read((char*)(histData+2), sizeof (float)*n_bins);
					histAry[j] = Histogram(histData, n_bins);

					free(histData);
				}

				free(headerBinary_nbins);
				free(headerBinary_minMaxV);

				return new ScalarDistrArray<Histogram>(histAry);
			}
		}
		else{
			myfile.close();
			throw NotImplementedException();
		}
	}

	DistrArray * readMixArray(ifstream & myfile, int n)
	{
		shared_ary<dist::Variant> pArray(new dist::Variant[n], n);

		//int nc;
		//myfile.read((char*)(&nc), sizeof (int));
		//if (nc == 1){
			float* gmData = new float[333 * 3]; // !!! this array is designed to avoid the time cost by multiple memory allocation and deletion. better check if the size is big enough before each time using it !!!

			for (int j = 0; j < n; j++){
				int distrTypeNumber = 1;
				myfile.read((char*)(&distrTypeNumber), sizeof (int));
				if (distrTypeNumber == 1){ //gaussian
					
					int nGM;
					myfile.read((char*)(&nGM), sizeof (int));

					myfile.read((char*)(gmData), sizeof (float)*nGM * 3);

					dist::GMM new_gmm = dist::GMM(nGM);
					for (int m = 0; m<nGM; m++)
					{
						new_gmm.models[m].m = gmData[3 * m];
						new_gmm.models[m].v = gmData[3 * m + 1];
						new_gmm.models[m].w = gmData[3 * m + 2];
					}
					pArray[j] = new_gmm;
				}
				else if (distrTypeNumber == 2){ //histogram
					int nbins;
					myfile.read((char*)(&nbins), sizeof (int));

					myfile.read((char*)(gmData), sizeof (float));
					myfile.read((char*)(gmData+1), sizeof (float));			
					myfile.read((char*)(gmData+2), sizeof (float)*nbins);

					pArray[j] = Histogram(gmData, nbins);
				}
				else if (distrTypeNumber == 3){ //joint gmm
					int nVar, nComp;
					myfile.read((char*)(&nVar), sizeof (int));
					myfile.read((char*)(nComp), sizeof (int));

					//... to be done
				}
				else{
					myfile.close();
					throw NotImplementedException();
				}
			}

			delete[] gmData;

			return (DistrArray *)(new ScalarDistrArray<dist::Variant>(pArray));		
		//}
		//else{
		//	myfile.close();
		//	throw NotImplementedException();
		//}
	}
	


	template <typename T>
	shared_ptr<Dataset<T> > loadEddaDatasetTemplate(const string &edda_file)
	{
		ifstream myfile(edda_file.c_str(), ios::binary);

		streampos begin = myfile.tellg();
		myfile.seekg(0, ios::end);
		streampos end = myfile.tellg();
		int fileByteSize = end - begin;

		myfile.seekg(0, ios::beg);

		char eddaFileMark[4];
		myfile.read(eddaFileMark, sizeof(char)*4);
		if (eddaFileMark[0] != 'E' || eddaFileMark[1] != 'D' || eddaFileMark[2] != 'D' || eddaFileMark[3] != 'A'){
			myfile.close();
			//TODO: throw a proper exception
			exit(0);
		}

		char majorVersion, minorVersion;
		myfile.read(&majorVersion, sizeof(char));
		myfile.read(&minorVersion, sizeof(char));

		if (majorVersion == 0 && minorVersion == 1){
			int gridTypeNumber;
			myfile.read((char*)(&gridTypeNumber), sizeof(int));
			if (gridTypeNumber == 1){
				int dims[3];
				myfile.read((char*)(&dims), sizeof (int)* 3);
				float spacing[3];
				myfile.read((char*)(&spacing), sizeof (float)* 3);

				int distrTypeNumber = 1;
				myfile.read((char*)(&distrTypeNumber), sizeof (int));
				if (distrTypeNumber == 1){
					shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
						new RegularCartesianGrid(dims[0], dims[1], dims[2]),
						readGMMArray(myfile, dims[0]*dims[1]*dims[2]));
					return dataset;
				}
				else if (distrTypeNumber == 2){
					shared_ptr<Dataset<T> > dataset = make_Dataset<T>(
						new RegularCartesianGrid(dims[0], dims[1], dims[2]),
						readHistoArray(myfile, dims[0] * dims[1] * dims[2])
						);
					return dataset;
				}
				else{
					myfile.close();
					throw NotImplementedException();
				}
			}
			else{
				myfile.close();
				throw NotImplementedException();
			}
		}
		else if (majorVersion == 0 && minorVersion == 2){
			int gridTypeNumber;
			myfile.read((char*)(&gridTypeNumber), sizeof(int));
			if (gridTypeNumber == 1){
				int dims[3];
				myfile.read((char*)(&dims), sizeof (int)* 3);
				float spacing[3];
				myfile.read((char*)(&spacing), sizeof (float)* 3);

				
				int numDistrArray = 1;
				myfile.read((char*)(&numDistrArray), sizeof (int));
				
				std::vector<DistrArray *> dVec(numDistrArray);
				for (int i = 0; i < numDistrArray; i++){
					dVec[i] = readMixArray(myfile, dims[0] * dims[1] * dims[2]);
				}
				
				myfile.close();

				return make_shared<Dataset<T>>(new RegularCartesianGrid(dims[0], dims[1], dims[2]), dVec);
			}
			else{
				myfile.close();
				throw NotImplementedException();
			}
		}
		else{
			myfile.close();
			throw NotImplementedException();
		}

		printf("Read file from %s.\n", edda_file.c_str());

	}



	shared_ptr<Dataset<Real> > loadEddaScalarDataset_noneVTK(const string &edda_file)
	{
		return loadEddaDatasetTemplate<Real>(edda_file);
	}

	shared_ptr<Dataset<VECTOR3> > loadEddaVector3Dataset_noneVTK(const string &edda_file)
	{
		return loadEddaDatasetTemplate<VECTOR3>(edda_file);
	}
}