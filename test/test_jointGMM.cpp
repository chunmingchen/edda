#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "vector"
#include "bmp_image.h"

using namespace edda;
using namespace std;
using namespace edda::dist;

//using image to test joint Gaussian GMM
int main()
{
	//load testing image from disk
	BMPImage image("test_img.bmp");

	//down-sampled block size
	int blockSize = 20;
	int nVar = 3;//number of variables, rgb is 3 variables
	int nGmmComp = 2;//number of Gaussian compoenents
	

	//output image means and resampled image
	BMPImage optImage("test_img.bmp");//output image, randomly give a image file

	//number of row and col after down-sampleing
	int dsVs = image.height / blockSize;
	int dsUs = image.width / blockSize;
	
	//Joint GMM array
	shared_ary<JointGMM> array(new JointGMM[dsVs*dsUs], dsVs*dsUs);
	thrust::default_random_engine rng;//random engine for getJointSample()

	//loop: go through each block
	for (int dsV = 0; dsV < dsVs; dsV++){
		for (int dsU = 0; dsU < dsUs; dsU++){
			printf("%d %d\n", dsV, dsU);
			//prepare training vectors to OpenCV EM
			ublas_matrix trainSamples(blockSize* blockSize, nVar);
			int cnt = 0;
			for (int v = dsV*blockSize; v < (dsV + 1)*blockSize; v++) {
				for (int u = dsU*blockSize; u < (dsU + 1)*blockSize; u++) {
					for (int k = 0; k < nVar; k++){
						trainSamples(cnt, k) = (Real)(image.bitmapImage[(v*image.width + u) * nVar + k]);
					}
					cnt++;
				}
			}

			//EM in Edda
			array[dsV*dsUs + dsU] = eddaComputeJointGMM(nGmmComp, trainSamples);

			//resample this block and write to image		
			for (int i = 0; i < blockSize; i++){
				for (int j = 0; j < blockSize; j++){
					std::vector<Real> sample = getJointSample(array[dsV*dsUs + dsU], rng);

					int rawV = dsV*blockSize + i;//u and v in the original image resolution
					int rawU = dsU*blockSize + j;

					//trim value to 0-255, it possible to sample any value from gaussian
					//and store back to image buffer
					for (int k = 0; k < nVar; k++){
						if (sample[k] < 0)  sample[k] = 0;
						if (sample[k] > 255)sample[k] = 255;

						optImage.bitmapImage[(rawV*image.width + rawU) * 3 + k] = sample[k];
					}
				}
			}
		}
	}

	//write three images to disk
	optImage.writeImage(std::string("jointGMMTestOutput.bmp").c_str());//this is the output file name

	std::cout << "Press any key to finish" << std::endl;
	getchar();
	
}