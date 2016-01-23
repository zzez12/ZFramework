#include "ZMeshFilterManifold.h"
#include "ZMeshFilterGaussian.h"
#include "ZMeshFilterGaussian2.h"
#include "ZMeshFilterWeightedGaussian.h"
#include "ZMeshBilateralFilter.h"
#include "MeshTools.h"
#include "ZKernelFunc.h"

#include "ZFileHelper.h"
#include <fstream>
#include <sstream>

namespace ZMeshSpace
{
	ZMeshFilterManifold::ZMeshFilterManifold(Mesh3D* mesh)
		: ZMeshFilter(mesh, MESH_FILTER_MANIFOLD)
	{
		pGaussianFilter_ = NULL;
		pAnnSearch_ = NULL;
	}

	ZMeshFilterManifold::~ZMeshFilterManifold()
	{
		destroy();
	}

	void ZMeshFilterManifold::destroy()
	{
		SAFE_DELETE(pGaussianFilter_);
		SAFE_DELETE(pAnnSearch_);
		sum_w_ki_Psi_blur_.resize(0,0);
		sum_w_ki_Psi_blur_0_.resize(0);
		min_pixel_dist_to_manifold_squared_.resize(0);
		wki_Psi_blurs_.clear();
		wki_Psi_blur_0s_.clear();
	}

	void ZMeshFilterManifold::init(const Eigen::MatrixXf& input)
	{
		destroy();
		computeTreeHeight();
		initMinPixelDist2Manifold();
		// init search tool
		pAnnSearch_ = new AnnWarper_Eigen;
		//Eigen::MatrixXf faceCenters = MeshTools::getAllFaceCenters(pMesh_);
		//Eigen::MatrixXf faceNormals = MeshTools::getAllFaceNormals(pMesh_);
		pAnnSearch_->init(input.block(0, 0, input.rows(), spatialDim_));
		// init gaussian filter
		pGaussianFilter_ = new ZMeshBilateralFilter(pMesh_);
		ZMeshBilateralFilter* pFilter = (ZMeshBilateralFilter*)pGaussianFilter_;
		pFilter->setAnnSearchHandle(pAnnSearch_);
	}

	void ZMeshFilterManifold::setPara(const std::string& name, float value)
	{
		if (name.compare(ZMeshFilterParaNames::SpatialSigma)==0)
		{
			filterPara_.spatial_sigma = value;
		}
		else if (name.compare(ZMeshFilterParaNames::RangeSigma)==0)
		{
			filterPara_.range_sigma = value;
		}
	}

	void ZMeshFilterManifold::computeTreeHeight()
	{
		int Hs = floor(log(filterPara_.spatial_sigma))-1;
		float Lr = 1-filterPara_.range_sigma;
		int height = ceil(Hs*Lr);
	//	filterPara_.tree_height = height > 2 ? height : 2; 
		filterPara_.tree_height = 3; 
	}

	void ZMeshFilterManifold::initMinPixelDist2Manifold()
	{
		if (pMesh_==NULL) return;
		int vSize = pMesh_->get_num_of_vertex_list();
		min_pixel_dist_to_manifold_squared_.resize(vSize);
		min_pixel_dist_to_manifold_squared_.fill(FLT_MAX);
	}

	bool ZMeshFilterManifold::apply()
	{

		return true;
	}

	bool ZMeshFilterManifold::apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags)
	{
		init(input);

		int inputSize = input.rows();//pMesh_->get_num_of_vertex_list();
		Eigen::MatrixXf initAttributes = input;

		sum_w_ki_Psi_blur_.resize(inputSize, rangeDim_);
		sum_w_ki_Psi_blur_0_.resize(inputSize);
		min_pixel_dist_to_manifold_squared_.resize(inputSize);	// un-initilized
		min_pixel_dist_to_manifold_squared_.fill(FLT_MAX);
		sum_w_ki_Psi_blur_.fill(0);
		sum_w_ki_Psi_blur_0_.fill(0);
		verticeClusterIds.resize(inputSize);
		verticeClusterIds.fill(1);

		// 1. low-pass filtering
		//ZMeshFilterGaussian filterGaussian(pMesh_);
		float spatialSigma = filterPara_.spatial_sigma * 1.f;
		float rangeSigma = cos(filterPara_.range_sigma*Z_PI/180.f);
		pGaussianFilter_->setPara(ZMeshFilterParaNames::SpatialSigma, spatialSigma);
		pGaussianFilter_->setPara(ZMeshFilterParaNames::RangeSigma, rangeSigma);
		pGaussianFilter_->setKernelFunc(NULL); // only using spatial filter
		//pGaussianFilter_->apply(initPositions, std::vector<bool>(inputSize, true));
		CHECK_FALSE_AND_RETURN(pGaussianFilter_->apply(initAttributes, spatialDim_, rangeDim_));
		Eigen::MatrixXf eta1 = pGaussianFilter_->getResult();
		//std::cout << " input=\n" << input.block(0,0,10,rangeDim_+spatialDim_) << "\n";
		//std::cout << " eta1=\n" << eta1.block(0,0,10,rangeDim_+spatialDim_) << "\n";
		//ZFileHelper::saveEigenMatrixToFile("eta1.txt", eta1);
		std::vector<bool> cluster1(inputSize, true);

		int currentTreeLevel = 1;
		CHECK_FALSE_AND_RETURN(buildManifoldAndPerformFiltering(input, eta1, cluster1, spatialSigma, rangeSigma, currentTreeLevel));
		//std::cout << "sum_w_ki_Psi_blur=\n" << sum_w_ki_Psi_blur_.block(0,0,10,rangeDim_) << "\n";
		//std::cout << "Sum_w_ki_Psi_blur norm():" << sum_w_ki_Psi_blur_.norm() << "\n";
		//std::cout << "sum_w_ki_Psi_blur_0=\n" << sum_w_ki_Psi_blur_0_.head(10) << "\n";

		//std::vector<Eigen::Vector3f> tilde_g(vSize);
		Eigen::MatrixXf tilde_g(inputSize, rangeDim_);
		for (int i=0; i<inputSize; i++)
		{
			if (g_isZero(sum_w_ki_Psi_blur_0_(i)))
				std::cout << "Zero! " << i << "\n";
			for (int j=0; j<rangeDim_; j++)
			{
				tilde_g(i, j) = sum_w_ki_Psi_blur_(i, j)/sum_w_ki_Psi_blur_0_(i);
			}
		}
		//ZFileHelper::saveEigenMatrixToFile("SumWeight.txt", sum_w_ki_Psi_blur_0_);

		Eigen::MatrixXf sumNewRange(inputSize, rangeDim_);
		sumNewRange.fill(0);
		rangeWeights_.resize(wki_Psi_blurs_.size(), 1.f);
		for (int i=0; i<wki_Psi_blurs_.size(); i++)
		{
			//float weight = 1.f*(i*10+1);
			float weight = 1.f/(i+1);
			//float weight = 1.f - log(1.f*(i+1));
			//float weight = exp(1.f*(i+1));
			//float weight = i==0 ? 1.f : 0.f;
			sumNewRange += wki_Psi_blurs_[i]*weight;
		}
// 		for (int i=0; i<inputSize; i++)
// 		{
// 			//sumNewRange.row(i) /= sum_w_ki_Psi_blur_0_(i);
// 			sumNewRange.row(i).normalize();
// 		}
		//for (int i=0; i<wki_Psi_blurs_.size(); i++)
		//{
		//	std::stringstream ss;
		//	ss << "wkiPsiBlurs" << i << ".txt";
		//	//ZFileHelper::saveEigenMatrixToFile(ss.str().c_str(), wki_Psi_blurs_[i]);
		//	std::stringstream ss2;
		//	ss2 << "wkiPsiBlurs0" << i << ".txt";
		//	ZFileHelper::saveEigenVectorToFile(ss2.str().c_str(), wki_Psi_blur_0s_[i]);
		//}		
		//ZFileHelper::saveEigenMatrixToFile("sumNewRange.txt", sumNewRange);

		output_ = input;
//		output_.block(0, spatialDim_, inputSize, rangeDim_) = tilde_g;
//		std::cout << output_;
// 		ZFileHelper::saveEigenMatrixToFile("input.txt", input);
// 		ZFileHelper::saveEigenMatrixToFile("minPixelDist2Mani.txt", min_pixel_dist_to_manifold_squared_);
// 		ZFileHelper::saveEigenMatrixToFile("tildeG.txt", tilde_g);
//    		Eigen::VectorXf alpha(inputSize);
//    		alpha.fill(0);
//    		for (int i=0; i<inputSize; i++)
//    		{
//    			alpha[i] = exp(min_pixel_dist_to_manifold_squared_(i)*(-0.5)/filterPara_.range_sigma/filterPara_.range_sigma);
//    		}
//    		for (int i=0; i<inputSize; i++)
//    		{
//   			Eigen::VectorXf newRange = input.row(i).tail(rangeDim_) + (tilde_g.row(i)-input.row(i).tail(rangeDim_))*alpha(i);
//  			newRange.normalize();
//    			output_.row(i).tail(rangeDim_) = newRange;
//   		}
		output_.block(0, spatialDim_, inputSize, rangeDim_) = sumNewRange;

// 		ZFileHelper::saveEigenMatrixToFile("sumWkiPsiBlur.txt", sum_w_ki_Psi_blur_);
// 		ZFileHelper::saveEigenMatrixToFile("tilderG.txt", tilde_g);
// 		ZFileHelper::saveEigenMatrixToFile("output.txt", output_);
		std::cout << "Filter done!\n";
		//std::cout << "output=\n" << output_.block(0,0,10,rangeDim_+spatialDim_) << "\n";
		//std::cout << " Output-input=\n" << (output_-input).block(0,0,10,rangeDim_+spatialDim_) << "\n";

		std::cout << ZConsoleTools::red << " Error: " << (output_-input).norm() << "\n" 
				  << ZConsoleTools::white;

		return true;
	}

	bool ZMeshFilterManifold::buildManifoldAndPerformFiltering(const Eigen::MatrixXf& input, 
		const Eigen::MatrixXf& etaK, const std::vector<bool>& clusterK, 
		float sigma_s, float sigma_r, int currentTreeLevel)
	{
		int inputSize = input.rows();;

		static std::string fileName("etaK.txt");
		fileName += "0";
		ZFileHelper::saveEigenMatrixToFile(fileName.c_str(), etaK);

		// splatting: project the vertices onto the current manifold etaK

		// NOTE: the sigma should be recursive!!
		float sigmaR_over_sqrt_2 = sigma_r/sqrt(2.0);
		Eigen::MatrixXf diffX = input.block(0, spatialDim_, inputSize, rangeDim_) - etaK.block(0, spatialDim_, inputSize, rangeDim_);
		Eigen::VectorXf gaussianDistWeight(inputSize);
		//std::cout << "diffX=\n" << diffX.block(0,0,10,rangeDim_+spatialDim_) << "\n";
		for (int i=0; i<inputSize; i++)
		{
			gaussianDistWeight(i) = ZKernelFuncs::GaussianKernelFunc(diffX.row(i), sigmaR_over_sqrt_2);
		}

		Eigen::MatrixXf psiSplat(inputSize, spatialDim_+rangeDim_+1);
		Eigen::VectorXf psiSplat0 = gaussianDistWeight;
		//////////////////////////////////////////////////////////////////////////
		// for debug
		//std::stringstream ss;
		static int etaI = 1;
		//ss << "gaussianWeights" << etaI << ".txt";
		//std::cout << ZConsoleTools::green << etaI << "\n" << ZConsoleTools::white;
		etaI++;
		//ZFileHelper::saveEigenVectorToFile(ss.str().c_str(), gaussianDistWeight);
		//////////////////////////////////////////////////////////////////////////
		psiSplat.block(0,0,inputSize,spatialDim_) = input.block(0,0,inputSize,spatialDim_);
		for (int i=0; i<inputSize; i++)
		{
			//psiSplat.block(i, spatialDim_, 1, rangeDim_) = input.block(i, spatialDim_, 1, rangeDim_)*gaussianDistWeight(i);
			psiSplat.block(i, spatialDim_, 1, rangeDim_) = etaK.block(i, spatialDim_, 1, rangeDim_)*gaussianDistWeight(i);
		}
		psiSplat.col(spatialDim_+rangeDim_) = gaussianDistWeight;

		// save min distance to later perform adjustment of outliers -- Eq.(10)
		// UNDO
		// update min_pixel_dist_to_manifold_square_
		Eigen::VectorXf pixelDist2Manifold(inputSize);
		for (int i=0; i<inputSize; i++)
		{
			if (clusterK[i])
				pixelDist2Manifold(i) = diffX.row(i).squaredNorm();
			else
				pixelDist2Manifold(i) = FLT_MAX;
		}
		min_pixel_dist_to_manifold_squared_ = min_pixel_dist_to_manifold_squared_.cwiseMin(pixelDist2Manifold);

		// Blurring
		Eigen::MatrixXf wkiPsiBlur(inputSize, rangeDim_);
		Eigen::VectorXf wkiPsiBlur0(inputSize);
		//Eigen::MatrixXf psiSplatAnd0(inputSize, spatialDim_+rangeDim_+1);
		//psiSplatAnd0.block(0, 0, inputSize, spatialDim_+rangeDim_) = psiSplat.block(0, 0, inputSize, spatialDim_+rangeDim_);
		//psiSplatAnd0.col(spatialDim_+rangeDim_) = psiSplat0;
		//psiSplatAnd0.block(0, 0, inputSize, spatialDim_+rangeDim_) = input;
		//psiSplatAnd0.col(spatialDim_+rangeDim_).fill(1.f);
		//pGaussianFilter_->setKernelFunc(ZKernelFuncs::GaussianKernelFuncNormal3);
		pGaussianFilter_->setKernelFunc(ZKernelFuncs::GaussianKernelFuncA);
		//pGaussianFilter_->setKernelFunc(NULL);
		//pGaussianFilter_->setKernelFunc(ZKernelFuncs::GaussianKernelSphericalFunc);
		pGaussianFilter_->setPara(ZMeshFilterParaNames::SpatialSigma, sigma_s);
		pGaussianFilter_->setPara(ZMeshFilterParaNames::RangeSigma, sigmaR_over_sqrt_2);
		//pGaussianFilter_->apply(psiSplatAnd0, clusterK);
		//CHECK_FALSE_AND_RETURN(pGaussianFilter_->apply(psiSplatAnd0, spatialDim_, rangeDim_+1, Eigen::VectorXf(), clusterK));
		CHECK_FALSE_AND_RETURN(pGaussianFilter_->apply(psiSplat, spatialDim_, rangeDim_+1, gaussianDistWeight, clusterK));
		Eigen::MatrixXf wkiPsiBlurAnd0 = pGaussianFilter_->getResult();
		wkiPsiBlur = wkiPsiBlurAnd0.block(0, spatialDim_, inputSize, rangeDim_);
		wkiPsiBlur0 = wkiPsiBlurAnd0.col(spatialDim_+rangeDim_);
		//ZFileHelper::saveEigenMatrixToFile("")
// 		std::cout << "wkiPsiBlurAnd0=\n" << wkiPsiBlurAnd0.block(0,0,10,rangeDim_+spatialDim_+1) << "\n";
// 		std::cout << "wkiPsiBlur0=\n" << wkiPsiBlur0.head(10) << "\n";

		//std::cout << "pixelDist2Manifold=\n" << pixelDist2Manifold.head(10) << "\n";
		//std::cout << "min_pixel_dist_to_manifold_squared=\n" << min_pixel_dist_to_manifold_squared_.head(10) << "\n";

		// for debug
// 		for (int i=0; i<inputSize; i++)
// 		{
// 			if (!g_isInfinite(wkiPsiBlur(i,0)))
// 				std::cout << "(" << i << "," << 0 << ")\n";
// 			if (!g_isInfinite(wkiPsiBlur(i,1)))
// 				std::cout << "(" << i << "," << 1 << ")\n";
// 			//if (!g_isInfinite(wkiPsiBlur(i,2)))
// 			//	std::cout << "(" << i << "," << 2 << ")\n";
// 		}
		//std::cout << wkiPsiBlur.norm() << "\n";

		Eigen::VectorXf rangeDiff(inputSize);
		for (int i=0; i<inputSize; i++)
		{
			Eigen::VectorXf n0 = wkiPsiBlur.row(i);
			Eigen::VectorXf n1 = input.row(i).tail(rangeDim_);
			n0.normalize();
			n1.normalize();
			rangeDiff(i) = 1.f-n0.dot(n1);
		}
		static bool bSaved = false;
		if (!bSaved)
		{
			ZFileHelper::saveEigenVectorToFile("rangeDiff.txt", rangeDiff);
			ZFileHelper::saveEigenVectorToFile("gaussian.txt", gaussianDistWeight);
			ZFileHelper::saveEigenMatrixToFile("splat.txt", psiSplat);
			ZFileHelper::saveEigenMatrixToFile("blur.txt", wkiPsiBlur);
			bSaved = true;
		}

		// Slicing
		Eigen::VectorXf wki = gaussianDistWeight;
		for (int i=0; i<inputSize; i++)
		{
			if (!clusterK[i]) continue;
			sum_w_ki_Psi_blur_.row(i) += wkiPsiBlur.row(i)*wki(i);
			sum_w_ki_Psi_blur_0_(i) += wkiPsiBlur0(i)*wki(i);
		}
		//////////////////////////////////////////////////////////////////////////
		// for debug
		wki_Psi_blurs_.push_back(Eigen::MatrixXf(inputSize, rangeDim_));
		wki_Psi_blur_0s_.push_back(Eigen::VectorXf(inputSize));
		Eigen::MatrixXf& lastM = wki_Psi_blurs_[wki_Psi_blurs_.size()-1];
		lastM.fill(0);
		Eigen::VectorXf& lastV = wki_Psi_blur_0s_[wki_Psi_blur_0s_.size()-1];
		lastV.fill(0);
		for (int i=0; i<inputSize; i++)
		{
			if (!clusterK[i]) continue;
			lastM.row(i) = wkiPsiBlur.row(i)*wki(i);
			lastV(i) = wkiPsiBlur0(i)*wki(i);
		}
		std::cout << sum_w_ki_Psi_blur_.norm() << "\n";
		//////////////////////////////////////////////////////////////////////////

		// compute two new manifolds eta_minus and eta_plus

		// test stopping criterion
		if (currentTreeLevel<filterPara_.tree_height)
		{
			// algorithm 1, Step 2: compute the eigenvector v1
			Eigen::VectorXf v1 = computeMaxEigenVector(diffX, clusterK);

			// algorithm 1, Step 3: Segment vertices into two clusters
			std::vector<bool> clusterMinus(inputSize, false);
			std::vector<bool> clusterPlus(inputSize, false);
			int countMinus=0;
			int countPlus =0;
			for (int i=0; i<inputSize; i++)
			{
				float dot = diffX.row(i).dot(v1);
				if (dot<0 && clusterK[i]) {countMinus++; verticeClusterIds[i] =etaI+0.5; clusterMinus[i] = true;}
				if (dot>=0 && clusterK[i]) {countPlus++; verticeClusterIds[i] =etaI-0.5; clusterPlus[i] = true;}
			}
			std::cout << "Minus manifold: " << countMinus << "\n";
			std::cout << "Plus manifold: " << countPlus << "\n"; 
// 			Eigen::MatrixXf diffXManifold(inputSize, spatialDim_+rangeDim_);
// 			diffXManifold.block(0, 0, inputSize, spatialDim_) = input.block(0, 0, inputSize, spatialDim_);
// 			diffXManifold.block(0, spatialDim_, inputSize, rangeDim_) = diffX;

			// algorithm 1, Step 4: Compute new manifolds by weighted low-pass filtering  -- Eq. (7)(8)
			Eigen::VectorXf theta(inputSize);
			theta.fill(1);
			theta = theta - wki.cwiseProduct(wki);
			pGaussianFilter_->setKernelFunc(NULL);
			CHECK_FALSE_AND_RETURN(pGaussianFilter_->apply(input, spatialDim_, rangeDim_, theta, clusterMinus));
			Eigen::MatrixXf etaMinus = pGaussianFilter_->getResult();
			CHECK_FALSE_AND_RETURN(pGaussianFilter_->apply(input, spatialDim_, rangeDim_, theta, clusterPlus));
			Eigen::MatrixXf etaPlus = pGaussianFilter_->getResult();

			// algorithm 1, Step 5: recursively build more manifolds
			CHECK_FALSE_AND_RETURN(buildManifoldAndPerformFiltering(input, etaMinus, clusterMinus, sigma_s, sigma_r, currentTreeLevel+1));
			CHECK_FALSE_AND_RETURN(buildManifoldAndPerformFiltering(input, etaPlus, clusterPlus, sigma_s, sigma_r, currentTreeLevel+1));
		}

		return true;
	}

	Eigen::VectorXf ZMeshFilterManifold::computeMaxEigenVector(const Eigen::MatrixXf& inputMat, const std::vector<bool>& clusterK)
	{
		//Eigen::VectorXf ret(rangeDim_);
		Eigen::VectorXf ret(3);
		ret.fill(0);
		int count = 0;
		for (int i=0; i<clusterK.size(); i++) if (clusterK[i]) count++;
		Eigen::MatrixXf realInput(count, rangeDim_);
		//Eigen::MatrixXf realInput(count, 3);
		count = 0;
		for (int i=0; i<clusterK.size(); i++)
		{
			if (clusterK[i])
			{
				realInput.row(count) = inputMat.row(i);
				//realInput.row(count) = MATH::spherical2Cartesian(inputMat.row(i));
				count++;
			}
		}
		Eigen::MatrixXf mat = realInput.transpose()*realInput;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(mat);
		if (eigenSolver.info()!=Eigen::Success) {
			std::cerr << "Error in SVD decomposition!\n";
			return ret;
		}
		Eigen::VectorXf eigenValues = eigenSolver.eigenvalues();
		Eigen::MatrixXf eigenVectors = eigenSolver.eigenvectors();
		int maxIdx = 0;
		float maxValue = eigenValues(maxIdx);
		for (int i=1; i<eigenValues.size(); i++)
		{
			if (eigenValues(i)>maxValue)
			{
				maxIdx = i;
				maxValue = eigenValues(maxIdx);
			}
		}
		ret = eigenVectors.col(maxIdx);
		return ret;
// 		Eigen::VectorXf retSph = MATH::cartesian2Spherical(ret, 1);
// 		return retSph.head(rangeDim_);
	}

	void ZMeshFilterManifold::setRangeWeight(int level, float weight)
	{
		if (level>=0 && level<rangeWeights_.size())
			rangeWeights_[level] = weight;
	}

	void ZMeshFilterManifold::updateRange()
	{
		output_.fill(0);
		for (int i=0; i<wki_Psi_blurs_.size(); i++)
		{
			float weight = rangeWeights_[i];
			output_.block(0, getSpatialDim(), getPointSize(), getRangeDim()) += wki_Psi_blurs_[i]*weight;
		}
	}
}