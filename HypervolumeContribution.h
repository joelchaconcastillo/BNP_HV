/*!
 * 
 *
 * \brief Implements the frontend for the HypervolumeContribution algorithms, including the approximations
 *
 *
 * \author     O.Krause
 * \date        2014-2016
 *
 *
 * \par Copyright 1995-2017 Shark Development Team
 * 
 * <BR><HR>
 * This file is part of Shark.
 * <http://shark-ml.org/>
 * 
 * Shark is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Shark is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Shark.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef SHARK_ALGORITHMS_DIRECTSEARCH_HYPERVOLUMECONTRIBUTION_H
#define SHARK_ALGORITHMS_DIRECTSEARCH_HYPERVOLUMECONTRIBUTION_H

#include "HypervolumeContribution2D.h"
#include "HypervolumeContribution3D.h"
//#include "HypervolumeContributionMD.h"
//#include "HypervolumeContributionApproximator.h"


namespace shark {
/// \brief Frontend for hypervolume contribution algorithms in m dimensions.
///
///  Depending on the dimensionality of the problem, one of the specialized algorithms is called.
///  For large dimensionalities for which there are no specialized fast algorithms,
///  the exponential time algorithm is called. 
///  Also a log-transformation of points is supported
struct HypervolumeContribution {

	/// \brief Default c'tor.
	HypervolumeContribution() : m_useApproximation(false) {}
	
	///\brief True if the hypervolume approximation is to be used in dimensions > 3.
//	void useApproximation(bool useApproximation){
//		m_useApproximation = useApproximation;
//	}
	
///	double approximationEpsilon()const{
///		return m_approximationAlgorithm.epsilon();
///	}
///	double& approximationEpsilon(){
///		return m_approximationAlgorithm.epsilon();
///	}
///	
///	double approximationDelta()const{
///		return m_approximationAlgorithm.delta();
///	}
///	
///	double& approximationDelta(){
///		return m_approximationAlgorithm.delta();
///	}
	/// \brief Returns the index of the points with smallest contribution as well as their contribution.
	///
	/// \param [in] points The set \f$S\f$ of points from which to select the smallest contributor.
	/// \param [in] k The number of points to select.
	/// \param [in] referencePoint The reference Point\f$\vec{r} \in \mathbb{R}^2\f$ for the hypervolume calculation, needs to fulfill: \f$ \forall s \in S: s \preceq \vec{r}\f$.
//	template<class Set, typename VectorType>
	std::vector<std::pair<double,std::size_t> > smallest(std::vector<std::vector<double> > const& points, std::size_t k, std::vector<double> const& ref)const{
		std::size_t numObjectives = ref.size();
		if(numObjectives == 2){
			HypervolumeContribution2D algorithm;
			return algorithm.smallest(points, k, ref);
		}
		else if(numObjectives == 3){
			HypervolumeContribution3D algorithm;
			return algorithm.smallest(points, k, ref);
		}
	//	else if(m_useApproximation){
	//		return m_approximationAlgorithm.smallest(points, k, ref);
	//	}else{
	//		HypervolumeContributionMD algorithm;
	//		return algorithm.smallest(points, k, ref);
	//	}
	}
	std::vector<std::pair<double,std::size_t> > lastIdxContributor(std::vector<std::vector<double> > const& points, std::size_t k, std::vector<double> const& ref)const{
		std::size_t numObjectives = ref.size();
		if(numObjectives == 2){
			HypervolumeContribution2D algorithm;
			return algorithm.lastIdxContributor(points, k, ref);
		}
		else if(numObjectives == 3){
			HypervolumeContribution3D algorithm;
			return algorithm.lastIdxContributor(points, k, ref);
		}
	//	else if(m_useApproximation){
	//		return m_approximationAlgorithm.smallest(points, k, ref);
	//	}else{
	//		HypervolumeContributionMD algorithm;
	//		return algorithm.smallest(points, k, ref);
	//	}
	}
	/// \brief Returns the index of the points with largest contribution as well as their contribution.
	///
	/// \param [in] points The set \f$S\f$ of points from which to select the largest contributor.
	/// \param [in] referencePoint The reference Point\f$\vec{r} \in \mathbb{R}^2\f$ for the hypervolume calculation, needs to fulfill: \f$ \forall s \in S: s \preceq \vec{r}\f$.
	std::vector<std::pair<double,std::size_t> > largest(std::vector<std::vector<double> > const& points, std::size_t k, std::vector<double> const& ref)const{
		std::size_t numObjectives = ref.size();
		if(numObjectives == 2){
			HypervolumeContribution2D algorithm;
			return algorithm.largest(points, k, ref);
		}
		else if(numObjectives == 3){
			HypervolumeContribution3D algorithm;
			return algorithm.largest(points, k, ref);
		}
	//	else{
	//		SHARK_RUNTIME_CHECK(!m_useApproximation, "Largest not implemented for approximation algorithm");
	//		HypervolumeContributionMD algorithm;
	//		return algorithm.largest(points, k, ref);
	//	}
	}

	/// \brief Returns the index of the points with smallest contribution as well as their contribution.
	///
	/// As no reference point is given, the extremum points can not be computed and are never selected.
	///
	/// \param [in] points The set \f$S\f$ of points from which to select the smallest contributor.
	/// \param [in] k The number of points to select.
	/// \param [in] referencePoint The reference Point\f$\vec{r} \in \mathbb{R}^2\f$ for the hypervolume calculation, needs to fulfill: \f$ \forall s \in S: s \preceq \vec{r}\f$.
	std::vector<std::pair<double,std::size_t> > smallest(std::vector< std::vector<double> >const& points, std::size_t k)const{

	 //	//In the Author's code (an the original paper in section 2.3.1) the reference point is defined as max(obj)+1
		std::size_t numObjectives = points[0].size(), npoints=points.size();
		std::vector<double> ref(numObjectives, -DBL_MAX);
		  for(int j = 0; j < numObjectives; j++)
	 	  {
	            for(int i = 0; i <npoints; i++)
			ref[j] = std::max(ref[j], points[i][j]);
			ref[j] +=1.0;
		  }
		if(numObjectives == 2){
			HypervolumeContribution2D algorithm;
			return algorithm.smallest(points, k, ref);
		}
		else if(numObjectives == 3){
			HypervolumeContribution3D algorithm;
			return algorithm.smallest(points, k, ref);
		}
	//	else if(m_useApproximation){
	//		return m_approximationAlgorithm.smallest(points, k, ref);
	//	}else{
	//		HypervolumeContributionMD algorithm;
	//		return algorithm.smallest(points, k, ref);
	//	}
	}
	std::vector<std::pair<double,std::size_t> > lastIdxContributor(std::vector< std::vector<double> >const& points, std::size_t k)const{

	 //	//In the Author's code (an the original paper in section 2.3.1) the reference point is defined as max(obj)+1
		std::size_t numObjectives = points[0].size(), npoints=points.size();
		std::vector<double> ref(numObjectives, -DBL_MAX);
		  for(int j = 0; j < numObjectives; j++)
	 	  {
	            for(int i = 0; i <npoints; i++)
			ref[j] = std::max(ref[j], points[i][j]);
			ref[j] +=1.0;
		  }
		if(numObjectives == 2){
			HypervolumeContribution2D algorithm;
			return algorithm.lastIdxContributor(points, k, ref);
		}
		else if(numObjectives == 3){
			HypervolumeContribution3D algorithm;
			return algorithm.lastIdxContributor(points, k, ref);
		}
	//	else if(m_useApproximation){
	//		return m_approximationAlgorithm.smallest(points, k, ref);
	//	}else{
	//		HypervolumeContributionMD algorithm;
	//		return algorithm.smallest(points, k, ref);
	//	}
	}
	/// \brief Returns the index of the points with largest contribution as well as their contribution.
	///
	/// As no reference point is given, the extremum points can not be computed and are never selected.
	///
	/// \param [in] points The set \f$S\f$ of points from which to select the smallest contributor.
	/// \param [in] referencePoint The reference Point\f$\vec{r} \in \mathbb{R}^2\f$ for the hypervolume calculation, needs to fulfill: \f$ \forall s \in S: s \preceq \vec{r}\f$.
	std::vector<std::pair<double,std::size_t> > largest(std::vector<std::vector<double> > const& points, std::size_t k)const{
		std::size_t numObjectives = points[0].size(), npoints=points.size();
		std::vector<double> ref(numObjectives, -DBL_MAX);
		  for(int j = 0; j < numObjectives; j++)
	 	  {
	            for(int i = 0; i <npoints; i++)
			ref[j] = std::max(ref[j], points[i][j]);
			ref[j] +=1.0;
		  }

		if(numObjectives == 2){
			HypervolumeContribution2D algorithm;
			return algorithm.largest(points, k);
		}
		else if(numObjectives == 3){
			HypervolumeContribution3D algorithm;
			return algorithm.largest(points, k);
		}
	//	else{
	//		HypervolumeContributionMD algorithm;
	//		return algorithm.largest(points, k);
	//	}
	}

private:
	bool m_useApproximation;
//	HypervolumeContributionApproximator m_approximationAlgorithm;
};

}
#endif
