#ifndef __PULSE_FINDER_H__
#define __PULSE_FINDER_H__ 1

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <list>
#include "SplitGaussian.h"

typedef std::vector<float> FVec;
typedef std::vector<double> DVec;
typedef std::vector<int> IVec;

typedef struct {
	double area;
	int index;
}AreaOrder;

double Overlap(SplitGaussian &lhs, SplitGaussian &rhs, double b1, double b2);
double Distance1(SplitGaussian *gl, SplitGaussian *gr);
double Distance2(SplitGaussian *gl, SplitGaussian *gr);
double Distance3(SplitGaussian *gl, SplitGaussian *gr);
double Distance4(SplitGaussian *gl, SplitGaussian *gr);

class PulseFinder {
	protected:
		// Supporting functions
		void MakeOrderList();
		void BuildPars(int start, int end, SplitGaussian &g);
		void RemoveAndUpdate(int index, bool left, bool right);
		
		double maxAmp;
		double maxGap;
		double avgAmp;
		double avgGap;
		std::list<AreaOrder> order_list;
		std::vector<bool> changed;
		
	public:
		// Parameters
		int kz_samples;
		int kz_iterations;
		double min_max_ratio;
		double max_threshold;
		int max_seperation;
		double nsigma;
		int pre_buffer_samples;
		int post_buffer_samples;
		double se_amp;
		double se_area_mean_phe;
		double se_area_sigma_phe;
		double se_width_mean_samples;
		double se_width_sigma_samples;
		double s1AMean,s1ASigma,s1AOffset;
		double s1HMean,s1HSigma,s1HOffset;
		double s2AMean,s2ASigma,s2AOffset;
		double s2HMean,s2HSigma,s2HOffset;
		double s1SDscale, s2SDscale;
		
		// Containers
		FVec wave;
		FVec smooth;
		IVec maxs;
		IVec mins;
		IVec borders;
		std::vector<SplitGaussian> pars;
		
		void KZFilter();
		void GetMaximums();
		void GetMinimums();
		void BuildSplitGaussians();
		void SEFinder();
		void IteritiveCluster();
		void BaselineTrim();
		double BDistance(int left, int right);
		PulseFinder() { 
			kz_samples=10; 
			kz_iterations=2; 
			min_max_ratio=0.6; 
			max_threshold=0.15;
			max_seperation=8; 
			nsigma=3.0; 
			pre_buffer_samples=0;
			post_buffer_samples=30;
			se_amp=0.281; 
			se_area_mean_phe=22.88;
			se_area_sigma_phe=5.673; 
			se_width_mean_samples=30.74;
			se_width_sigma_samples=5.635;
			s1AMean = 0.272;
			s1ASigma = 0.00474;
			s1AOffset = 0.255;
			s1HMean = 0.978;
			s1HSigma = 46.2;
			s1HOffset = -0.224;
			s1SDscale = 1.657;
			s2AMean = 0.774;
			s2ASigma = 1.23;
			s2AOffset = 0.132;
			s2HMean = 0.898;
			s2HSigma = 10.0;
			s2HOffset = 0.0938;
			s2SDscale = 0.1279;
			
		}
		void Initilize();
		void Execute() {
			Clear();
			KZFilter();
			GetMaximums();
			GetMinimums();
			BuildSplitGaussians();
			SEFinder();
			IteritiveCluster();
			BaselineTrim();
		}
		void Clear() {
			smooth.clear();
			maxs.clear();
			mins.clear();
			borders.clear();
			pars.clear();
			order_list.clear();
			changed.clear();
		}
};

#endif

