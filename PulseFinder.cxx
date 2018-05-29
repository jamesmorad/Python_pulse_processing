
#include <cmath>
#include <algorithm>
#include "PulseFinder.h"

#define DEBUG 0
#define SMALLNUMBER 0.0001

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Independent Functions
////////////////////////////////////////////////////////////////////////////////
double SigmoidFunction(double x, double mu, double sg, double offset)  {
	return (1.-offset)/(1.+exp(-6.*(x-mu)/sg))+offset;
}

double Distance1(SplitGaussian *gl, SplitGaussian *gr) {
	return gl->BhattacharyyaDistance(*gr);
}

double Distance2(SplitGaussian *gl, SplitGaussian *gr) {
	return gl->BhattacharyyaDistance2(*gr);
}

double Distance3(SplitGaussian *gl, SplitGaussian *gr) {
	double a=0, d1=0, d2=0, area=0, da=0;
	// Power law distribution
	d1 = (gr->mode - gl->mode) - 2.58*gr->lsigma;
	if (d1 < 0.5) d1 = 0.5;
	d2 = (gr->mode - gl->mode) + 2.58*gr->rsigma;
	if (d2 < 0.5) d2 = 0.5;
	a = 1./gl->rsigma;
	area = gl->area * fabs(pow(d1, -a) - pow(d2, -a));		
	da = (gr->area - area)/(sqrt(area) + SMALLNUMBER);
	//da = da *(1.-0.1341);
	if (da < 0) da = 0;
	return da;
}

double Distance4(SplitGaussian *gl, SplitGaussian *gr) {
	SplitGaussian g = *gl + *gr;
	double max = (gl->area > gr->area)? fabs(gl->area) : fabs(gr->area);
	double sigma = (gl->area > gr->area)? gl->sigma : gr->sigma;
	double mean = (gl->area > gr->area)? gl->mean : gr->mean;
	return fabs(g.mean-mean)/(sigma/sqrt(max));
}

double PulseFinder::BDistance(int l, int r) {
  
  double d1 = Distance1(&pars[l], &pars[r]);
  double d2 = Distance2(&pars[l], &pars[r]);
  double d3 = Distance3(&pars[l], &pars[r]);
  double d4 = Distance4(&pars[l], &pars[r]);
  
  double ar = (pars[l].area < pars[r].area)? pars[l].area/pars[r].area : pars[r].area/pars[l].area;
  double hr = (pars[l].amplitude < pars[r].amplitude)? pars[l].amplitude/pars[r].amplitude : pars[r].amplitude/pars[l].amplitude;
  
  double s1scale = 0.5*(SigmoidFunction(ar,s1AMean,s1ASigma,s1AOffset) + 
                        SigmoidFunction(hr,s1HMean,s1HSigma,s1HOffset));
  double s2scale = 0.5*(SigmoidFunction(ar,s2AMean,s2ASigma,s2AOffset) + 
                        SigmoidFunction(hr,s2HMean,s2HSigma,s2HOffset));
  
  double s1sd = (d1*s1scale + d3*(1.-s1scale))*s1SDscale;
  double s2sd = (d2*s2scale + d4*(1.-s2scale))*s2SDscale;
  double seLike = sqrt((pow(pars[l].sigma/se_width_mean_samples,2)+
                     pow(pars[r].sigma/se_width_mean_samples,2))/2.);
  double seTransition = SigmoidFunction(seLike, 0.35, 0.125, 0);
  
#if DEBUG>1
cout << "\t\tpar left " << pars[l].amplitude << " " << pars[l].mode << " "
     << pars[l].lsigma << " " << pars[l].rsigma << " "
     << pars[l].area << " " << pars[l].mean << " " << pars[l].sigma << endl;
cout << "\t\tpar right " << pars[r].amplitude << " " << pars[r].mode << " "
     << pars[r].lsigma << " " << pars[r].rsigma << " "
     << pars[r].area << " " << pars[r].mean << " " << pars[r].sigma << endl;
cout << "\t\tDistance " << d1 << " " << d2 << " " << d3 << " " << d4
     << "\t" << s1scale << " " << s2scale << " " << s1sd << " " << s2sd
     << "\t" << seLike << " " << seTransition << " "<< s2sd*seTransition + s1sd*(1.-seTransition) << "\n";
cout.flush();
#endif  
  return s2sd*seTransition + s1sd*(1.-seTransition);
}

////////////////////////////////////////////////////////////////////////////////
// Pulse Finder : Main Operations Functions
////////////////////////////////////////////////////////////////////////////////
void PulseFinder::Initilize() {
	
	int entries = se_area_mean_phe;
	double width = se_width_mean_samples * sqrt(12.);
	double spWidth = 1./se_amp;
	
	double *a = new double[entries];
	double sep, amp, tMaxSep, tAvgSep, tMaxAmp, tAvgAmp;
	double maxSepMean=0, maxSepSigma=0, avgSepMean=0, avgSepSigma=0;
	double maxAmpMean=0, maxAmpSigma=0, avgAmpMean=0, avgAmpSigma=0;
	int k, l, kCount, lCount;
	
	for (int i=0; i<1000; i++) {
		
		for (int j=0; j<entries; j++) {
			a[j] = (rand() % 1000) * width/1000.;
		}
		std::sort(a, a+entries);
		
		tMaxSep = tAvgSep = 0;
		tMaxAmp = tAvgAmp = amp = se_amp;
		k = l = kCount = lCount = 0;
		for (int j=1; j<entries; j++) {
			sep = a[j] - a[j-1] - spWidth;
			if (a[j] - a[k] < spWidth) {
				amp += se_amp * (1.- (a[j]-a[k])/spWidth);
			}else {
				k = j;
				tAvgAmp += amp;
				kCount++;
				amp = se_amp;
			}
			if (a[j]-a[l] > max_seperation) {
				tAvgSep += a[j]-a[l] - spWidth;
				lCount++;
				l = j;
			}
			if (tMaxSep < sep) tMaxSep = sep;
			if (tMaxAmp < amp) tMaxAmp = amp;
		}
		if (lCount) tAvgSep /= lCount;
		if (kCount) tAvgAmp /= kCount;
		
		maxSepMean  += tMaxSep;
		maxSepSigma += tMaxSep*tMaxSep;
		avgSepMean  += tAvgSep;
		avgSepSigma += tAvgSep*tAvgSep;
		
		maxAmpMean  += tMaxAmp;
		maxAmpSigma += tMaxAmp*tMaxAmp;
		avgAmpMean  += tAvgAmp;
		avgAmpSigma += tAvgAmp*tAvgAmp;
	}
	delete[] a;
	maxSepMean /= 1000;
	maxSepSigma = sqrt(fabs(maxSepSigma-maxSepMean*maxSepMean*1000.)/999.);
	avgSepMean /= 1000;
	avgSepSigma = sqrt(fabs(avgSepSigma-avgSepMean*avgSepMean*1000.)/999.);
	
	maxAmpMean /= 1000;
	maxAmpSigma = sqrt(fabs(maxAmpSigma-maxAmpMean*maxAmpMean*1000.)/999.);
	avgAmpMean /= 1000;
	avgAmpSigma = sqrt(fabs(avgAmpSigma-avgAmpMean*avgAmpMean*1000.)/999.);
	
	double error = sqrt(pow(se_width_sigma_samples/se_width_mean_samples,2) + 
	                    pow(se_area_sigma_phe/se_area_mean_phe,2));
	maxAmp = (maxAmpMean + maxAmpSigma*nsigma)*(1+error);
  maxGap = (maxSepMean + maxSepSigma*nsigma)*(1+error);
  avgAmp = (avgAmpMean + avgAmpSigma*nsigma)*(1+error);
  avgGap = (avgSepMean + avgSepSigma*nsigma)*(1+error);
}
//______________________________________________________________________________
void PulseFinder::KZFilter() {
	// KZ filter is a repeatition of a normalized box filter 
  int lower = floor(-0.5*kz_samples);
  int upper = floor(0.5*(kz_samples+1));
  int length = wave.size();
  //Divide by 2 since LUX has good baseline. XABER does not
  //double weight, total, ty, threshold = max_threshold/2.;
  double weight, total, ty, threshold = max_threshold;
  FVec tmp;
  // remove majority of the baseline 
  smooth.push_back(0);
  for (int i=1; i<length-1; i++) {
  	//if (wave[i]>threshold || wave[i+1]>threshold || wave[i-1]>threshold) {
  		smooth.push_back(wave[i]);
  	//} 
  	//else smooth.push_back(0);
  }
  //smooth.push_back(0);
  smooth.push_back(0);
  // smooth out to minimize other baseline 
  for (int i=0; i<kz_iterations; i++) {    
    tmp = smooth;
    for (int j=0; j<length; j++) { 
      total = 0;
      ty = 0;
      for (int k=lower; k<=upper; k++) {
        weight = 0.5*(1.0-cos((2.0*M_PI*(k-lower))/(kz_samples)));
        total += weight;
        if (j+k >= 0 && j+k < length) {
          ty += tmp[j+k] * weight;
        }
      }
      smooth[j] = ty/total;
    }
  }

  smooth[0] = 0;
  for (int i=1; i<length-1; i++) {
    if (!(smooth[i] >threshold || smooth[i+1]>threshold || smooth[i-1]>threshold)) {
      smooth[i] =0;
    }
  }

  smooth[length-1] =0;



}
//______________________________________________________________________________
void PulseFinder::GetMaximums() {
	
  int width=floor(0.5*max_seperation), max_index=0, length=smooth.size();
  bool max_found, min_found;
  double min_value, max_value=0, threshold = max_threshold;
  
  for (int i=width; i<length-width; i++) {
    max_found = true;
    // 1st condition: above the max threshold
    if (smooth[i] < threshold) continue; 
    
    // 2nd condition: two maximums are seperated by a specified distance
    if (maxs.size()) {
      if (i-max_index < max_seperation) {
        if (smooth[i] > smooth[max_index]) {
          max_index = i;
          maxs[int(maxs.size())-1] = max_index;
          max_value = smooth[i];
        }
        continue;
      } 
    }
    
    // 3rd condition: highest value in the range 
    for (int j=1; j<width && max_found; j++) {
      max_found = max_found && smooth[i] >= smooth[i+j] 
                            && smooth[i] >= smooth[i-j];
    }
    if (!max_found) continue;
    
    //cout << "max potential " << i << " " << smooth[i] << endl;
    
    // 4th condition: a valid minimum exist between two maximums 
    if (maxs.size()) {    
      // Find the minimum between the current candidate and the previous max.
      // Using the square of the value to avoid sudden charge dips due to 
      // large pulses. 
      min_value = smooth[i];
      for (int j=i-1; j>max_index; j--)
        if (min_value > smooth[j]) {
          min_value = smooth[j];
        }
      // a valid minimum is when the min lower by a specified percent than the 
      // max and the difference between the max and min is above the noise.
      // This applies to both maximums. 
      min_found = min_max_ratio > fabs(min_value/smooth[i])
                  && min_max_ratio > fabs(min_value/max_value)
                  && min_value < smooth[i] - threshold
                  && min_value < max_value - threshold;
            
      if (!min_found) {
        // If no minimum is found, then one of the maximums is invalid.
        // So the larger max is kept. 
        max_found = false;
        if (max_value < smooth[i]) {
          max_index = i;
          max_value = smooth[max_index];
          maxs[int(maxs.size())-1] = max_index;
        }
      }
    } // end of looking for minimums
    
    if (max_found) {
      max_index = i;
      max_value = smooth[max_index];
      maxs.push_back(max_index);
      //i += max_seperation-1;
    }
  }
#if DEBUG>0
	cout << "Maximum threshold used " << max_threshold << endl;
  cout << "Maximum:: ";
	for (size_t i=0; i<maxs.size(); i++) 
		cout << maxs[i] << " ";
	cout << endl;
	cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::GetMinimums() {
	// The minimums include the first and last points 
	mins.push_back(0);
  double min=1e6, value;
  int idx_big, idx_small, min_idx, nmaxs=maxs.size(), step;
  int width = floor(0.5*max_seperation);
  
  for (int i=1; i<nmaxs; i++) {
    // want to search from smaller maximum towards the bigger maximum 
    if (smooth[maxs[i]] > smooth[maxs[i-1]]) {
    	idx_big = maxs[i]-width;
    	idx_small = maxs[i-1]+width;
    	step = 1;
    }else {
    	idx_big = maxs[i-1]+width;
    	idx_small = maxs[i]-width;
    	step = -1;
    }
    // find the average minimum of over the sample of points. 
    min_idx = idx_small;
    for (int j=idx_small; j!=idx_big; j+=step) {
    	value = smooth[j]+wave[j]; 
    	for (int k=1; k<width; k++) {
    		value += smooth[j+k]+wave[j+k] + smooth[j-k]+wave[j-k];
    	}
    	
    	if (min>value || j==idx_small) {
    		min = value;
    		min_idx = j;
    	}    	
    }
    mins.push_back(min_idx);
  }
  mins.push_back(int(smooth.size()));
  
  // assume the minimums give correct borders 
	borders = mins;
#if DEBUG>0
	cout << "Minimum:: ";
	for (size_t i=0; i<mins.size(); i++) 
		cout << mins[i] << " ";
	cout << endl;
	cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::BuildSplitGaussians() {
  //if (borders.size() < 2) return;
  pars.resize(borders.size()-1);
  changed.resize(borders.size()-1);

  for (size_t i=0; i<pars.size(); i++) {
    pars[i].Fit(smooth.begin()+borders[i], smooth.begin()+borders[i+1], borders[i]);
    //pars[i].Derive();
    changed[i] = true;
#if DEBUG>0
	  cout << "Build SP: " 
	       << i << " [" << borders[i] << "," << borders[i+1] << "]\t"
			   << "par(a,m,ls,rs,mm,s,sk)=(" << pars[i].area << "," << pars[i].mode 
			   << "," << pars[i].lsigma << "," << pars[i].rsigma << ","
			   << pars[i].mean << "," << pars[i].sigma << ","
			   << pars[i].skew << ")\n";
		cout.flush();
#endif
  }
}
//______________________________________________________________________________
void PulseFinder::SEFinder() {
  
  SplitGaussian g;
  size_t startIdx, endIdx;
  double distance, d1, d2, d3, d4, d5, gap;
  double prev_distance, cur_distance;
  bool pass;
  double sq3 = sqrt(3.);
  
  // loop through the SplitGaussian
  for (size_t i=0; i+1<pars.size(); i++) {
    // Set initial conditions
    startIdx = endIdx = i;
    g = pars[startIdx];
    d3 = pars[startIdx].amplitude;
    d5 = d4 = d2 = d1 =  0;
    prev_distance = cur_distance = 300;
    
    // loop through possibilities 
    for (size_t j=startIdx+1; j<pars.size(); j++) {
      g = g + pars[j];
      gap = (pars[j].mode-pars[j].lsigma*sq3) - (pars[j-1].mode+pars[j-1].rsigma*sq3);
      
      d1 = (g.area-se_area_mean_phe)/se_area_sigma_phe;
      d2 = (g.sigma-se_width_mean_samples)/se_width_sigma_samples;
      d3 = d3 < pars[j].amplitude ? pars[j].amplitude : d3;
      d4 = d4 < gap ? gap : d4;
      d5 += gap;
      distance = sqrt(pow(d1, 2) + pow(d2, 2));
      pass = g.amplitude<avgAmp && d4<maxGap && d3<maxAmp;
      pass = pass && d5/(j-startIdx)+(maxGap-avgGap)/(j-startIdx+0.5)<avgGap;
      
      if (distance<nsigma && g.amplitude<avgAmp && 
          d4<maxGap && d3<maxAmp && d5/(j-startIdx)<avgGap) {
      	prev_distance = cur_distance;
        cur_distance = distance;
        endIdx = j;
      }else if (endIdx != startIdx) { // check if the previous SE was found 
        prev_distance = cur_distance;
        cur_distance = 300;
        endIdx = j;
      }
      if (d1>=nsigma || d2>=nsigma) break;
      if (cur_distance>=nsigma && prev_distance<nsigma) break;
      
#if DEBUG>2
	  cout << "  SE Loop: "
	       << i << " " << j << " "
			   << "par(a,mm,s)=(" << g.area << "," << g.mean << "," << g.sigma << ") "
			   << "d(1,2,3,4,5)=(" << d1 << "," << d2 << "," << d3 << "," << d4 << "," 
			   << d5 << ")\n";
		cout.flush();
#endif 
    }
    
    // select the minimum distance and combine 
    if (prev_distance<nsigma || cur_distance<nsigma) {
      if (cur_distance>=nsigma) {
        endIdx--;
      }
      if (startIdx == endIdx) endIdx++;
      
      // combine 
      borders.erase(borders.begin()+startIdx+1, borders.begin()+endIdx+1);
      pars.erase(pars.begin()+startIdx+1, pars.begin()+endIdx+1);
      changed.erase(changed.begin()+startIdx+1, changed.begin()+endIdx+1);
      BuildPars(borders[startIdx], borders[startIdx+1], pars[startIdx]);
      
#if DEBUG>1
	  cout << "SEFinder: "
	       << "[" << borders[startIdx] << "," << borders[startIdx+1] << "] "
	       << startIdx << " " << endIdx << "\t"
	       << "distance(p,c,)=(" << prev_distance << ","
	       << cur_distance << ")\n";
	  cout.flush();
#endif      
      i--;
    }
  }
  
#if DEBUG>0
	cout << "SEFinder Results\n";
	for (size_t i=0; i+1<borders.size(); i++) 
	cout << i << " [" << borders[i] << "," << borders[i+1] << "]\t"
			   << "par(a,m,ls,rs)=(" << pars[i].area << "," << pars[i].mode 
			   << "," << pars[i].lsigma << "," << pars[i].rsigma << ")\n";
  cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::IteritiveCluster() {
	
	// helpful variables 
	size_t idx;
	double llDis, rrDis, lD, rD;
	bool lCombine, rCombine, lExist, llExist, rExist, rrExist;
	
	// Make list of smallest to largest pulse candidates 
	MakeOrderList();
	list<AreaOrder>::iterator it = order_list.begin();
	list<AreaOrder>::iterator it_last = order_list.end();
	
	// Loop through regions 
	for (; it != it_last && pars.size()>1;) {
		
		// ignore pulse candiates which had nothing happen around them 
		if (!changed[it->index]) continue;
		
		// setup  
		idx = it->index;
		lExist = rExist = llExist = rrExist = false;
		if (idx>0) if (pars[idx].area>pars[idx-1].area) lExist = true;
		if (idx>1) if (pars[idx-2].area>pars[idx-1].area) llExist = true; 
		if (idx+1<pars.size()) if (pars[idx].area>pars[idx+1].area) rExist = true;
		if (idx+2<pars.size()) if (pars[idx+2].area>pars[idx+1].area) rrExist = true; 
		lD = rD = llDis = rrDis = 100*nsigma;
		lCombine = rCombine = false;
		
#if DEBUG>1
	  cout << "IteritiveCluster: "
	       << idx << " [" << borders[idx] << "," << borders[idx+1] << "]\t"
			   << "par(a,m,ls,rs)=(" << pars[idx].area << "," << pars[idx].mode 
			   << "," << pars[idx].lsigma << "," << pars[idx].rsigma << ")\n";
		cout.flush();
#endif
		
		// get some quick stats 
		if (lExist) {
			lD = BDistance(idx-1, idx);
			if (llExist) llDis = BDistance(idx-2, idx-1);
			lCombine = lD<nsigma && llDis>=lD;
#if DEBUG>1
	  cout << "\tleft "
	       << lCombine << " " << lD << " " << llDis << "\n";
	  cout.flush();
#endif
		}		
		
		if (rExist) {
			rD = BDistance(idx, idx+1);
			if (rrExist) rrDis = BDistance(idx+1, idx+2);
			rCombine = rD<nsigma && rrDis>=rD;
#if DEBUG>1
	  cout << "\tright "
	       << rCombine << " " << rD << " " << rrDis << "\n";
#endif
		}
    
    // combine pulse candidates
		if (lCombine || rCombine) {
			RemoveAndUpdate(idx, lCombine, rCombine);
			it = order_list.begin();
			it_last = order_list.end();
			continue;
		}
		
		// set that this pulse candidate cannot combine with its neighbors
		changed[idx] = false;
		it++;
	}
	
	
#if DEBUG>0
	cout << "IteritiveCluster Results\n";
	for (size_t i=0; i+1<borders.size(); i++) 
	cout << i << " [" << borders[i] << "," << borders[i+1] << "]\t"
			   << "par(a,m,ls,rs)=(" << pars[i].area << "," << pars[i].mode 
			   << "," << pars[i].lsigma << "," << pars[i].rsigma << ")\n";
	cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::BaselineTrim() {
	
	// find the standard deviation of the baseline
	double threshold = max_threshold;
	int border, cur=0, last=int(borders.size()-1);
	
	// the new maximums between the borders 
	vector<int> newmaxs;
	for (size_t j=1; j<borders.size(); j++) {
		double max=0;
		int max_location = -1; 
		for (size_t k=0; k<maxs.size(); k++) {
			if (maxs[k]>borders[j]) break;
			if (maxs[k]>borders[j-1] && max<smooth[maxs[k]]) {
				max = smooth[maxs[k]];
				max_location = maxs[k];
			}
		}
		if (max_location == -1) { // if no maximum is found 
			for (int k=borders[j-1]; k<borders[j]; k++) {
				if (max<wave[k]) {
					max = wave[k];
					max_location = k;
				}
			}
		}
		newmaxs.push_back(max_location);
	}	
	
	threshold = max_threshold;
	IVec newborders;
	
	// remove prebuffer baseline  
	border = borders[0];
	for (int j=borders[0]+1; j<newmaxs[0]; j++) {
		border = j-1;
		if (wave[j] > threshold) break;
	}
	if (border-pre_buffer_samples < 0) newborders.push_back(0);
	else newborders.push_back(border-pre_buffer_samples);
	
	// any borders in between the original waveform
	for (int i=1; i<last; i++) {
		border = borders[i];
		for (int j=borders[i]-2; j>newmaxs[i-1]; j--) {
			border = j+2;
			if (wave[j] > threshold) break;
		}
		newborders.push_back(border);
		
		border = borders[i];
		for (int j=borders[i]+1; j<newmaxs[i]; j++) {
			border = j-1;
			if (wave[j] > threshold) break;
		}
		newborders.push_back(border);
		
		// adjust the boundaries to accomodate the buffer
		cur = newborders.size() - 1;
		if (border-pre_buffer_samples > newborders[cur-1]) 
			newborders[cur] = border-pre_buffer_samples;
		else newborders[cur] = newborders[cur-1]+1;
		if (newborders[cur-1]+post_buffer_samples < newborders[cur])
			newborders[cur-1] = newborders[cur-1]+post_buffer_samples;
		else newborders[cur-1] = newborders[cur]-1;
		
	}
	// remove postbuffer baseline
	border = borders[last];
	for (int j=borders[last]-2; j>newmaxs[last-1]; j--) {
		border = j+2;
		if (wave[j] > threshold) break;
	}
	border = (border+pre_buffer_samples<borders[last]) ? border+pre_buffer_samples : borders[last];
	newborders.push_back(border);
	
	borders = newborders;
	
#if DEBUG>0
	cout << "BaselineTrim: ";
	for (size_t i=0; i<borders.size(); i++) 
		cout << borders[i] << " ";
	cout << endl;
	cout.flush();
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Pulse Finder : Supporting Functions
////////////////////////////////////////////////////////////////////////////////
bool SortCompareFunc(AreaOrder a, AreaOrder b) {
	return a.area > b.area;
}
void PulseFinder::MakeOrderList() {
	order_list.clear();
	AreaOrder a;
	for (unsigned int j=0; j<pars.size(); j++) {
		if (!changed[j]) continue;
		a.area = pars[j].area; a.index = j;
		order_list.push_back(a);
	}
	order_list.sort(SortCompareFunc);
}
//______________________________________________________________________________
void PulseFinder::BuildPars(int start, int end, SplitGaussian &par) {
	par.Fit(smooth.begin()+start, smooth.begin()+end, start);
}
//______________________________________________________________________________
void PulseFinder::RemoveAndUpdate(int index, bool left, bool right) {
	
	
	if (left && right) {
		// erase pulse candidate
		borders.erase(borders.begin()+index+1);
		pars.erase(pars.begin()+index+1);
		changed.erase(changed.begin()+index+1);
		
		borders.erase(borders.begin()+index);
		pars.erase(pars.begin()+index-1);
		changed.erase(changed.begin()+index-1);
    
		index--;
		
	}else if (left) {
    // erase pulse candidate
		borders.erase(borders.begin()+index);
		pars.erase(pars.begin()+index-1);
		changed.erase(changed.begin()+index-1);
    // update changed array 
		index--;
	
	}else if (right) {
    // erase pulse candidate
		borders.erase(borders.begin()+index+1);
		pars.erase(pars.begin()+index+1);
    changed.erase(changed.begin()+index+1);
	}
	
	// update changed array 
	changed[index] = true;
	BuildPars(borders[index], borders[index+1], pars[index]);
	if (index>0) changed[index-1] = true;
	if (index>1) changed[index-2] = true;
	if (index+1<int(changed.size())) changed[index+1] = true;
	if (index+2<int(changed.size())) changed[index+2] = true;
	MakeOrderList();
}

