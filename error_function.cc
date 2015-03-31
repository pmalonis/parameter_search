#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <gsl/gsl_statistics.h>

double max_xcorr(const double x[], const double y[], 
                 const size_t n, const int max_offset) {
  /* Finds maximum cross-correlation coefficient between two arrays */

    assert(n > max_offset);
    double* data1 = new double[n];
    double* data2 = new double[n];
    double* offset = new double[n];
    double* corrs = new double[max_offset + 1];
    //    std::vector<double> corrs;
    //corrs.reserve(max_offset + 1);
    for (int i = 0; i <= max_offset; i++){
        memset(data1, 0, n);
        memset(data2, 0, n);
        memcpy(data1 + i, x, (n-i)*sizeof(double));
        memcpy(data2, y + i, (n-i)*sizeof(double));
        //corrs.push_back(gsl_stats_correlation(data1, 1, data2, 1, n));
        corrs[i] = gsl_stats_correlation(x + i, 1, y, 1, n-i);
    }   
    
    delete[] data1;
    delete[] data2; 
    delete[] offset;
    delete[] corrs;
    return *std::max_element(corrs, corrs + max_offset + 1);
}

// double* get_subthreshold(const double v[], double (&sub_threshold)[],
//                          const size_t n, const double threshold){
//     /* returns array with spikes removed 
//      Arguments:
//      v: Voltage data
//      n: size of v
//      thresh: threshold of local maximum to be considered a spike
//      window: window (in samples) around spikes to remove*/
    
//     std::vector<double> sub_threshold;
//     //    std::vector<int>::const_iterator it;
//     for (int i = 0; i < n; i++) {
//         if (v[i] < threshold) 
//             sub_threshold.push_back(v[i]);
//     }
//     return sub_threshold.data();
// }

// double compare_subthreshold(const double x[], const double y[],const double threshold, 
//                             const size_t n, const int max_offset){
//     /* Compares subthreshold voltage by finding the maximum cross correlation*/
//     double x_subthreshold = get_subthreshold(x, n, threshold);
//     double y_subthreshold[] = get_subthreshold(y, n, threshold);
//     return max_xcorr(x_subthreshold, y_subthreshold, n, max_offset);
// }

std::vector<int> find_spikes(const double v[], const int n, const double threshold) {
    std::vector<int> spikes;
    for (int i = 1; i < n-1 ; i++){
        if ((v[i] > v[i-1]) & (v[i] > v[i+1]) & (v[i] > threshold)) {
            spikes.push_back(i);
        }
    }
    return spikes; 
}

// std::vector<std::vector<double> >
// spike_waveforms(const double v[], const int n) {
//     size_t half_window = 20;
//     std::vector<int> spikes = find_spikes(v, n);
//     std::vector<double> temp_waveform;
//     std::vector<std::vector<double> > waveforms;
//     for (std::vector<int>::const_iterator it = spikes.begin(); 
//          it != spikes.end(); it++) {
//         temp_waveform.assign(it-half_window, it+half_window);
//         waveforms.push_back(temp_waveform);
//     }

//     return waveforms;
// }

// int main() {
//     double x[] = {0,5,2,3,4};
//     double y[] = {5,2,3,4,0};
//     double r = max_xcorr(x, y, 5, 3);
//     //spike_waveforms(x, 5);
      
//     return 0;

// }
