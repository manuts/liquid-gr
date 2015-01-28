#include "alamouti.h"

namespace liquid {
  namespace alamouti {
    int channel_estimator::find_corr_index(std::complex<float> * buff_ptr,
                                           unsigned int approximate_index)
    {
      std::complex<float> corr;

      glo_max = 0;
      for(int i = -RANGE; i < RANGE + 1; i++) {
        float rxy_max = 0;
        for (unsigned int k = 0; k < m; k++) {
          dotprod_cccf_execute(dp[k], buff_ptr + i + approximate_index, &corr);
          rxy[k] = std::abs(corr);
          if (rxy[k] > rxy_max) {
            rxy_max = rxy[k];
            loc_max_index = k;
          }
        }
        if(rxy[loc_max_index] > glo_max) {
          glo_max = rxy[loc_max_index];
          glo_max_index = i;
          glo_max_freq_index = loc_max_index;
        }
      }
      return glo_max_index;
    }

    channel_estimator::channel_estimator(std::complex<float> * _training_seq,
                                         unsigned int _training_seq_len,
                                         float _threshold,
                                         float _dphi_max)
    {
      training_seq_len = _training_seq_len;
      training_seq = (std::complex<float> *)malloc(sizeof(std::complex<float>)*training_seq_len);
      memmove(training_seq, _training_seq, sizeof(std::complex<float>)*training_seq_len);
      threshold = _threshold;
      dphi_max = _dphi_max;

      m = NUM_FILTS;
      dphi_step = 0.8f * M_PI / float(training_seq_len);

      dp   = (dotprod_cccf*) malloc(m*sizeof(dotprod_cccf));

      std::complex<float> sconj[training_seq_len];

      for(unsigned int k = 0; k < m; k++)
      {
        dphi[k] = ((float)k - (float)(m - 1)/2) * dphi_step;
        for(unsigned int i = 0; i < training_seq_len; i++)
        {
          sconj[i] = std::conj(training_seq[i]) * std::exp(liquid::math::I*dphi[k]*(float)i);
        }
        dp[k] = dotprod_cccf_create(sconj, training_seq_len);
      }
    }

    channel_estimator::~channel_estimator()
    {
      free(training_seq);
      for(unsigned int i = 0; i < m; i++)
        dotprod_cccf_destroy(dp[i]);
      free(dp);
    }
  }
}
