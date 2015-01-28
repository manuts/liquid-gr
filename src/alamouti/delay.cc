#include "alamouti.h"

namespace liquid {
  namespace alamouti {

    tx_buff_transformer::tx_buff_transformer()
    {
      items_txed = 0;
      hist = 0;
    }

    tx_buff_transformer::~tx_buff_transformer()
    {
      if(initialized) {
        free(history[0]);
        free(history[1]);
      }
    }

    void tx_buff_transformer::set_usrp_buff_len(unsigned int i)
    {
      usrp_buff_len = i;
      history[0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*usrp_buff_len);
      history[1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*usrp_buff_len);
      initialized = true;
    }

    void tx_buff_transformer::set_packer_buff_len(unsigned int i)
    {
      packer_buff_len = i;
    }

    void tx_buff_transformer::set_packer_buff(std::complex<float> ** _packer_buff)
    {
      packer_buff[0] = _packer_buff[0];
      packer_buff[1] = _packer_buff[1];
    }

    void tx_buff_transformer::set_usrp_buff(std::complex<float> ** _usrp_buff)
    {
      usrp_buff[0] = _usrp_buff[0];
      usrp_buff[1] = _usrp_buff[1];
    }

    void tx_buff_transformer::copy_buff()
    {
      memmove(usrp_buff[0], history[0], sizeof(std::complex<float>)*hist);
      memmove(usrp_buff[1], history[1], sizeof(std::complex<float>)*hist);
      items_txed = (items_txed + hist)%packer_buff_len;
      memmove(usrp_buff[0] + hist, packer_buff[0] + items_txed, sizeof(std::complex<float>)*(usrp_buff_len - hist));
      memmove(usrp_buff[1] + hist, packer_buff[1] + items_txed, sizeof(std::complex<float>)*(usrp_buff_len - hist));
      items_txed = (items_txed + usrp_buff_len - hist)%packer_buff_len;
      if(packer_buff_len - items_txed < usrp_buff_len) {
        hist = packer_buff_len - items_txed;
        memmove(history[0], packer_buff[0] + items_txed, sizeof(std::complex<float>)*(hist));
        memmove(history[1], packer_buff[1] + items_txed, sizeof(std::complex<float>)*(hist));
      }
      else
      {
        hist = 0;
      }
    }

    bool tx_buff_transformer::fill_buff()
    {
      return((packer_buff_len - items_txed) < usrp_buff_len);
    }

    void rx_buff_transformer::copy_buff()
    {
      memmove(unpacker_buff[0] + items_rxed, history[0], sizeof(std::complex<float>)*hist);
      memmove(unpacker_buff[1] + items_rxed, history[0], sizeof(std::complex<float>)*hist);
      items_rxed = (items_rxed + hist)%unpacker_buff_len;
      if(unpacker_buff_len - items_rxed < usrp_buff_len) {
        memmove(unpacker_buff[0] + items_rxed, usrp_buff[0], sizeof(std::complex<float>)*usrp_buff_len);
        memmove(unpacker_buff[1] + items_rxed, usrp_buff[1], sizeof(std::complex<float>)*usrp_buff_len);
        items_rxed = (items_rxed + usrp_buff_len)%unpacker_buff_len;
        hist = 0;
      }
      else {
        memmove(unpacker_buff[0] + items_rxed, usrp_buff[0], sizeof(std::complex<float>)*(unpacker_buff_len - items_rxed));
        memmove(unpacker_buff[1] + items_rxed, usrp_buff[1], sizeof(std::complex<float>)*(unpacker_buff_len - items_rxed));
        hist = usrp_buff_len - (unpacker_buff_len - items_rxed);
        memmove(history[0], usrp_buff[0] + unpacker_buff_len - items_rxed, sizeof(std::complex<float>)*hist);
        memmove(history[1], usrp_buff[1] + unpacker_buff_len - items_rxed, sizeof(std::complex<float>)*hist);
        items_rxed = 0;
      }
    }

    bool rx_buff_transformer::dump_buff()
    {
      return((items_rxed + usrp_buff_len) > unpacker_buff_len);
    }

    rx_buff_transformer::rx_buff_transformer()
    {
      items_rxed = 0;
      hist = 0;
    }

    rx_buff_transformer::~rx_buff_transformer()
    {
      if(initialized)
        free(history);
    }

    void rx_buff_transformer::set_usrp_buff_len(unsigned int i)
    {
      usrp_buff_len = i;
      history[0] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*usrp_buff_len);
      history[1] = (std::complex<float> *)malloc(sizeof(std::complex<float>)*usrp_buff_len);
      initialized = true;
    }

    void rx_buff_transformer::set_unpacker_buff_len(unsigned int i)
    {
      unpacker_buff_len = i;
    }

    void rx_buff_transformer::set_unpacker_buff(std::complex<float> ** _packer_buff)
    {
      unpacker_buff[0] = _packer_buff[0];
      unpacker_buff[1] = _packer_buff[1];
    }

    void rx_buff_transformer::set_usrp_buff(std::complex<float> ** _usrp_buff)
    {
      usrp_buff[0] = _usrp_buff[0];
      usrp_buff[1] = _usrp_buff[1];
    }
    
    delay::delay(unsigned int _d)
    {
      d = _d;
      hist = (std::complex<float> *)malloc(sizeof(std::complex<float>)*d);
      for(unsigned int i = 0; i < d; i++)
        hist[i] = 0.0f;
    }

    delay::~delay()
    {
      free(hist);
    }

    void delay::set_delay(unsigned int _d)
    {
      free(hist);
      hist = (std::complex<float> *)malloc(sizeof(std::complex<float>)*d);
      for(unsigned int i = 0; i < d; i++)
        hist[i] = 0.0f;
    }

    void delay::work(std::complex<float> * in, unsigned int num_items)
    {
      std::complex<float> * temp = 
        (std::complex<float> *)malloc(sizeof(std::complex<float>)*d);
      memmove(temp, in + num_items - d, sizeof(std::complex<float>)*d);
      memmove(in + d, in, sizeof(std::complex<float>)*(num_items - d));
      memmove(in, hist, sizeof(std::complex<float>)*d);
      memmove(hist, temp, sizeof(std::complex<float>)*d);
      free(temp);
    }
  }
}
