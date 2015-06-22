/*
 * Copyright (c) 2014, 2015 Manu T S
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <uhd/usrp/multi_usrp.hpp>
#include <liquid/ofdm.h>
#include <iostream>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <ctime>

// define default configurations
// USRP identities
#define N200_12             "addr=134.147.118.212"
#define N200_15             "addr=134.147.118.215"
#define X300A               "addr=134.147.118.216"
#define X300B               "addr=134.147.118.217"
#define B210A               "serial=308F955"
#define B210B               "serial=308F965"

// subdevice specifications
#define X300_SUBDEV_SPEC    "A:0 B:0"
#define N200_SUBDEV_SPEC    "A:0"
#define B210_SUBDEV_SPEC    "A:A A:B"

// tx/rx streamer configurations
#define CPU                 "fc32"
#define WIRE                "sc16"

// RF front end configurations
#define CENTER_FREQUENCY    5100e6
#define SAMPLING_RATE       1e6
#define TX_FRONTEND_GAIN    45.0
#define RX_FRONTEND_GAIN    45.0

// device synchronization configurations
typedef enum
{
  CLOCK_SOURCE_NONE = 0,
  CLOCK_SOURCE_EXTERNAL
} clock_source_type;

typedef enum
{
  TIME_SOURCE_NONE = 0,
  TIME_SOURCE_EXTERNAL
} time_source_type;
#define CLOCK_SOURCE        CLOCK_SOURCE_NONE
#define TIME_SOURCE         TIME_SOURCE_NONE

// OFDM configurations
#define NUM_SUBCARRIERS     64
#define CP_LENGTH           16
#define BASEBAND_GAIN       0.25
#define TAPER_LEN           0

// misc configurations
#define VERBOSITY           true
#define RUNTIME             10.0
#define SIMULATION          false

namespace po = boost::program_options;

// structure to hold command line inputs
typedef struct
{
  double cent_freq;         // center frequency of transmission
  double samp_rate;         // usrp samping rate
  float dsp_gain;           // dsp gain
  double txgain;            // tx frontend gain
  double rxgain;            // rx frontend gain
  unsigned int M;           // number of subcarriers
  unsigned int cp_len;      // length of cyclic prefix
  bool verbose;             // verbosity of the app
} options;

// structure to hold usrp configurations
typedef struct
{
  double cent_freq;               // center frequency
  double samp_rate;               // sampling ratea
  float rf_gain;                  // sampling ratea
  clock_source_type clock_source; // clock source TODO List the sources
  time_source_type time_source;   // time source TODO List the sources 
} usrp_config;

// structure to load data to simulator
typedef struct
{
  liquid::ofdm::modulator * mod;
  liquid::ofdm::demodulator * dem;
  void * usr_data;
} simulator_data;

// structure to load data to tx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * tx;
  liquid::ofdm::modulator * mod;
  time_t * tx_begin;
  time_t * tx_end;
  unsigned int * pid;
  float runtime;
} tx_thread_data;

// structure to load data to rx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * rx;
  liquid::ofdm::demodulator * dem;
  time_t * rx_begin;
  time_t * rx_end;
  float runtime;
  bool * streamer_error;
  uhd::rx_metadata_t * rxmd_to_main;
  bool online_decoding;
  float * time_for_offline_decoding;
} rx_thread_data;

// global counters
unsigned int num_frames_detected;
unsigned int num_valid_headers_received;
unsigned int num_valid_bytes_received;
unsigned int num_valid_packets_received;

int callback(unsigned char *  _header,
             int              _header_valid,
             unsigned char *  _payload,
             unsigned int     _payload_len,
             int              _payload_valid,
             framesyncstats_s _stats,
             void *           _userdata)
{
  // update global counters
  num_frames_detected++;

  if (_header_valid)
      num_valid_headers_received++;

  if (_payload_valid) {
      num_valid_packets_received++;
      num_valid_bytes_received += _payload_len;
  }
  return 0;
}

// function to read commandline options
// argc: argument count
// argv: argument list
// parameter options: pointer to a struct_options object
void read_options(int       argc,
                  char  **  argv,
                  options * d_options)
{
  //set the operating parameters
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "help message")
    ("freq,f",
     po::value<double>(&(d_options->cent_freq))->default_value(
                                            CENTER_FREQUENCY),
     "RF center frequency in Hz")
    ("rate,r",
     po::value<double>(&(d_options->samp_rate))->default_value(
                                            SAMPLING_RATE),
     "USRP Sampling rate")
    ("dsp_gain",
     po::value<float>(&(d_options->dsp_gain))->default_value(
                                          BASEBAND_GAIN),
     "TX DSP gain")
    ("tx_gain",
     po::value<double>(&(d_options->txgain))->default_value(
                                         TX_FRONTEND_GAIN),
     "TX Front end gain")
    ("rx_gain",
     po::value<double>(&(d_options->rxgain))->default_value(
                                         RX_FRONTEND_GAIN),
     "RX Front end gain")
    ("num_subcarriers",
     po::value<unsigned int>(&(d_options->M))->default_value(
                                          NUM_SUBCARRIERS),
     "Number of OFDM subcarriers")
    ("cp_len",
     po::value<unsigned int>(&(d_options->cp_len))->default_value(
                                               CP_LENGTH),
     "Cyclic Prefix Length")
    ("verbose,v",
     "Verbose")
    ("quite,q",
     "Quite")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  //print the help message
  if (vm.count("help")) {
    std::cout << boost::format("ofdmtxrx %s") % desc << std::endl;
    std::cout
      << std::endl
      << "Basic OFDM P2P Link\n"
      << std::endl;
    exit(0);
  }

  //check sanity of options
  if (vm.count("verbose") && vm.count("quite")) {
    std::cout << "Conflicting Options Verbose and Quite."
      << " Please use only one of those."
      << std::endl;
    exit(0);
  }

  if (vm.count("verbose"))
    d_options->verbose = true;
}

void init_usrp(uhd::usrp::multi_usrp::sptr * u,
               usrp_config * configs,
               bool is_tx)
{
  size_t num_chans;
  switch(configs->clock_source) {
    case(CLOCK_SOURCE_NONE):
      break;
    case(CLOCK_SOURCE_EXTERNAL):
      (*u)->set_clock_source("external");
      break;
    default:
      std::cout << "Clock source not known. Exiting\n";
      exit(1);
  }
  switch(configs->time_source) {
    case(TIME_SOURCE_NONE):
      break;
    case(TIME_SOURCE_EXTERNAL):
      (*u)->set_time_source("external");
      break;
    default:
      std::cout << "Time source not known. Exiting\n";
      exit(1);
  }

  if(is_tx) {
    // set subdev specs
    if((*u)->get_mboard_name() == "X300") {
      (*u)->set_tx_subdev_spec(
          uhd::usrp::subdev_spec_t(X300_SUBDEV_SPEC),
          uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if(((*u)->get_mboard_name() == "N200") || 
            ((*u)->get_mboard_name() == "N200r4")){
      (*u)->set_tx_subdev_spec(
          uhd::usrp::subdev_spec_t(N200_SUBDEV_SPEC),
          uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if((*u)->get_mboard_name() == "B210") {
      (*u)->set_tx_subdev_spec(
          uhd::usrp::subdev_spec_t(B210_SUBDEV_SPEC),
          uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else {
      std::cout << "TX Motherboard not compatible\n"
                << "Subdevice specification for "
                << (*u)->get_mboard_name()
                << " not known. Exiting\n";
      exit(1);
    }
    num_chans = (*u)->get_tx_num_channels();
    for (size_t chan = 0; chan < num_chans; chan++) {
      (*u)->set_tx_rate(configs->samp_rate, chan);
      (*u)->set_tx_gain(configs->rf_gain, chan);
      (*u)->set_tx_antenna("TX/RX", chan);
    }
    uhd::time_spec_t cmd_time = (*u)->get_time_now() + 
                                uhd::time_spec_t(0.1);
    (*u)->set_command_time(cmd_time);
    uhd::tune_request_t tx_tune_request(configs->cent_freq);
    for (size_t chan = 0; chan < num_chans; chan++) 
      (*u)->set_tx_freq(tx_tune_request, chan);
    (*u)->clear_command_time();
  }
  else {
    // set subdev specs
    if((*u)->get_mboard_name() == "X300") {
      (*u)->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
            X300_SUBDEV_SPEC),
            uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if(((*u)->get_mboard_name() == "N200") ||
            ((*u)->get_mboard_name() == "N200r4")) {
      (*u)->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
            N200_SUBDEV_SPEC),
            uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if((*u)->get_mboard_name() == "B210") {
      (*u)->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
            B210_SUBDEV_SPEC),
            uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else {
      std::cout << "RX Motherboard not compatible\n"
                << "Subdevice specification for "
                << (*u)->get_mboard_name()
                << " not known. Exiting\n";
      exit(1);
    }
    num_chans = (*u)->get_rx_num_channels();
    for (size_t chan = 0; chan < num_chans; chan++) {
      (*u)->set_rx_rate(configs->samp_rate, chan);
      (*u)->set_rx_gain(configs->rf_gain, chan);
      (*u)->set_rx_antenna("TX/RX", chan);
    }
    uhd::time_spec_t cmd_time = (*u)->get_time_now() + 
                                uhd::time_spec_t(0.1);
    (*u)->set_command_time(cmd_time);
    uhd::tune_request_t rx_tune_request(configs->cent_freq);
    for (size_t chan = 0; chan < num_chans; chan++) 
      (*u)->set_rx_freq(rx_tune_request, chan);
    (*u)->clear_command_time();
  }
}

void * simulator (void * _data)
{
  std::cout << "Running simulator\n";
  unsigned int num_frames = 100;
  simulator_data * data = (simulator_data *)_data;
  // assertions
  assert((data->mod)->get_M() 
      == (data->dem)->get_M());
  assert((data->mod)->get_frame_len() 
      == (data->dem)->get_frame_len());
  assert((data->mod)->get_payload_mod_len() 
      == (data->dem)->get_payload_mod_len());
  assert((data->mod)->get_payload_enc_len() 
      == (data->dem)->get_payload_enc_len());
  assert((data->mod)->get_cp_len() 
      == (data->dem)->get_cp_len());

  unsigned int frame_len;
  unsigned int M;
  unsigned int cp_len;
  unsigned int sig_buff_len;
  std::complex<float> * sig_buff;
  // payload len
  unsigned int payload_len;
  // header len
  unsigned int header_len;
  // payload data
  unsigned char * payload;
  // header data
  unsigned char * header;
  unsigned int pid;

  frame_len = (data->mod)->get_frame_len();
  M = (data->mod)->get_M();
  cp_len = (data->mod)->get_cp_len();
  sig_buff_len = (M + cp_len)*frame_len;
  header_len = (data->mod)->get_h_usr_len();
  payload_len = (data->mod)->get_payload_dec_len();
  payload = (unsigned char *) malloc 
    (payload_len*sizeof(unsigned char));
  header = (unsigned char *) malloc 
    (header_len*sizeof(unsigned char));

  sig_buff = (std::complex<float> *) malloc 
                 (sizeof(std::complex<float>)*sig_buff_len);

  unsigned int i;
  for(pid = 0; pid < num_frames; pid++)
  {
    printf("tx packet id: %6u\n", pid);
    header[0] = (pid >> 8) & 0xff;
    header[1] = (pid     ) & 0xff;
    for(i = 2; i < header_len; i++)
      header[i] = rand() & 0xff;

    for(i = 0; i < payload_len; i++)
      payload[i] = rand() & 0xff;

    (data->mod)->assemble_frame(header, payload);
    (data->mod)->assemble_output_samples(sig_buff);
    (data->dem)->demodulate_samples(sig_buff, sig_buff_len);
  }
  free(sig_buff);
  return NULL;
}

void * rx_worker (void * _data)
{
  // vector of pointers to the modulator buffers.
  // TODO modifications for multichannel modulators
  std::complex<float> ** demodulator_buffer;
  // usrp buffer
  std::vector<std::complex<float> *> usrp_buffer;
  // usrp buffer len
  size_t usrp_buffer_len;
  // number of channels available on USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples received in each call of recv
  unsigned int num_samples_rcvd;
  // timeout
  float timeout = 0.2;
  // file to save data;
  FILE * f_rx_sig;
  // total number of samples received
  unsigned int num_accumulated_samples;

  rx_thread_data * data = (rx_thread_data *)_data;

  num_channels = (*(data->rx))->get_rx_num_channels();

  if(!(data->online_decoding))
    f_rx_sig = fopen("/tmp/rx_sig", "wb");

  for(size_t chan = 0; chan < num_channels; chan++)
    channels.push_back(chan);
  uhd::stream_args_t rx_stream_args(CPU, WIRE);
  rx_stream_args.channels = channels;
  uhd::rx_streamer::sptr rx_stream = 
    (*(data->rx))->get_rx_stream(rx_stream_args);

  usrp_buffer_len = rx_stream->get_max_num_samps();
  demodulator_buffer = (std::complex<float> **) 
    malloc(num_channels*sizeof(std::complex<float> *));
  for(size_t chan = 0; chan < num_channels; chan++) {
    usrp_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*usrp_buffer_len));
    demodulator_buffer[chan] = usrp_buffer[chan];
  }

  uhd::rx_metadata_t rxmd;
  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.stream_now = false;
  stream_cmd.time_spec = (*(data->rx))->get_time_now() 
                       + uhd::time_spec_t(0.1);
  rx_stream->issue_stream_cmd(stream_cmd);
  *(data->streamer_error) = false;
  *(data->rx_begin) = time(NULL);
  num_accumulated_samples = 0;

  while(time(NULL) - *(data->rx_begin) < data->runtime)
  {
    num_samples_rcvd = rx_stream->recv(usrp_buffer,
                                       usrp_buffer_len,
                                       rxmd,
                                       timeout);
    if(rxmd.error_code) {
      *(data->streamer_error) = true;
      memmove(data->rxmd_to_main, &rxmd,
              sizeof(uhd::rx_metadata_t));
      break;
    }
    timeout = 0.1;
    num_accumulated_samples += num_samples_rcvd;
    if(data->online_decoding) {
      (data->dem)->demodulate_samples(demodulator_buffer[0],
                                      num_samples_rcvd);
    }
    else {
      assert(fwrite((void *)demodulator_buffer[0],
                    sizeof(std::complex<float>),
                    num_samples_rcvd,
                    f_rx_sig) ==
             num_samples_rcvd);
    }
  }

  *(data->rx_end) = time(NULL);
  std::cout << "Exiting rx thread\n";
  // stop rx streaming
  stream_cmd.stream_mode = 
    uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_stream->issue_stream_cmd(stream_cmd);

  if(!(data->online_decoding)) {
    fclose(f_rx_sig);
    unsigned int samples_remaining = num_accumulated_samples;
    unsigned int num_samples_to_process;
    f_rx_sig = fopen("/tmp/rx_sig", "rb");
    time_t decode_begin = time(NULL);
    while(samples_remaining)
    {
      num_samples_to_process = 
        (samples_remaining > usrp_buffer_len) ?
        usrp_buffer_len :
        samples_remaining;
      assert(
          fread((void *)demodulator_buffer[0],
                sizeof(std::complex<float>),
                num_samples_to_process,
                f_rx_sig) == 
                num_samples_to_process);
      (data->dem)->demodulate_samples(demodulator_buffer[0],
                                      num_samples_to_process);
      samples_remaining -= num_samples_to_process;
    }
    *(data->time_for_offline_decoding) = 
      time(NULL) - decode_begin;
  }

  for(size_t chan = 0; chan < num_channels; chan++) {
    free(usrp_buffer[chan]);
  }
  free(demodulator_buffer);
  pthread_exit(NULL);
}

void * tx_worker (void * _data)
{
  // number of subcarriers
  unsigned int M;
  // cyclic prefix length
  unsigned int cp_len;
  // number of OFDM symbols per packet
  unsigned int frame_len;
  // number of samples per frame at
  // the output of modulator
  unsigned int modulator_buffer_len;
  // vector of pointers to the modulator buffers.
  // TODO modifications for multichannel modulators
  std::complex<float> ** modulator_buffer;
  // usrp buffer
  std::vector<std::complex<float> *> usrp_buffer;
  // payload len
  unsigned int payload_len;
  // header len
  unsigned int header_len;
  // payload data
  unsigned char * payload;
  // header data
  unsigned char * header;
  // number of channels available with USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples sent in each call of send
  unsigned int num_samples_sent;
  
  tx_thread_data * data = (tx_thread_data *)_data;

  M = (data->mod)->get_M();
  cp_len = (data->mod)->get_cp_len();
  frame_len = (data->mod)->get_frame_len();
  modulator_buffer_len = (M + cp_len)*frame_len;
  header_len = (data->mod)->get_h_usr_len();
  payload_len = (data->mod)->get_payload_dec_len();
  num_channels = (*(data->tx))->get_tx_num_channels();

  payload = (unsigned char *) malloc 
    (payload_len*sizeof(unsigned char));
  header = (unsigned char *) malloc 
    (header_len*sizeof(unsigned char));
  modulator_buffer = (std::complex<float> **) 
    malloc(num_channels*sizeof(std::complex<float> *));
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    modulator_buffer[chan] = (std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*modulator_buffer_len);
    usrp_buffer.push_back(modulator_buffer[chan]);
  }

  uhd::stream_args_t tx_stream_args(CPU, WIRE);
  tx_stream_args.channels = channels;
  uhd::tx_streamer::sptr tx_stream = 
    (*(data->tx))->get_tx_stream(tx_stream_args);

  uhd::tx_metadata_t txmd;
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.1);
  *(data->tx_begin) = time(NULL);

  unsigned int i;
  *(data->pid) = 0;
  while(time(NULL) - *(data->tx_begin) < data->runtime)
  {
//    printf("tx packet id: %6u\n", *(data->pid));
    header[0] = (*(data->pid) >> 8) & 0xff;
    header[1] = (*(data->pid)     ) & 0xff;
    for(i = 2; i < header_len; i++)
      header[i] = rand() & 0xff;

    for(i = 0; i < payload_len; i++)
      payload[i] = rand() & 0xff;

    (data->mod)->assemble_frame(header, payload);
    (data->mod)->assemble_output_samples(modulator_buffer[0]);
    num_samples_sent = tx_stream->send(usrp_buffer,
                                       modulator_buffer_len,
                                       txmd);
    assert(num_samples_sent == modulator_buffer_len);
    txmd.start_of_burst = false;
    txmd.has_time_spec = false;
    (*(data->pid))++;
  }

  // transmission ends
  *(data->tx_end) = time(NULL);
  std::cout << "Exiting tx thread\n";
  // send the last packet
  txmd.end_of_burst = true;
  tx_stream->send(usrp_buffer,
                  0,
                  txmd);

  for(size_t chan = 0; chan < num_channels; chan++) {
    free(modulator_buffer[chan]);
  }
  free(modulator_buffer);

  pthread_exit(NULL);
}

int UHD_SAFE_MAIN(int argc, char **argv)
{
  uhd::set_thread_priority_safe();

  options d_options;
  usrp_config * tx_config;
  usrp_config * rx_config;
  simulator_data * modem;
  tx_thread_data * tx_data;
  rx_thread_data * rx_data;
  uhd::usrp::multi_usrp::sptr tx;
  uhd::usrp::multi_usrp::sptr rx;
  pthread_t tx_thread, rx_thread;

  read_options(argc, argv, &d_options);
  std::complex<float> tx_gain(d_options.dsp_gain, 0.0);
  bool rx_streamer_error;
  uhd::rx_metadata_t rxmd_from_worker;
  float time_for_offline_decoding;

  unsigned char * p = NULL;
  void * user_data = (void *)&(d_options.samp_rate);

  time_t tx_begin, tx_end, rx_begin, rx_end;
  unsigned int pid;
  bool online_decoding = true;

  OFDMFRAME_STRUCT * frame_struct;
  frame_struct = (OFDMFRAME_STRUCT *) malloc 
    (sizeof(OFDMFRAME_STRUCT));
  frame_struct->VERSN =  104;
  frame_struct->H_USR =  2;
  frame_struct->H_LEN =  frame_struct->H_USR + 6;
  frame_struct->H_EC1 =  LIQUID_FEC_GOLAY2412;
  frame_struct->H_EC2 =  LIQUID_FEC_NONE;
  frame_struct->H_CRC =  LIQUID_CRC_32;
  frame_struct->H_MOD =  LIQUID_MODEM_BPSK;
  frame_struct->P_LEN =  512;
  frame_struct->P_EC1 =  LIQUID_FEC_NONE;
  frame_struct->P_EC2 =  LIQUID_FEC_NONE;
  frame_struct->P_CRC =  LIQUID_CRC_NONE;
  frame_struct->P_MOD =  LIQUID_MODEM_QPSK;

  liquid::ofdm::modulator mod(d_options.M,
                              d_options.cp_len,
                              TAPER_LEN,
                              p,
                              frame_struct);
  mod.set_tx_gain(tx_gain);
  liquid::ofdm::demodulator dem(d_options.M,
                                d_options.cp_len,
                                TAPER_LEN,
                                p,
                                callback,
                                user_data,
                                frame_struct);

  num_frames_detected         = 0;
  num_valid_headers_received  = 0;
  num_valid_bytes_received    = 0;
  num_valid_packets_received  = 0;

  if(SIMULATION) {
    modem = (simulator_data *) malloc (sizeof(simulator_data));
    modem->mod = &mod;
    modem->dem = &dem;
    simulator((void *)modem);
    free(modem);
  } 
  else {
    rx = uhd::usrp::multi_usrp::make(
         uhd::device_addr_t(B210B));
    tx = uhd::usrp::multi_usrp::make(
         uhd::device_addr_t(B210A));
    tx_config = (usrp_config *) malloc 
                (sizeof(usrp_config));
    rx_config = (usrp_config *) malloc 
                (sizeof(usrp_config));

    // assign tx and rx sampling rate freq etc
    tx_config->cent_freq = d_options.cent_freq;
    rx_config->cent_freq = d_options.cent_freq;
    tx_config->samp_rate = d_options.samp_rate;
    rx_config->samp_rate = d_options.samp_rate;
    tx_config->rf_gain   = d_options.txgain;
    rx_config->rf_gain   = d_options.rxgain;
    // TODO Add options for clock and time source
    tx_config->clock_source = CLOCK_SOURCE;
    rx_config->clock_source = CLOCK_SOURCE;
    tx_config->time_source = TIME_SOURCE;
    rx_config->time_source = TIME_SOURCE;

    init_usrp(&tx, tx_config, true);
    init_usrp(&rx, rx_config, false);

    tx_data = (tx_thread_data *) malloc (sizeof(tx_thread_data));
    rx_data = (rx_thread_data *) malloc (sizeof(rx_thread_data));
    tx_data->tx = &tx;
    rx_data->rx = &rx;
    tx_data->mod = &mod;
    rx_data->dem = &dem;
    tx_data->tx_begin = &tx_begin;
    rx_data->rx_begin = &rx_begin;
    tx_data->tx_end = &tx_end;
    rx_data->rx_end = &rx_end;
    tx_data->pid = &pid;
    tx_data->runtime = RUNTIME;
    rx_data->runtime = RUNTIME;
    rx_data->streamer_error = &rx_streamer_error;
    rx_data->rxmd_to_main = &rxmd_from_worker;
    rx_data->online_decoding = online_decoding;
    rx_data->time_for_offline_decoding = &time_for_offline_decoding;
    if(pthread_create(&tx_thread, NULL, tx_worker, (void *)tx_data)){
      std::cout << "Error invoking tx thread\n";
      return 1;
    }
    if(pthread_create(&rx_thread, NULL, rx_worker, (void *)rx_data)){
      std::cout << "Error invoking rx thread\n";
      return 1;
    }
    pthread_join(tx_thread, NULL);
    pthread_join(rx_thread, NULL);
    free(tx_data);
    free(rx_data);
    free(frame_struct);
    free(tx_config);
    free(rx_config);
  }
  printf("    frames transmitted      : %6u\n", pid);
  printf("    frames detected         : %6u\n", num_frames_detected);
  printf("    valid headers           : %6u\n", num_valid_headers_received);
  printf("    valid packets           : %6u\n", num_valid_packets_received);
  printf("    bytes received          : %6u\n", num_valid_bytes_received);
  std::cout << "    tx runtime              : "
            << tx_end - tx_begin
            << std::endl;
  std::cout << "    rx runtime              : "
            << rx_end - rx_begin
            << std::endl;
  printf("    bit rate                : %e bps\n",
         float(num_valid_bytes_received*8)/(rx_end - rx_begin));
  if(!(online_decoding))
    printf("    decoding time           : %3.3f sec\n",
           time_for_offline_decoding);
  if(rx_streamer_error)
  {
    std::cout << "    rx thread exit with     : "
              << rxmd_from_worker.strerror() << "\n";
  }
  return 0;
}
