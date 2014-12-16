#include "test.h"
#include <string.h>

namespace liquid {
  namespace test {
    copy::copy(unsigned int item_size)
    {
      d_item_size = item_size;
    }

    void copy::work(void *out_buff, void *in_buff, unsigned int num_items)
    {
      memmove(out_buff, in_buff, d_item_size*num_items);
    }
  }
}
