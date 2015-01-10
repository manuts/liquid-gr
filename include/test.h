#ifndef TEST_H
#define TEST_H

namespace liquid {
  namespace test {
    class copy
    {
      public:
        copy(unsigned int item_size);
        void work(void *out_buff, void *in_buff, unsigned int num_items);
      private:
        unsigned int d_item_size;
    };
  }
}

#endif
