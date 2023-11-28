#include <iostream>
#include <vector>
#include <fstream>
#include <ios>
#include <cstring>
#include <cstdlib>
#include <execution>

// 37 bytes before the data begins; this is the FAISS
// header. 
// https://github.com/facebookresearch/faiss/blob/main/faiss/impl/io.cpp 
struct flat_header {
  uint32_t fourcc;  // 'IFxI'
  int32_t  dim;     // Vector dimensionality
  int64_t  ntotal;  // Number of vectors stored
  int64_t  dummy_a; // A dummy value, 1 << 20
  int64_t  dummy_b; // A dummy value, 1 << 20
  bool     trained; // Is the model trained? Meaningless here.
  uint32_t   metric;  // Metric type; actually an enum in FAISS, meaningless here.

  // Loads into self from an istream
  void load(std::istream& in) {
    in.read(reinterpret_cast<char *>(&fourcc), sizeof(uint32_t));
    in.read(reinterpret_cast<char *>(&dim), sizeof(int32_t));
    in.read(reinterpret_cast<char *>(&ntotal), sizeof(int64_t));
    in.read(reinterpret_cast<char *>(&dummy_a), sizeof(int64_t));
    in.read(reinterpret_cast<char *>(&dummy_b), sizeof(int64_t));
    in.read(reinterpret_cast<char *>(&trained), sizeof(bool));
    in.read(reinterpret_cast<char *>(&metric), sizeof(uint32_t));
  }

  // Back to disk
  void write(std::ostream& out) {
    out.write(reinterpret_cast<const char *>(&fourcc), sizeof(uint32_t));
    out.write(reinterpret_cast<char *>(&dim), sizeof(int32_t));
    out.write(reinterpret_cast<char *>(&ntotal), sizeof(int64_t));
    out.write(reinterpret_cast<char *>(&dummy_a), sizeof(int64_t));
    out.write(reinterpret_cast<char *>(&dummy_b), sizeof(int64_t));
    out.write(reinterpret_cast<char *>(&trained), sizeof(bool));
    out.write(reinterpret_cast<char *>(&metric), sizeof(uint32_t));
  }

  // Some basic info from the header
  void info() const {
    std::cout << "Dim    = " << dim << "\n"
              << "Ntotal = " << ntotal << "\n";
  }

};

// Everything we need to store/recover vectors without FAISS
class vector_data_32 {

  public:
    vector_data_32(size_t dim, size_t ntotal) : m_dimensions(dim), m_num_vectors(ntotal) {}

    // Loads into self from an istream and checks vector length
    void load(std::istream& in) {
      size_t vals_to_read;
      in.read(reinterpret_cast<char *>(&vals_to_read), sizeof(size_t));
      m_codes.resize(vals_to_read); // This is in terms of bytes
      in.read(reinterpret_cast<char *>(&m_codes[0]), vals_to_read * sizeof(float));
      if (!in) {
        std::cout << "Error: only " << in.gcount() << " bytes could be read\n";
        std::exit(-1);
      }
    }

    // New format
    void write(std::ostream& out) const {
      size_t elements = m_dimensions * m_num_vectors;
      out.write(reinterpret_cast<const char *>(&m_dimensions), sizeof(size_t));
      out.write(reinterpret_cast<const char *>(&m_num_vectors), sizeof(size_t));
      out.write(reinterpret_cast<const char *>(&m_codes[0]), elements * sizeof(uint32_t)); // Assume 32 bits
    }

    void sort() {
      std::sort(std::execution::par_unseq, m_codes.begin(), m_codes.end());
    }

  private:
    size_t               m_dimensions;  // How large are the strides?
    size_t               m_num_vectors; // Where does the data end?
    std::vector<float>   m_codes;       // The data itself 
};


// Assume 4-byte (floats)
const size_t UNIT_BYTES = 4;

int main(int argc, char **argv) {

  if (argc != 3) {
    std::cerr << "Usage " << argv[0] << "<path_to_flat_FAISS_index> <out_index>\n";
    return -1;
  }

  // Load the FAISS flat index
  std::ifstream ifs(argv[1], std::ios::binary);
  flat_header fh;
  fh.load(ifs);
  vector_data_32 idx(fh.dim, fh.ntotal);
  idx.load(ifs);

  // Sort the numbers for quantization later
  idx.sort();

  // Dump the data as an `sidx` file
  std::ofstream ofs(argv[2], std::ios::binary);
  idx.write(ofs);
 
}
