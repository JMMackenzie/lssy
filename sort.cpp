#include <iostream>
#include <vector>
#include <fstream>
#include <ios>
#include <cstring>
#include <cstdlib>
#include <execution>

/* The fourcc code for indexes.
// Flat gets "IFxI"
// Stolen from FAISS: https://github.com/facebookresearch/faiss/blob/main/faiss/impl/io.cpp 
static uint32_t _fourcc (const char sx[4]) {
    assert(4 == strlen(sx));
    const unsigned char *x = (unsigned char*)sx;
    return x[0] | x[1] << 8 | x[2] << 16 | x[3] << 24;
}
*/

//https://stackoverflow.com/questions/15685181/how-to-get-the-sign-mantissa-and-exponent-of-a-floating-point-number
typedef union {
  float f;
  struct {
    uint32_t mantissa : 23;
    uint32_t exponent : 8;
    uint32_t sign : 1;
  } parts;
} float_cast_32;


// Prunes the least significant b bits from a 32
uint32_t prune_lsb(uint32_t value, uint32_t b) {
  uint32_t mask = (-1 << b);
  return value & mask;
}


// 37 bytes before the data begins
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

  void info() const {
    std::cout << "Dim    = " << dim << "\n"
              << "Ntotal = " << ntotal << "\n";
  }

};

// Everything we need to store/recover vectors
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

    void truncate_bits(uint32_t bits) {
      
      // Work through each float, build streams
      for (size_t i = 0; i < m_codes.size(); ++i) {
        float_cast_32 f = { .f = m_codes[i] };
        //std::cerr << f.parts.mantissa << "\n";
        f.parts.mantissa = prune_lsb(f.parts.mantissa, bits);
        //std::cerr << "Mantissa now " << f.parts.mantissa << "\n";
        m_codes[i] = f.f;
      }
    }

    void sort() {
      std::sort(std::execution::par_unseq, m_codes.begin(), m_codes.end());
    }

    // New format
    void write(std::ostream& out) const {
      size_t elements = m_dimensions * m_num_vectors;
      out.write(reinterpret_cast<const char *>(&m_dimensions), sizeof(size_t));
      out.write(reinterpret_cast<const char *>(&m_num_vectors), sizeof(size_t));
      out.write(reinterpret_cast<const char *>(&m_codes[0]), elements * sizeof(uint32_t)); // Assume 32 bits
    }

    void peel_and_write(std::string outname, uint32_t bits) const {

      // Stream of 32s with one active bit
      std::vector<uint32_t> sign_stream;
      sign_stream.reserve(m_codes.size());

      // Stream of 32s with 8 active bits
      std::vector<uint32_t> expo_stream;
      expo_stream.reserve(m_codes.size());

      // Stream of 32s with 23 active bits
      std::vector<uint32_t> mant_stream;
      mant_stream.reserve(m_codes.size());

      // Work through each float, build streams
      for (size_t i = 0; i < m_codes.size(); ++i) {
        float_cast_32 f = { .f = m_codes[i] };
        //std::cerr << f.f << " " << f.parts.sign << " " << f.parts.exponent << " " << f.parts.mantissa << "\n";
        //sign_stream.push_back(f.parts.sign);
        //expo_stream.push_back(f.parts.exponent);
        mant_stream.push_back(prune_lsb(f.parts.mantissa, bits));
      }

      //std::ofstream out_sign(outname + ".sign", std::ios::binary);
      //std::ofstream out_expo(outname + ".exponent", std::ios::binary);
      std::ofstream out_mant(outname + ".mantissa." + std::to_string(bits), std::ios::binary);
   
      //out_sign.write(reinterpret_cast<char *>(&sign_stream[0]), sizeof(int32_t) * sign_stream.size());
      //out_expo.write(reinterpret_cast<char *>(&expo_stream[0]), sizeof(int32_t) * expo_stream.size());
      out_mant.write(reinterpret_cast<char *>(&mant_stream[0]), sizeof(int32_t) * mant_stream.size());

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

  idx.sort();
  std::ofstream ofs(argv[2], std::ios::binary);
  //fh.write(ofs);
  idx.write(ofs);
 
}
