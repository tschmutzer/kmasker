//
//  multi_db.cpp
//  
//
//  Created by Chris Ulpinnis (Uni) on 22.08.16.
//
//

#include <stdio.h>


#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_heap.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/cpp_array.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <stdio.h>
#include <stdlib.h>

namespace err = jellyfish::err;

using jellyfish::file_header;
using jellyfish::RectangularBinaryMatrix;
using jellyfish::mer_dna;
using jellyfish::cpp_array;
using namespace std;
typedef std::unique_ptr<binary_reader>           binary_reader_ptr;
typedef std::unique_ptr<text_reader>             text_reader_ptr;

int IntDivRoundUp(int n, int d) {
    // If n and d are the same sign ...
    if ((n < 0) == (d < 0)) {
        // If n (and d) are negative ...
        if (n < 0) {
            n = -n;
            d = -d;
        }
        // Unsigned division rounds down.  Adding d-1 to n effects a round up.
        return (((unsigned) n) + ((unsigned) d) - 1)/((unsigned) d);
    }
    else {
        return n/d;
    }
}

struct file_info {
    std::ifstream is;
    file_header   header;
    
    file_info(const char* path) :
    is(path),
    header(is)
    { }
};


struct common_info {
    unsigned int            key_len;
    size_t                  max_reprobe_offset;
    size_t                  size;
    unsigned int            out_counter_len;
    std::string             format;
    RectangularBinaryMatrix matrix;
    
    common_info(RectangularBinaryMatrix&& m) : matrix(std::move(m))
    { }
};

common_info read_headers(int argc, char* input_files[], cpp_array<file_info>& files) {
    // Read first file
    files.init(0, input_files[0]);
    if(!files[0].is.good())
        err::die(err::msg() << "Failed to open input file '" << input_files[0] << "'");
    
    file_header& h = files[0].header;
    common_info res(h.matrix());
    res.key_len            = h.key_len();
    res.max_reprobe_offset = h.max_reprobe_offset();
    res.size               = h.size();
    res.format = h.format();
    size_t reprobes[h.max_reprobe() + 1];
    h.get_reprobes(reprobes);
    res.out_counter_len = h.counter_len();
    
    // Other files must match
    for(int i = 1; i < argc; i++) {
        files.init(i, input_files[i]);
        file_header& nh = files[i].header;
        if(!files[i].is.good())
            err::die(err::msg() << "Failed to open input file '" << input_files[i] << "'");
        if(res.format != nh.format())
            err::die(err::msg() << "Can't compare files with different formats (" << res.format << ", " << nh.format() << ")");
        if(res.key_len != nh.key_len())
            err::die(err::msg() << "Can't compare hashes of different key lengths (" << res.key_len << ", " << nh.key_len() << ")");
        if(res.max_reprobe_offset != nh.max_reprobe_offset())
            err::die("Can't compare hashes with different reprobing strategies");
        if(res.size != nh.size())
            err::die(err::msg() << "Can't compare hash with different size (" << res.size << ", " << nh.size() << ")");
        if(res.matrix != nh.matrix())
            err::die("Can't compare hash with different hash function");
    }
    
    return res;
}

/*template<typename reader_type>
void output_counts(cpp_array<file_info>& files) {
    cpp_array<reader_type> readers(files.size());
    typedef jellyfish::mer_heap::heap<mer_dna, reader_type> heap_type;
    typedef typename heap_type::const_item_t heap_item;
    heap_type heap(files.size());
    
    // Prime heap
    for(size_t i = 0; i < files.size(); ++i) {
        readers.init(i, files[i].is, &files[i].header);
        if(readers[i].next())
            heap.push(readers[i]);
    }
    
    heap_item          head      = heap.head();
    mer_dna            key;
    const int          num_files = files.size();
    const reader_type* base      = &readers[0];
    uint64_t           counts[num_files];
    while(heap.is_not_empty()) {
        key = head->key_;
        memset(counts, '\0', sizeof(uint64_t) * num_files);
        do {
            counts[head->it_ - base] = head->val_;
            heap.pop();
            if(head->it_->next())
                heap.push(*head->it_);
            head = heap.head();
        } while(head->key_ == key && heap.is_not_empty());
        std::cout << key;
        for(int i = 0; i < num_files; ++i)
            std::cout << " " << counts[i];
        std::cout << "\n";
    }
}*/

template<typename reader_type, typename PathIterator, typename Database>
void compare_counts(cpp_array<file_info>& PathIterator file_begin, PathIterator file_end ) {
    cpp_array<reader_type> readers(files.size());
    typedef jellyfish::mer_heap::heap<mer_dna, reader_type> heap_type;
    typedef typename heap_type::const_item_t heap_item;
    heap_type heap(files.size());
    jellyfish::stream_manager<PathIterator> streams(file_begin, file_end);
    sequence_parser                         parser(4, 100, 1, streams);

    
    // Prime heap
    for(size_t i = 0; i < files.size(); ++i) {
        readers.init(i, files[i].is, &files[i].header);
        if(readers[i].next())
            heap.push(readers[i]);
    }
    
    heap_item          head      = heap.head();
    mer_dna            key;
    const int          num_files = files.size();
    const reader_type* base      = &readers[0];
    uint64_t           counts[num_files];
    while(heap.is_not_empty()) {
        key = head->key_;
        memset(counts, '\0', sizeof(uint64_t) * num_files);
        do {
            counts[head->it_ - base] = head->val_;
            heap.pop();
            if(head->it_->next())
                heap.push(*head->it_);
            head = heap.head();
        } while(head->key_ == key && heap.is_not_empty());
        std::cout << key;
        for(int i = 0; i < num_files; ++i)
            std::cout << " " << counts[i];
        std::cout << "\n";
    }
}

int main(int argc, char *argv[])
{
    // Check number of input files
    if(argc < 3)
        err::die(err::msg() << "Usage: " << argv[0] << " FASTA DB1 DB2 ...");
    char* fasta[1];
    fasta[0]=fasta[1];
    
    // Read the header of each input files and do sanity checks.
    cpp_array<file_info> files(argc - 2);
    common_info cinfo = read_headers(argc - 2, argv + 2, files);
    mer_dna::k(cinfo.key_len / 2);
    
    if(cinfo.format == binary_dumper::format)
        output_counts<binary_reader>(files, fasta, fasta+1);
    else if(cinfo.format == text_dumper::format)
        output_counts<text_reader>(files, fasta, fasta+1);
    else
        err::die(err::msg() << "Format '" << cinfo.format << "' not supported\n");
    
    return 0;
}
