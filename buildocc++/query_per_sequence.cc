/*  Modified by Chris Ulpinnis for KMasker project to generate occ files.
        
    This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <algorithm>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/jellyfish.hpp>
//#include "sequence_mers.hpp"

namespace err = jellyfish::err;

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;

template<typename PathIterator, typename Database>
void query_from_sequence(PathIterator file_begin, PathIterator file_end, const Database& db,
                         bool canonical) {
  jellyfish::stream_manager<PathIterator> streams(file_begin, file_end);
  sequence_parser                         parser(4, 100, 1, streams);
  //sequence_mers                           mers(canonical);
  //const sequence_mers                     mers_end(canonical);
  while(true) {
    sequence_parser::job j(parser);
    if(j.is_empty()) break;
    for(size_t i = 0; i < j->nb_filled; ++i) {
      std::cout << ">" << j->data[i].header << "\n";
        /* What we trying to archive with the code modification is:
         -get the first kmer of the sequence string
         -then get the next char from the sequence string by using a pointer to the chars of that string
         -using the shift operation of the sequence-mer-class to append this char */
        //mers = j->data[i].seq;
        //lets start with the first kmer
        int mershift = 0; //how often we have to shift to get rid of the N
        jellyfish::mer_dna mer;
        std::string firstmer = j->data[i].seq.substr(0,mer.k());
        std::replace (firstmer.begin(), firstmer.end(), 'N', 'A'); //replace all N with A for jellyfish::mer_dna
        mer = firstmer; //use the string to initalize the mer object
        //put the position of the right-most N into mershift and check if there is an N at all
        if((mershift = j->data[i].seq.substr(0,mer.k()).find_last_of('N')) != std::string::npos) {
            /*std::cout << "-1" << " : " << mer.to_str()<< " : " << mershift;
            std::cout << "\n";*/
            std::cout << "-1";
        }
        else {
            //std::cout << db.check(mer) << " : " << mer.to_str()<< " : " << mershift;
            mershift = 0; //necessary because mershift is std::string::npos after the above condition
            std::cout << db.check(mer);
        }
        //now calculate how many times the mer must be shifted
        int k = 1;
        bool whitespace = true;
        for (char& c: j->data[i].seq.substr(mer.k())) {
            if (whitespace == true) {
                std::cout << " ";
            }
            if (c == 'N') { //got a new N, have to shift mer.k()-1 times from now
                mershift = mer.k() - 1;
                mer.shift_left('A');
                /*std::cout << "-1" << " : " << mer.to_str() << " : " << mershift;
                std::cout << "\n";*/
                std::cout << "-1";

            }
            else {
                if(mershift > 0) { //still needs shifting because there is an N somewhere in the kmer
                    mershift--;
                    mer.shift_left(c);
                    /*std::cout << "-1" << " : " << mer.to_str()<< " : " << mershift;
                    std::cout << "\n";*/
                    std::cout << "-1";

                }
                else { //allright, no Ns in the kmer
                    mer.shift_left(c);
                    /*std::cout << db.check(mer) << " : " << mer.to_str()<< " : " << mershift;
                    std::cout << "\n";*/
                    std::cout << db.check(mer);

                }
                
            }
            k++;
            whitespace = true; //we need a whitespace in the next beginning of the loop
            if(k == 40) { //check if we have to make a line break
                std::cout << "\n";
                whitespace=false;
                k=1;
            }
        }
        for (int a = 1; a < mer.k(); a++, k++) { //print last bases which have no kmer
            if(whitespace == false) {
                std::cout << "0";
                whitespace = true;
            }
            else {
                std::cout << " " << "0";
            }
            if(k == 40) {
                std::cout << "\n";
                whitespace=false;
                k=1;
            }
        }
        std::cout << "\n";
     }
  }
}

int main(int argc, char *argv[])
{
  if(argc < 3)
    err::die(err::msg() << "Usage: " << argv[0] << "db.jf file.fa [...]");

  std::ifstream in(argv[1], std::ios::in|std::ios::binary);
  jellyfish::file_header header(in);
  if(!in.good())
    err::die(err::msg() << "Failed to parse header of file '" << argv[1] << "'");
  mer_dna::k(header.key_len() / 2);
  if(header.format() == "bloomcounter") {
    jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
    mer_dna_bloom_counter filter(header.size(), header.nb_hashes(), in, fns);
    if(!in.good())
      err::die("Bloom filter file is truncated");
    in.close();
    query_from_sequence(argv + 2, argv + argc, filter, header.canonical());
  } else if(header.format() == binary_dumper::format) {
    jellyfish::mapped_file binary_map(argv[1]);
    binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
                    header.size() - 1, binary_map.length() - header.offset());
    query_from_sequence(argv + 2, argv + argc, bq, header.canonical());
  } else {
    err::die(err::msg() << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list.");
  }

  return 0;
}
