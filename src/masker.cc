/*  Modified by Chris Ulpinnis for cMasker project to mask fasta files with jellyfish.
        
    Based on jellyfishs example code (query_per_sequence.cc)
 
    Jellyfish/cMasker are free software: you can redistribute it and/or modify
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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fstream>
#include <libgen.h>
#include <string>
using namespace std;
//#include "sequence_mers.hpp"
namespace err = jellyfish::err;
using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;

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
} //http://stackoverflow.com/questions/17005364/dividing-two-integers-and-rounding-up-the-result-without-using-floating-pont

typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;

template<typename PathIterator, typename Database>
void query_from_sequence(PathIterator file_begin, PathIterator file_end, const Database& db,
                         bool canonical, bool occfile, int rt, int norm) {
    jellyfish::stream_manager<PathIterator> streams(file_begin, file_end);
    sequence_parser                         parser(4, 100, 1, streams);
    ofstream occfilestream;
    ofstream occnormfilestream;
    char* path = *file_begin;
    string file = basename(path);
    string dir = dirname(path);
    //cout << dir << "/" << file << "\n";
    string occnormoutname = dir + "/" + file + "_N" + to_string(norm) +"normalized.occ";
    string occoutname = dir + "/" + file +".occ";
    if(occfile == true){
        occfilestream.open(occoutname);
        occnormfilestream.open(occnormoutname);
    }
    ofstream fastaout;
    //cout << "input was: " << string(dirname(*file_begin)) << "/" << string(basename(*file_begin)) << "\n";
    //cout << string(*file_begin) << "\n";
    //cout << string(*file_begin) << "\n";
    string fastaoutname = dir + "/freakmasked_RT" + to_string(rt) + "." + file;
    cout << "Out is: " << fastaoutname << "\n";
    fastaout.open(fastaoutname);
  //sequence_mers                           mers(canonical);
  //const sequence_mers                     mers_end(canonical);
  while(true) {
    sequence_parser::job j(parser);
    if(j.is_empty()) break;
    for(size_t i = 0; i < j->nb_filled; ++i) {
      fastaout << ">" << j->data[i].header << "\n";
        if(occfile) {
            occfilestream << ">" << j->data[i].header << "\n";
            occnormfilestream << ">" << j->data[i].header << "\n";
        }
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
            //std::cout << "-1";
            if (occfile) {
                occfilestream << "0";
                occnormfilestream << "0";
            }
            fastaout << j->data[i].seq.substr(0,1);
        }
        else {
            //std::cout << db.check(mer) << " : " << mer.to_str()<< " : " << mershift;
            mershift = 0; //necessary because mershift is std::string::npos after the above condition
            if(occfile){
                occnormfilestream << IntDivRoundUp(db.check(mer),norm);
                occfilestream << db.check(mer);
            }
            if(IntDivRoundUp(db.check(mer),norm) > rt) {
                fastaout << "X";
            }
            else{
                fastaout << j->data[i].seq.substr(0,1);
            }
        }
        //now calculate how many times the mer must be shifted
        int k = 2;
        bool whitespace = true;
        string* seq = &j->data[i].seq;
        int length = seq->length();
        //for (char& c: j->data[i].seq.substr(mer.k()), s: j->data[i].seq.substr(1)) {
        for(int cc = mer.k(),cs=1 ; cc < length; ++cc, ++cs) {
        //for(char& c: j->data[i].seq.substr(1)) {
            char c =  seq->at(cc);
            char s =  seq->at(cs);
            if (whitespace == true) {
                if (occfile) {
                    occfilestream << " ";
                    occnormfilestream << " ";
                }
            }
            if (c == 'N') { //got a new N, have to shift mer.k()-1 times from now
                mershift = mer.k() - 1;
                mer.shift_left('A');
                /*std::cout << "-1" << " : " << mer.to_str() << " : " << mershift;
                std::cout << "\n";*/
                //std::cout << "-1";
                if(occfile) {
                    occfilestream << "0";
                    occnormfilestream << "0";
                }
                fastaout << s;

            }
            else {
                if(mershift > 0) { //still needs shifting because there is an N somewhere in the kmer
                    mershift--;
                    mer.shift_left(c);
                    /*std::cout << "-1" << " : " << mer.to_str()<< " : " << mershift;
                    std::cout << "\n";*/
                    //std::cout << "-1";
                    if(occfile) {
                        occfilestream << "0";
                        occnormfilestream << "0";
                    }
                    fastaout << s;
                }
                else { //allright, no Ns in the kmer
                    mer.shift_left(c);
                    /*std::cout << db.check(mer) << " : " << mer.to_str()<< " : " << mershift;
                    std::cout << "\n";*/
                    if (occfile) {
                        occnormfilestream << IntDivRoundUp(db.check(mer),norm);
                        occfilestream << db.check(mer);
                    }
                    if (IntDivRoundUp(db.check(mer),norm) > rt) {
                        fastaout << "X";
                    }
                    else {
                        fastaout << s;
                    }

                }
                
            }
            k++;
            whitespace = true; //we need a whitespace in the next beginning of the loop
            if(k == 26) { //check if we have to make a line break
                if (occfile) {
                    occfilestream << "\n";
                    occnormfilestream << "\n";
                }
                whitespace=false;
                k=1;
            }
        }
        for (int a = 1; a < mer.k(); a++, k++) { //print last bases which have no kmer
            if(whitespace == false) {
                if (occfile) {
                    occfilestream << "0";
                    occnormfilestream << "0";
                }
                whitespace = true;
            }
            else {
                if (occfile) {
                    occfilestream << " " << "0";
                    occnormfilestream << " " << "0";
                }
            }
            if(k == 26) {
                if (occfile) {
                    occfilestream << "\n";
                    occnormfilestream << "\n";
                }
                whitespace=false;
                k=1;
            }
        }
        //now print the bases in fasta
        fastaout << j->data[i].seq.substr(j->data[i].seq.length() - (mer.k() - 1));
        
        if (occfile) {
            occfilestream << "\n";
            occnormfilestream << "\n";
        }
        fastaout << "\n";
     }
  }
    fastaout.close();
    if (occfile) {
        occfilestream.close();
        occnormfilestream.close();
    }
}

int main(int argc, char *argv[])
{
    bool occflag = false;
    int norm = 1;
    int rt = 40;
    int opt;
    char* fasta = NULL;
    char* jellydb = NULL;
    int index;
    
    
    opterr = 0;
    
    while ((opt = getopt(argc, argv, "of:j:hn:r:")) != EOF) {
        switch (opt) {
            case 'o':
                occflag = true;
                break;
            case 'f':
                fasta = optarg;
                cout << "Input is: " << fasta << "\n";
		break;
            case 'j':
                jellydb = optarg;
                break;
            case 'n':
                norm = atoi(optarg);
                break;
            case 'r':
                rt = atoi(optarg);
                break;
            case 'h':
                cout << "Usage: " << argv[0] << "\n\t-h\t Shows this help\n\t-f\tFASTA Input\n\t-j\tJellfish Database\n\t-o\tCreate OCC output\n\t-n\tNormalize Value\n\t-r\tRT Value for masking threshold\n";
                return 1;
                break;
            case '?':
                if (optopt == 'f' || optopt == 'j' || optopt == 'n' || optopt == 'r')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",
                             optopt);
                return 1;
            default:
                abort ();
        }
    }
    
    for (index = optind; index < argc; index++)
        printf ("Non-option argument %s\n", argv[index]);
    
  /*if(argc < 3)
    err::die(err::msg() << "Usage: " << argv[0] << "db.jf file.fa [...]");*/
    char* file[1]; //I have no idea how this stream iterator thing works in jellyfish. For this reason I just put the input file into an array, which seems to be strange...
    file[0]=fasta;
    std::ifstream in(jellydb, std::ios::in|std::ios::binary);
  jellyfish::file_header header(in);
  if(!in.good())
    err::die(err::msg() << "Failed to parse header of file '" << jellydb << "'");
  mer_dna::k(header.key_len() / 2);
  if(header.format() == "bloomcounter") {
    jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
    mer_dna_bloom_counter filter(header.size(), header.nb_hashes(), in, fns);
    if(!in.good())
      err::die("Bloom filter file is truncated");
    in.close();
    query_from_sequence(file, file+1, filter, header.canonical(),occflag, rt, norm);
  } else if(header.format() == binary_dumper::format) {
    jellyfish::mapped_file binary_map(jellydb);
    binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
                    header.size() - 1, binary_map.length() - header.offset());
    query_from_sequence(file, file+1 , bq, header.canonical(), occflag ,rt, norm);
  } else {
    err::die(err::msg() << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list.");
  }

  return 0;
}
