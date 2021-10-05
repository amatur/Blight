#include "blight.h"
using namespace std;

vector<string> string_split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;
    while (getline (ss, item, delim)) {
        result.push_back (item);
    }
    return result;
}

int main(int argc, char** argv) {
    string input_kff_outstr_with_count=string(argv[1]);
	string input_spss_fasta_without_count=string(argv[2]);
	string output_spss_fasta_with_count=string(argv[3]);



    //SOME VARIABLES TO PLAY WITH
    //int kmer_size(32);
	int kmer_size = strtol(argv[4], NULL, 10);
    int core_number(31);
    int minimizer_size(10);
    int file_number_exponent(4);
    int subsampling_bits(0);

    //INDEX INITIZALIZATION
    //Index Initialization with a given kmer size
    kmer_Set_Light blight_index_1(kmer_size);

    //Index Initialization allowing the use  of multiple thread for faster construction and queries (default is 1)
    kmer_Set_Light blight_index_2(kmer_size, core_number);

    //Index Initialization with a given minimizer size (default is 10)
    kmer_Set_Light blight_index_5(kmer_size, core_number, minimizer_size);

    //Index Initialization allowing a custom  amount of temporary file (default is 4 for 256 files)
    kmer_Set_Light blight_index_4(kmer_size, core_number, minimizer_size, file_number_exponent);
    //Blight will create 4^file_number_exponent temporary files (1->4 2>16 3->64 4->256 5->1024)

    //Index Initialization using position subsampling to reduce the memory usage of the index (default is 0)
    kmer_Set_Light blight_index_3(kmer_size, core_number, minimizer_size, file_number_exponent, subsampling_bits);
    //A value of N will save up to N bits per kmer but each query can lead up to 2^N comparisons to find a kmer in the index



    //INDEX CONSTRUCTION
    //Construct the index from an input FASTA file, temporary files are put in folder wdir that MUST exist (!)
    blight_index_5.construct_index(kff_fasta_with_count, "wdir");
    //Blight handles .gz and .lz4 files but do not handle FASTQ or multiline FASTA

    cout<<"The index contains "<<blight_index_5.get_kmer_number()<<" distinct kmers"<<endl;

    //We allocate an  integer for each  kmer of the index and set them to 0
    int abundance[blight_index_5.get_kmer_number()]={0};

    vector<int64_t> hash_vector;
    ifstream f_sset;
    f_sset.open (kff_fasta_with_count);
    string line;
    while ( getline (f_sset,line) )
    {
        stringstream ss(line);
        string kmer = "";
        ss >> kmer;
        int count=0;
        ss >> count;
            //We print the identifier of the query performed on the loaded index
            hash_vector=blight_index_5.get_hashes_query(kmer);
            for(int i(0);i<hash_vector.size();++i){
                //cout<<kmer<< " " << hash_vector[i]<<' ';
                int indice(hash_vector[i]);
                //If the kmer is in the index its indice will be a position in the abundance vector
                if(indice!=-1){
                    //We increment the associated counter
                    abundance[indice]=count;
                }
            }
            //cout<<endl;
     }
     f_sset.close();

    ifstream f_spss(spss_fasta_without_count);
    ofstream f_out("output_spss_fasta_with_count");
    while ( getline (f_spss,line) )
    {
        stringstream ss(line);
        string spell = "";
        ss >> spell;
        if(spell.at(0)=='>'){
            continue;
        }

            //We print the identifier of the query performed on the loaded index
            hash_vector=blight_index_5.get_hashes_query(spell);
            f_out<<spell<<" ";
            for(int i(0);i<hash_vector.size();++i){
                //cout<<kmer<< " " << hash_vector[i]<<' ';
                 int indice(hash_vector[i]);
                //If the kmer is in the index its indice will be a position in the abundance vector
                if(indice!=-1){
                    if(i!=0){
                        f_out<< ",";
                    }
                    f_out<< abundance[indice];
                }
            }
            f_out<<endl;

     }
     f_out.close();
     f_sset.close();
     cout<<"Complete!"<<endl;
     return 0;
}
