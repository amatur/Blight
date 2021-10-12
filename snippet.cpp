#include "blight.h"
using namespace std;


string reverseComplement(string base) {
    size_t len = base.length();
    char* out = new char[len + 1];
    out[len] = '\0';
    for (int i = 0; i < len; i++) {
        if (base[i] == 'A') out[len - i - 1] = 'T';
        else if (base[i] == 'C') out[len - i - 1] = 'G';
        else if (base[i] == 'G') out[len - i - 1] = 'C';
        else if (base[i] == 'T') out[len - i - 1] = 'A';
    }
    string outString(out);
    free(out);
    return outString;
}


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
    // string input_kff_outstr_with_count=string(argv[1]);
	// string input_spss_fasta_without_count=string(argv[2]);
	// string output_spss_fasta_with_count=string(argv[4]);
    // string output_spss_fasta_with_count=string(argv[5]);

    int kmer_size = (int) strtol(argv[1], NULL, 10);
    string basename = string(argv[2]);

    //argv1 = k
    //argv2 = basename
    // string input_kff_outstr_with_count="example/m_88095.outstr";
	// string input_spss_fasta_without_count="example/m_88095.fa.essd";
	// string output_spss_fasta_with_count="example/m_88095.instr";
    // string input_minimizer_file="example/m_88095.minimizer";
    // string minimizer("AAATAACACA");
    // string output_position_file = "example/m_88095.pos";
    // string output_maxlen_file = "example/m_88095.maxlen";
    //int kmer_size(32);

    string input_kff_outstr_with_count=basename+".outstr";
	string input_spss_fasta_without_count=basename+".fa.essd";
	string output_spss_fasta_with_count=basename+".instr";
    string input_minimizer_file=basename+".minimizer";
    string output_position_file = basename+".pos";
    string output_maxlen_file = basename+".maxlen";

    //read minimizer
    std::string minimizer;
    std::ifstream infile_minimizer;
    infile_minimizer.open(input_minimizer_file);
    infile_minimizer >> minimizer;
    infile_minimizer.close();


	//SOME VARIABLES TO PLAY WITH
    int core_number(40);
    int minimizer_size(minimizer.size()-1);
    int file_number_exponent(4);
    int subsampling_bits(0);
	

    //INDEX INITIZALIZATION
    //Index Initialization with a given kmer size
    //kmer_Set_Light blight_index_5(kmer_size);
	
    //Index Initialization allowing the use  of multiple thread for faster construction and queries (default is 1)
    //kmer_Set_Light blight_index_2(kmer_size, core_number);

    //Index Initialization with a given minimizer size (default is 10)
    kmer_Set_Light blight_index_5(kmer_size, core_number, minimizer_size);

    //Index Initialization allowing a custom  amount of temporary file (default is 4 for 256 files)
    //kmer_Set_Light blight_index_4(kmer_size, core_number, minimizer_size, file_number_exponent);
    //Blight will create 4^file_number_exponent temporary files (1->4 2>16 3->64 4->256 5->1024)

    //Index Initialization using position subsampling to reduce the memory usage of the index (default is 0)
    //kmer_Set_Light blight_index_3(kmer_size, core_number, minimizer_size, file_number_exponent, subsampling_bits);
    //A value of N will save up to N bits per kmer but each query can lead up to 2^N comparisons to find a kmer in the index



    //INDEX CONSTRUCTION
    //Construct the index from an input FASTA file, temporary files are put in folder wdir that MUST exist (!)
    blight_index_5.construct_index(input_spss_fasta_without_count, "wdir");
    //Blight handles .gz and .lz4 files but do not handle FASTQ or multiline FASTA

    cout<<"The index contains "<<blight_index_5.get_kmer_number()<<" distinct kmers"<<endl;
cout<<"Complete1!"<<endl;
    //We allocate an  integer for each  kmer of the index and set them to 0
    int abundance[blight_index_5.get_kmer_number()]={0};
    cout<<"Complete2!"<<endl;
    vector<int64_t> hash_vector;
    ifstream f_sset;
    f_sset.open (input_kff_outstr_with_count);
    string line;
        int count=0;
    while ( getline (f_sset,line) )
    {
        stringstream ss(line);
        string kmer = "";
        ss >> kmer;
        //cout<<kmer<<endl;
            //We print the identifier of the query performed on the loaded index
            hash_vector=blight_index_5.get_hashes_query(kmer);
            for(int i(0);i<hash_vector.size();++i){
                //cout<<kmer<< " " << hash_vector[i]<<' ';
                int indice(hash_vector[i]);
                //If the kmer is in the index its indice will be a position in the abundance vector
                if(indice!=-1){
                    //We increment the associated counter
        			ss >> count;
                    abundance[indice]=count;
                }
            }
            //cout<<endl;
     }
     f_sset.close();
    ifstream f_spss(input_spss_fasta_without_count);
    ofstream f_out(output_spss_fasta_with_count);
    ofstream f_pos(output_position_file);
     ofstream f_maxlen(output_maxlen_file);
    uint64_t max_len=0;
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
           

            //f_out<<spell<<" "; print actual seq
            //mod: print reverse complement and without mini
            size_t minimizer_pos=spell.find(minimizer);
            if(minimizer_pos == std::string::npos){
                spell = reverseComplement(spell);
                minimizer_pos=spell.find(minimizer);
            }
            if(minimizer_pos == std::string::npos){
                cerr<<"minimizer error"<<endl;
                exit(2);
            }
            
             f_pos<<minimizer_pos<<endl;
            
            if((uint64_t) spell.length()>max_len){
                max_len=(uint64_t) spell.length();
                cout<<max_len<<endl;
            }
            
            //removing minimizer: spell=spell.substr(0, minimizer_pos)+spell.substr(minimizer_pos+minimizer.length(), spell.length()-minimizer_pos-minimizer.length());
            
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
     f_pos.close();
     f_sset.close();

     max_len=max_len-kmer_size+1;
     f_maxlen<<max_len<<endl;
     f_maxlen.close();

     cout<<"Complete3!"<<endl;
     return 0;
}
