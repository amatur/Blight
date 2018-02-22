
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include "bbhash.h"
#include "ksl.h"



using namespace std;



#define kmer uint64_t



typedef boomphf::SingleHashFunctor<kmer>  hasher;
typedef boomphf::mphf<  kmer, hasher  > MPHF;




uint64_t nuc2int(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	cout<<"bug"<<c<<"!"<<endl;
	exit(0);
	return 0;
}



uint64_t nuc2intrc(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		case 'T': return 0;
	}
	cout<<"bug"<<c<<"!"<<endl;
	exit(0);
	return 0;
}

string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str){
	return (min(str,revComp(str)));
}


kmer str2num(const string& str){
	kmer res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}


uint32_t xs(uint32_t y){
	y^=(y<<13); y^=(y>>17);y=(y^=(y<<15)); return y;
}


string min_accordingtoXS(const string& seq1,const string& seq2){
	uint32_t u1(str2num(seq1)),u2(str2num(seq2));
	if(xs(u1)<xs(u2)){
		return seq1;
	}
	return seq2;
}



kmer min_accordingtoXS(const kmer& u1,const kmer& u2){
	if(xs(u1)<xs(u2)){
		return u1;
	}
	return u2;
}



kmer rcb(kmer min,uint n){
	kmer res(0);
	kmer offset(1);
	offset<<=(2*n-2);
	for(uint i(0); i<n;++i){
		res+=(3-(min%4))*offset;
		min>>=2;
		offset>>=2;
	}
	return res;
}

void kmer_Set_Light::updateK(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}



kmer min_k (const kmer& k1,const kmer& k2){
	if(k1<=k2){
		return k1;
	}
	return k2;
}



void kmer_Set_Light::updateRCK(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
}



void kmer_Set_Light::create_super_buckets(const string& input_file){
	ifstream inUnitigs(input_file);
	if( not inUnitigs.good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	vector<ofstream*> out_files;
	for(uint i(0);i<number_superbuckets;++i){
		auto out =new ofstream("_out"+to_string(i));
		out_files.push_back(out);
	}
	string ref,useless;
	string super_minimizer,kmer,minimizer,mmer;
	while(not inUnitigs.eof()){
		getline(inUnitigs,useless);
		getline(inUnitigs,ref);
		//~ cout<<ref<<endl;
		//FOREACH UNITIG
		if(not ref.empty() and not useless.empty()){
			super_minimizer=kmer=minimizer=mmer="";
			uint last_position(0);
			//FOREACH KMER
			uint i(0);
			for(;i+k<=ref.size();++i){
				kmer=getCanonical(ref.substr(i,k));
				minimizer=(getCanonical(kmer.substr(0,m)));
				//COMPUTE KMER MINIMIZER
				for(uint j(1);j+m<kmer.size();++j){
					mmer=getCanonical(kmer.substr(j,m));
					minimizer=min_accordingtoXS(minimizer,mmer);
				}
				if(super_minimizer==""){
					super_minimizer=minimizer;
				}else{
					if(super_minimizer!=minimizer){
						*(out_files[xs(str2num(super_minimizer))%number_superbuckets])<<">"+(super_minimizer)+"\n"<<ref.substr(last_position,i-last_position+k-1)<<"\n";
						last_position=i;
						super_minimizer=minimizer;
					}
				}
			}
			if(ref.size()-last_position>k-1){
				//~ cout<<"NADINE"<<endl;
				//~ cout<<ref.substr(last_position)<<endl;
				//~ cout<<minimizer<<endl;
				*(out_files[xs(str2num(super_minimizer))%number_superbuckets])<<">"+(super_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
			}
		}
	}
	for(uint i(0);i<number_superbuckets;++i){
		out_files[i]->close();
		delete(out_files[i]);
	}
}



void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);
	uint number_superbuckets_by_buckets(minimizer_number/number_superbuckets);
	#pragma omp parallel num_threads(8)
	{
		string useless,line;
		#pragma omp for
		for(uint SBC=0;SBC<number_superbuckets;++SBC){
			//~ cin.get();cout<<SBC<<endl;
			ifstream in(input_file+to_string(SBC));
			while(not in.eof()){
				getline(in,useless);
				getline(in,line);
				if(not useless.empty()){
					useless=useless.substr(1);
					uint minimizer(str2num(useless));
					//~ cout<<minimizer<<endl;
					buckets[minimizer]+=line;
					#pragma omp atomic
					total_size+=line.size();
				}
			}
		}
	}
	cout<<"Total size of the partitionned graph: "<<intToString(total_size)<<endl;
}



void kmer_Set_Light::create_mphf(){
	uint64_t max_size(0);
	#pragma omp parallel for num_threads(1)
	for(uint BC=(0);BC<buckets.size();++BC){
		if(buckets[BC].size()<100000000){
		//~ if(buckets[BC].size()==0){
			buckets_size[BC]=0;
			continue;
		}
		auto anchors=new vector<kmer>;
		kmer seq(str2num(buckets[BC].substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
		anchors->push_back(canon);
		for(uint j(0);j+k<buckets[BC].size();++j){
			if(buckets[BC][j+k]=='\n'){
				cout<<"endl"<<endl;
				j+=k;
				seq=(str2num(buckets[BC].substr(j,k))),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
				anchors->push_back(canon);
			}else{
				updateK(seq,buckets[BC][j+k]);
				updateRCK(rcSeq,buckets[BC][j+k]);
				canon=(min_k(seq, rcSeq));
				anchors->push_back(canon);
			}

		}
		//~ cout<<anchors->size()<<endl;
		max_size=max(max_size,anchors->size());
		auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
		kmer_MPHF[BC]= boomphf::mphf<kmer,hasher>(anchors->size(),data_iterator3,8,gammaFactor,false);
		buckets_size[BC]=anchors->size();

		delete anchors;
	}
	cout<<"Largest bucket: "<<intToString(max_size)<<endl;
}



void int_to_bool(uint n_bits_to_encode,uint32_t X, uint32_t pos,vector<bool>& res){
	//~ cout<<"ITB"<<endl;

	for(uint i(0);i<n_bits_to_encode;++i){
		//~ cout<<X<<endl;
		res[i+pos*n_bits_to_encode]=X%2;
		//~ cout<<X%2<<endl;
		//~ cout<<"pos"<<i+pos*n_bits_to_encode<<endl;
		X>>=1;
	}
}


uint32_t bool_to_int(uint n_bits_to_encode,uint pos,const vector<bool>& V){
	uint32_t res(0);
	uint32_t acc(1);
	for(uint i(0);i<n_bits_to_encode;++i){
		//~ cout<<i+pos*n_bits_to_encode<<endl;
		if(V[i+pos*n_bits_to_encode]){
			res+=acc;
		}
		acc<<=1;
	}
	return res;
}



void kmer_Set_Light::fill_positions(){
	//TODO MULTITHREAD OMG
	//~ cout<<"FILL POSITION"<<endl;
	#pragma omp parallel for num_threads(8)
	for(uint BC=(0);BC<buckets.size();++BC){
		if(buckets_size[BC]==0){
			continue;
		}
		uint n_bits_to_encode(log2(buckets_size[BC])+1);
		//~ cout<<n_bits_to_encode<<endl;
		positions[BC].resize(buckets_size[BC]*n_bits_to_encode,0);
		kmer seq(str2num(buckets[BC].substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));

		//~ positions[BC][kmer_MPHF[BC].lookup(canon)].assign(temp.begin(),temp.end());
		//~ for(uint ii(0);ii<positions[BC].size();++ii){
				//~ cout<<positions[BC][ii];
			//~ }
			//~ cout<<endl;
		for(uint j(0);j+k<buckets[BC].size();++j){
			updateK(seq,buckets[BC][j+k]);
			updateRCK(rcSeq,buckets[BC][j+k]);
			canon=(min_k(seq, rcSeq));
			int_to_bool(n_bits_to_encode,j+1,kmer_MPHF[BC].lookup(canon),positions[BC]);
			//~ cout<<"encode"<<j+1<<" "<<kmer_MPHF[BC].lookup(canon)<<endl;
			//~ for(uint ii(0);ii<positions[BC].size();++ii){
				//~ cout<<positions[BC][ii];
			//~ }
			//~ cout<<endl;
		}
	}
}



bool kmer_Set_Light::exists(const string& query2){
	//~ cout<<endl<<"GO QUERY"<<endl;
	string query=getCanonical(query2);
	string minimizer, mmer;
	for(uint j(0);j<query.size();++j){
		minimizer=(getCanonical(query.substr(0,m)));
		//COMPUTE KMER MINIMIZER
		for(uint j(0);j+m<query.size();++j){
			mmer=getCanonical(query.substr(j,m));
			minimizer=min_accordingtoXS(minimizer,mmer);
		}
	}
	uint32_t mini((str2num(minimizer)));
	//~ cout<<kmer_MPHF[mini].lookup(str2num(query))<<endl;
	int32_t hash(kmer_MPHF[mini].lookup(str2num(query)));
	if(hash<0){
		return 0;
	}
	uint pos(bool_to_int( ((uint32_t)log2(buckets_size[mini])+1), hash, positions[mini]));
	//~ uint32_t pos(positions[mini][]);
	//~ cout<<buckets[mini]<<endl;
	if(buckets[mini].substr(pos,k)==query){
		return true;
	}
	if(buckets[mini].substr(pos,k)==revComp(query)){
		return true;
	}
	//~ cout<<buckets[mini].substr(pos,k)<<endl;
	//~ cout<<query<<endl;
	//~ cout<<revComp(query)<<endl;
	//~ cin.get();

	return false;
}



void kmer_Set_Light::multiple_query(const string& query_file){
	ifstream in(query_file);
	//TODO MEGA SLOW QUERY

	uint TP(0),FP(0);
	#pragma omp parallel num_threads(8)
	while(not in.eof()){
		string query,Q;
		#pragma omp critical(dataupdate)
		{
			getline(in,query);
			getline(in,query);
		}
		for(uint i(0);i+k<=query.size();++i){
			Q=query.substr(i,k);
			//~ cout<<"aller"<<endl;
			//~ cout<<Q<<endl;
			//~ cout<<revComp(Q)<<endl;
			if(exists(Q)){
				TP++;
			}else{
				FP++;
			}
		}
		//~ cout<<TP<<" "<<FP<<endl;
	}

	cout<<"Good kmer: "<<TP<<endl;
	cout<<"Erroneous kmers: "<<FP<<endl;
	cout<<"I am glad you are here with me. Here at the end of all things."<<endl;
}









