#include "utils.h"
#include "xor.h"
#include <fstream>
#include <cstring>
#include <unordered_set>



using namespace std;


// cmdline: ./bmean [<genome.fasta>] [--xxhash]
int main(int argc, char ** argv){
	srand (time(NULL));
	uint Nelement(1*1000*100);
	vector<uint64_t> originalkmers(Nelement);
	uint k(31);
	string seq;
	uint P(20);
	for(uint i(0);i<Nelement;++i){
		seq = randomSeq(k);
		originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
	}
	unordered_set<uint64_t> set;
	for(uint i(0);i<Nelement;++i){
		set.insert(originalkmers[i]);
	}
	double gamma(0.5);
	uint count(3);
	int hashmode = 0;
	uint nBitUsed(0);
	vector<uint64_t> kmers=originalkmers;
	vector<vector<bool>> approxSet;
	vector<vector<bool>> approxSet2;
	
	for(uint L(0);kmers.size()>1000;++L){
		uint casesNumber(kmers.size()*gamma+1);
		vector<uint8_t> counting(casesNumber,0);
		vector<uint8_t> counting2(casesNumber,0);
		for(uint i(0);i<kmers.size();++i){
			uint64_t h(iterHash64(kmers[i],2*L,hashmode)%casesNumber);
			++counting[h];
			uint64_t h2(iterHash64(kmers[i],2*L+1,hashmode)%casesNumber);
			++counting2[h2];
		}
		vector<uint64_t> unplaced;
		vector<bool> bitSet(casesNumber,false);
		vector<bool> bitSet2(casesNumber,false);
		for(uint i(0);i<kmers.size();++i){
			uint64_t h(iterHash64(kmers[i],2*L,hashmode)%casesNumber);
			uint64_t h2(iterHash64(kmers[i],2*L+1,hashmode)%casesNumber);
			if(counting[h]<count or counting2[h2]<count){
			//~ if(counting[h]!=1){
				unplaced.push_back(kmers[i]);
			}else{
				bitSet[h]=true;
				bitSet2[h2]=true;
			}
		}
		kmers=unplaced;
		nBitUsed+=2*casesNumber;
		approxSet.push_back(bitSet);
		approxSet2.push_back(bitSet2);
		P=L;
	}
	
	
	uint FP(0);
	for(uint test(0);test<Nelement;){
		seq = randomSeq(k);
		uint64_t kmer(getRepresent(seq2intStranded(seq),k));
		if(set.count(kmer)==0){
			++test;
			for(uint L(0);L<P;++L){
				uint64_t h(iterHash64(kmer,2*L,hashmode)%approxSet[L].size());
				uint64_t h2(iterHash64(kmer,2*L+1,hashmode)%approxSet2[L].size());
				if(approxSet[L][h] and approxSet2[L][h2]){
					++FP;
					goto lol2;
				}
			}
		}
		lol2:;
		//~ cout<<"testé"<<test<<endl;
	}
	
	uint TP(0),FN(0);
	for(uint i(0);i<Nelement;++i){
		uint64_t kmer(originalkmers[i]);
		for(uint L(0);L<P;++L){
			if(approxSet[L].size()!=0){
				uint64_t h(iterHash64(kmer,2*L,hashmode)%approxSet[L].size());
				uint64_t h2(iterHash64(kmer,2*L+1,hashmode)%approxSet2[L].size());
				if(approxSet[L][h] and approxSet2[L][h2]){
					++TP;
					goto lol;
				}
			}else{
				goto lol;
			}
			
		}
		++FN;
		lol:;
	}
	cout<<"nummber level: "<<P<<endl;
	cout<<"unplaced element rate: "<<(double)kmers.size()/Nelement<<endl;
	cout<<"bit used: "<<nBitUsed<<endl;
	cout<<"bit/element: "<<((double)nBitUsed/Nelement)<<endl;
	cout<<"FP rate: "<<((double)FP/Nelement)<<endl;
	//~ cout<<FP<<" "<<Nelement<<endl;
	cout<<"FN: "<<FN-kmers.size()<<endl;
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
    //~ double unplaced_tolerance(0.02); // (corresponds to 2%); used to be 0.1%
	//~ // double gamma(2);
	//~ // uint Ncase(gamma*Nelement);
	
	//~ double minbitsbyelts(5);
	//~ // cout<<"gamma :"<<gamma<<endl;
	//~ vector<uint64_t> originalkmers(Nelement);
    //~ string genome = "";
    //~ if (argc>1 && argv[1][0] != '-') { // whether a genome is passed as the 1st argument
        //~ genome = loadGenome(string(argv[1]));
        //~ if (genome.size() < Nelement-k+1) { cout << "not enough genomic kmers in " << string(argv[1]) << endl; exit(1);}
    //~ }
    //~ int hashmode = 0;
    //~ for (int i = 1; i < argc; i++)
        //~ if (strcmp(argv[i],"--xxhash")==0) { cout << "using xxhash" << endl; hashmode = 1; }
	//~ for(uint i(0);i<Nelement;++i){
        //~ string seq = genome.size() > 0 ? genomeKmer(genome,i,k) : randomSeq(k);
		//~ originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
	//~ }

	//~ // for(uint H(2);H<257;H*=2){
	//~ // 	vector<uint64_t> kmers=originalkmers;
	//~ // 	vector<vector<uint64_t>> map(Ncase);
	//~ // 	for(uint n(0);n<H-1;++n){
	//~ // 		for(uint i(0);i<kmers.size();++i){
	//~ // 			uint64_t h(iterHash64(kmers[i],n)%Ncase);
	//~ // 			map[h].push_back(kmers[i]);
	//~ // 		}
	//~ // 		kmers={};
	//~ // 		for(uint i(0);i<Ncase;++i){
	//~ // 			if(map[i].size()>1){
	//~ // 				for(uint j(0);j<map[i].size();++j){
	//~ // 					kmers.push_back(map[i][j]);
	//~ // 				}
	//~ // 				map[i]={};
	//~ // 			}
	//~ // 		}
	//~ // 	}
	//~ // 	cout<<"Number of hash function:    "<<H-1<<endl;
	//~ // 	cout<<"Percent Unplaced element:  "<<(double)100*kmers.size()/Nelement<<" %"<<endl<<endl;
	//~ // }
	//~ for(double gamma(0.8);gamma<=4;gamma+=0.2){
		//~ for(uint P(1);P<=40;P*=2){
			//~ for(uint H(2);H<3;H*=2){
				//~ vector<uint64_t> kmers=originalkmers;
				//~ uint elts(0);
				//~ for(uint phase(0); phase<P;++phase){
					//~ uint Ncase(kmers.size()*gamma);
					//~ vector<vector<uint64_t>> map(Ncase);
					//~ elts+=Ncase;
					//~ for(uint n(0);n<H-1;++n){
						//~ for(uint i(0);i<kmers.size();++i){
							//~ uint64_t h(iterHash64(kmers[i],n,hashmode)%Ncase);
							//~ // cout<<h<<endl;
							//~ map[h].push_back(kmers[i]);
						//~ }
						//~ kmers={};
						//~ for(uint i(0);i<Ncase;++i){
							//~ if(map[i].size()>1){
								//~ for(uint j(0);j<map[i].size();++j){
									//~ kmers.push_back(map[i][j]);
								//~ }
								//~ map[i]={};
							//~ }
						//~ }
					//~ }
				//~ }
				//~ if((double)elts*nbits(H-1)/Nelement<minbitsbyelts and kmers.size()<unplaced_tolerance*Nelement){
					//~ cout
						//~ <<"Gamma: "<<gamma
						//~ <<"\tH: "<<H-1
						//~ <<"\tP: "<<P
						//~ <<"\tBits/elts: "<<(double)elts*nbits(H-1)/Nelement
						//~ <<"\t% Unplaced element:   "<<(double)100*kmers.size()/Nelement<<" %"<<endl;
					//~ }

			//~ }
			//~ // cout<<endl;
		//~ }
		//~ // cout<<endl<<endl;;
	//~ }

	//~ //
	//~ // for(uint i(0);i<Nelement;++i){
	//~ // 	uint64_t h(hash64(kmers[i])%Nelement);
	//~ // 	// cout<<h<<endl;
	//~ // 	map[h].push_back(kmers[i]);
	//~ // }
	//~ // kmers={};
	//~ //
	//~ // for(uint i(0);i<kmers.size();++i){
	//~ // 	if(map[i].size()>1){
	//~ // 		for(uint j(0);j<map[i].size();++j){
	//~ // 			kmers.push_back(map[i][j]);
	//~ // 		}
	//~ // 		map[i]={};
	//~ // 	}
	//~ // }
	//~ // cout<<"Phase "<<n<<endl;
	//~ // cout<<kmers.size()<<endl;
	//~ // cout<<100*kmers.size()/Nelement<<" %"<<endl;




	return 0;
}
