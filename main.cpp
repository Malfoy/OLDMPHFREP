#include "utils.h"
#include "xor.h"
#include <fstream>
#include <cstring>



using namespace std;


// cmdline: ./bmean [<genome.fasta>] [--xxhash]
int main(int argc, char ** argv){
	srand (time(NULL));
	uint Nelement(1*1000*1000);
    double unplaced_tolerance(0.02); // (corresponds to 2%); used to be 0.1%
	// double gamma(2);
	// uint Ncase(gamma*Nelement);
	uint k(31);
	double minbitsbyelts(5);
	// cout<<"gamma :"<<gamma<<endl;
	vector<uint64_t> originalkmers(Nelement);
    string genome = "";
    if (argc>1 && argv[1][0] != '-') { // whether a genome is passed as the 1st argument
        genome = loadGenome(string(argv[1]));
        if (genome.size() < Nelement-k+1) { cout << "not enough genomic kmers in " << string(argv[1]) << endl; exit(1);}
    }
    int hashmode = 0;
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i],"--xxhash")==0) { cout << "using xxhash" << endl; hashmode = 1; }
	for(uint i(0);i<Nelement;++i){
        string seq = genome.size() > 0 ? genomeKmer(genome,i,k) : randomSeq(k);
		originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
	}

	// for(uint H(2);H<257;H*=2){
	// 	vector<uint64_t> kmers=originalkmers;
	// 	vector<vector<uint64_t>> map(Ncase);
	// 	for(uint n(0);n<H-1;++n){
	// 		for(uint i(0);i<kmers.size();++i){
	// 			uint64_t h(iterHash64(kmers[i],n)%Ncase);
	// 			map[h].push_back(kmers[i]);
	// 		}
	// 		kmers={};
	// 		for(uint i(0);i<Ncase;++i){
	// 			if(map[i].size()>1){
	// 				for(uint j(0);j<map[i].size();++j){
	// 					kmers.push_back(map[i][j]);
	// 				}
	// 				map[i]={};
	// 			}
	// 		}
	// 	}
	// 	cout<<"Number of hash function:    "<<H-1<<endl;
	// 	cout<<"Percent Unplaced element:  "<<(double)100*kmers.size()/Nelement<<" %"<<endl<<endl;
	// }
	for(double gamma(0.8);gamma<=4;gamma+=0.2){
		for(uint P(1);P<=40;P*=2){
			for(uint H(2);H<3;H*=2){
				vector<uint64_t> kmers=originalkmers;
				uint elts(0);
				for(uint phase(0); phase<P;++phase){
					uint Ncase(kmers.size()*gamma);
					vector<vector<uint64_t>> map(Ncase);
					elts+=Ncase;
					for(uint n(0);n<H-1;++n){
						for(uint i(0);i<kmers.size();++i){
							uint64_t h(iterHash64(kmers[i],n,hashmode)%Ncase);
							// cout<<h<<endl;
							map[h].push_back(kmers[i]);
						}
						kmers={};
						for(uint i(0);i<Ncase;++i){
							if(map[i].size()>1){
								for(uint j(0);j<map[i].size();++j){
									kmers.push_back(map[i][j]);
								}
								map[i]={};
							}
						}
					}
				}
				if((double)elts*nbits(H-1)/Nelement<minbitsbyelts and kmers.size()<unplaced_tolerance*Nelement){
					cout
						<<"Gamma: "<<gamma
						<<"\tH: "<<H-1
						<<"\tP: "<<P
						<<"\tBits/elts: "<<(double)elts*nbits(H-1)/Nelement
						<<"\t% Unplaced element:   "<<(double)100*kmers.size()/Nelement<<" %"<<endl;
					}

			}
			// cout<<endl;
		}
		// cout<<endl<<endl;;
	}

	//
	// for(uint i(0);i<Nelement;++i){
	// 	uint64_t h(hash64(kmers[i])%Nelement);
	// 	// cout<<h<<endl;
	// 	map[h].push_back(kmers[i]);
	// }
	// kmers={};
	//
	// for(uint i(0);i<kmers.size();++i){
	// 	if(map[i].size()>1){
	// 		for(uint j(0);j<map[i].size();++j){
	// 			kmers.push_back(map[i][j]);
	// 		}
	// 		map[i]={};
	// 	}
	// }
	// cout<<"Phase "<<n<<endl;
	// cout<<kmers.size()<<endl;
	// cout<<100*kmers.size()/Nelement<<" %"<<endl;




	return 0;
}
