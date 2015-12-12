#include "utils.h"
#include "xor.h"
#include <fstream>



using namespace std;



int main(int argc, char ** argv){
	srand (time(NULL));
	uint Nelement(1000*1000);
	double gamma(1.1);
	uint Ncase(gamma*Nelement);
	uint k(31);
	cout<<"gamma :"<<gamma<<endl;
	vector<uint64_t> originalkmers(Nelement);
	for(uint i(0);i<Nelement;++i){
		originalkmers[i]=(getRepresent(seq2intStranded(randomSeq(k)),k));
		// cout<<randomSeq(k)<<endl;
		// cout<<kmer[i]<<endl;
	}

	for(uint H(2);H<257;H*=2){
		vector<uint64_t> kmers=originalkmers;
		vector<vector<uint64_t>> map(Ncase);

		for(uint n(0);n<H-1;++n){
			for(uint i(0);i<kmers.size();++i){
				uint64_t h(iterHash64(kmers[i],n)%Ncase);
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
			// cout<<"Phase "<<n<<endl;
			// cout<<kmers.size()<<endl;
			// cout<<100*kmers.size()/Nelement<<" %"<<endl;
		}
		cout<<"Number of hash function:    "<<H-1<<endl;
		cout<<"Percent Unplaced element:  "<<(double)100*kmers.size()/Nelement<<" %"<<endl<<endl;
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
