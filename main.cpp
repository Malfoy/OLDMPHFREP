#include "utils.h"
#include "xor.h"
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <algorithm>    // std::count
#include <math.h>




using namespace std;



void nadine(const double gamma, const uint Count, uint H,uint expectedSizeBucket,uint Nessai){
	
	uint Nelement(1*100*1000);
	uint minSizeBucket(Nelement);
	uint maxSizeBucket(0);

	vector<uint64_t> originalkmers(Nelement);
	uint k(31);
	string seq;
	//~ uint P(1);
	for(uint i(0);i<Nelement;++i){
		seq = randomSeq(k);
		originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
	}
	unordered_set<uint64_t> set;
	for(uint i(0);i<Nelement;++i){
		set.insert(originalkmers[i]);
	}
	
	vector<vector<uint64_t>> Buckets(Nelement/expectedSizeBucket);
	
	vector<vector<bool>> finalStruct(Buckets.size());
	vector<uint> HashUsed(Buckets.size());
	int hashmode = 0;
	uint nBitUsed(0);
	for(uint i(0);i<Nelement;++i){
		uint64_t h(iterHash64(originalkmers[i],0,hashmode)%Buckets.size());
		Buckets[h].push_back(originalkmers[i]);
	}
	cout<<"element distribued"<<endl;
	
	for(uint subSetNumber(0);subSetNumber<Buckets.size();++subSetNumber){
		vector<uint64_t> kmers=Buckets[subSetNumber];
		if(kmers.size()>maxSizeBucket){
			maxSizeBucket=kmers.size();
		}
		if(kmers.size()<minSizeBucket){
			minSizeBucket=kmers.size();
		}
		//~ cout<<kmers.size()<<endl;
		uint minOne(kmers.size()*gamma);
		
		for(uint essai(1);essai<=Nessai;++essai){
			vector<bool> approxSet(10*gamma*(kmers.size()/10),false);
			for(uint i(0);i<kmers.size();++i){
				for(uint hashNumber(0);hashNumber<H;++hashNumber){
					uint64_t h(iterHash64(kmers[i],hashNumber+essai*H,hashmode)%approxSet.size());
					approxSet[h]=true;
				}
			}
			uint one(0),zero;
			for(uint i(0);i<approxSet.size();++i){
				if(approxSet[i]){
					one++;
				}else{
					zero++;
				}
			}
			//~ cout<<one<<" "<<zero<<endl;cin.get();
			if(one<minOne){
				//~ cout<<"happen"<<endl;
				finalStruct[subSetNumber]=approxSet;
				HashUsed[subSetNumber]=essai*H;
				minOne=one;
			}
		}
	}
	
	cout<<"set constructed"<<endl;
	
	//~ vector<uint> count(gamma*Nelement,0);
	//~ for(uint i(0);i<Nelement;++i){
		//~ for(uint hashNumber(0);hashNumber<H;++hashNumber){
			//~ uint64_t h(iterHash64(kmers[i],hashNumber,hashmode)%count.size());
				//~ ++count[h];
		//~ }
	//~ }
	
	//~ for(uint i(0);i<count.size();++i){
		//~ if(count[i]>=Count){
			//~ approxSet[i]=true;
		//~ }
	//~ }
	
	
	uint FP(0);
	for(uint test(0);test<Nelement;){
		seq = randomSeq(k);
		uint64_t kmer(getRepresent(seq2intStranded(seq),k));
		if(set.count(kmer)==0){
			++test;
			uint vote(0);
			uint setChosen(iterHash64(kmer,0,hashmode)%Buckets.size());
			if(not finalStruct[setChosen].empty()){
				for(uint hashNumber(0);hashNumber<H;++hashNumber){
					
					uint64_t h(iterHash64(kmer,hashNumber+HashUsed[setChosen],hashmode)%finalStruct[setChosen].size());
					//~ cout<<h<<" "<<finalStruct[setChosen].size()<<endl;
					if(finalStruct[setChosen][h]){
						++vote;
					}
				}
			}
			if(vote==H){
				FP++;
			}
		}
	}
	cout<<"fp estimated"<<endl;
	
	
	//ESTIMATING FALSE NEGATIVE
	uint TP(0),FN(0);
	//~ for(uint i(0);i<Nelement;++i){
		//~ uint64_t kmer(originalkmers[i]);
		//~ uint vote(0);
		//~ for(uint setNumber(0);setNumber<H;++setNumber){
			//~ uint64_t h(iterHash64(kmer,setNumber,hashmode)%approxSet.size());
			//~ if(approxSet[h]){
				//~ ++vote;
			//~ }
		//~ }
		//~ if(vote<P){
			//~ ++FN;
		//~ }
	//~ }
	
	//~ uint one(0),total(0);
	//~ for(uint i(0);i<finalSet.size();++i){
		//~ if(finalSet[i]){
			//~ one++;
		//~ }
		//~ total++;
	//~ }
	uint one(0),zero(0);
	for(uint i(0);i<finalStruct.size();++i){
		for(uint j(0);j<finalStruct[i].size();++j){
			if(finalStruct[i][j]){
				one++;
			}else{
				zero++;
			}
		}
	}
	
	
	double bits((double)(one+zero)/Nelement+(double)Buckets.size()*(log2(Nessai)+3)/Nelement);
	double fprate((double)FP/Nelement);
	cout<<"Number try: "<<Nessai<<endl;
	cout<<"expected bucket size: "<<expectedSizeBucket<<endl;
	cout<<"min/max bucketsize: "<<minSizeBucket<<" "<<maxSizeBucket<<endl;

	//~ cout<<(double)one/total<<endl;;
	cout<<"H: "<<H<<endl;
	//~ cout<<"P: "<<P<<endl;
	//~ cout<<"Count: "<<Count<<endl;
	cout<<"bit to 1 in the sets: "<<(double)one/(one+zero)<<endl;
	cout<<"bit/element: "<<bits<<endl;
	cout<<"FP rate: "<<fprate<<endl;
	//~ cout<<FP<<" "<<Nelement<<endl;
	cout<<"FN rate: "<<(double)FN/Nelement<<endl;
	cout<<"coef: "<<exp(log(fprate)/bits)<<endl;
	cout<<endl;
}


	
	/*
	//RANDOM SET CREATION
	//~ cout<<gamma<<" "<<Count<<endl;
	for(H=2;H<100;H++){
		for(P=H/2;P<=H;P+=1){
		uint Nelement(1*100*1000);
		vector<uint64_t> originalkmers(Nelement);
		uint k(31);
		string seq;
		//~ uint P(1);
		for(uint i(0);i<Nelement;++i){
			seq = randomSeq(k);
			originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
		}
		unordered_set<uint64_t> set;
		for(uint i(0);i<Nelement;++i){
			set.insert(originalkmers[i]);
		}
		int hashmode = 0;
		uint nBitUsed(0);
		vector<uint64_t> kmers=originalkmers;
		
		vector<vector<bool>> approxSet;
		
		//~ cout<<1<<endl;
		//BUILDING THE GBF
		for(uint setNumber(0);setNumber<H;++setNumber){
			uint casesNumber(kmers.size()*gamma);
			vector<uint8_t> counting(casesNumber,0);
			for(uint i(0);i<kmers.size();++i){
				uint64_t h(iterHash64(kmers[i],setNumber,hashmode)%casesNumber);
				++counting[h];
			}
			vector<uint64_t> unplaced;
			vector<bool> bitSet(casesNumber,false);
			for(uint i(0);i<kmers.size();++i){
				uint64_t h(iterHash64(kmers[i],setNumber,hashmode)%casesNumber);
				if(counting[h]<Count){
					unplaced.push_back(kmers[i]);
				}else{
					bitSet[h]=true;
				}
			}
			//~ kmers=unplaced;
			approxSet.push_back(bitSet);
		}
		
		//~ cout<<2<<endl;
		//ESTIMATING FALSE POSITIVE
		uint FP(0);
		for(uint test(0);test<Nelement;){
			seq = randomSeq(k);
			uint64_t kmer(getRepresent(seq2intStranded(seq),k));
			if(set.count(kmer)==0){
				++test;
				uint vote(0);
				for(uint setNumber(0);setNumber<H;++setNumber){
					if(approxSet[setNumber].size()!=0){
						uint64_t h(iterHash64(kmer,setNumber,hashmode)%approxSet[setNumber].size());
						if(approxSet[setNumber][h]){
							++vote;
						}
					}
				}
				if(vote>=P){
					FP++;
				}
			}
		}
		
		//~ cout<<3<<endl;
		//ESTIMATING FALSE NEGATIVE
		uint TP(0),FN(0);
		for(uint i(0);i<Nelement;++i){
			uint64_t kmer(originalkmers[i]);
			uint vote(0);
			for(uint setNumber(0);setNumber<H;++setNumber){
				if(approxSet[setNumber].size()!=0){
					uint64_t h(iterHash64(kmer,setNumber,hashmode)%approxSet[setNumber].size());
					if(approxSet[setNumber][h]){
						++vote;
					}
				}
			}
			if(vote<P){
				++FN;
			}
		}
		
		//~ cout<<4<<endl;
		//NUMBER OF ONE CHEKING FOR ENTROPY
		uint one(0),totalbit(0);
		for(uint setNumber(0);setNumber<H;++setNumber){
			//~ for(uint test(0);test<Nelement;){
				for(uint i(0);i<approxSet[setNumber].size();++i){
					if(approxSet[setNumber][i]){
						one++;
					}
					++totalbit;	
				}
			//~ }
		}
		double FPrate((double)FP/Nelement),FNrate((double)FN/Nelement), FPMax(0.01),FNMax(0.01);
		if(FPrate<FPMax and ((double)gamma*Nelement*H/Nelement+FNrate*10) <12){
			cout<<"P: "<<P<<endl;
			cout<<"H: "<<H<<endl;
			cout<<"Count: "<<Count<<endl;
			cout<<"bit to 1 in the sets: "<<(double)one/totalbit<<endl;
			//~ cout<<"unplaced element rate: "<<(double)kmers.size()/Nelement<<endl;
			//~ cout<<"bit used: "<<gamma*Nelement*H <<endl;
			cout<<"bit/element: "<<((double)gamma*Nelement*H/Nelement)+FNrate*10<<endl;
			cout<<"FP rate: "<<((double)FP/Nelement)<<endl;
			//~ cout<<FP<<" "<<Nelement<<endl;
			cout<<"FN rate: "<<(double)FN/Nelement<<endl;
			cout<<endl;
			//~ return;
		}
	}
	}
}
*/
	

// cmdline: ./bmean [<genome.fasta	xxhash]
int main(int argc, char ** argv){
	srand (time(NULL));
	double gamma(1.2);
	uint Count(1);
	uint lolmin(80);
	uint lolmax(80);
	for(uint lol(lolmin);lol<=lolmax;lol+=1){
		nadine(gamma,Count,1,lol,255);
	}
	
	//~ cout<<"count: "<<Count<<endl;
	//~ nadine(gamma,Count);
	
	//~ gamma=0.59;
	//~ Count=2;
	//~ nadine(gamma,Count);
	
		//~ gamma=0.37;
		//~ Count=3;
		//~ nadine(gamma,Count);
		
		//~ gamma=0.21;
		//~ Count=5;
		//~ nadine(gamma,Count);
	
	//~ gamma=0.103;
	//~ Count=10;
	//~ nadine(gamma,Count);
	
	//~ gamma=0.051;
	//~ Count=20;
	//~ nadine(gamma,Count,H);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*
	srand (time(NULL));
	uint Nelement(1*1000*1000);
	vector<uint64_t> originalkmers(Nelement);
	uint k(31);
	string seq;
	for(uint i(0);i<Nelement;++i){
		seq = randomSeq(k);
		originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
	}
	unordered_set<uint64_t> set;
	for(uint i(0);i<Nelement;++i){
		set.insert(originalkmers[i]);
	}
	
	uint nbSet(256);
	uint finger(8);
	uint hashSize(16);
	uint hashMax(256*256);	
	for(uint H(1);H<=1;H*=2){
		cout<<"H: "<<H<<endl;
		double gamma(1);
		int hashmode(0);
		vector<vector<uint8_t>> setSet(256);
		vector<vector<bool>> setBool(256);
		vector<vector<uint64_t>> buckets(256);
		vector<uint64_t> unplaced;
		cout<<1<<endl;
		//foreach element
		for(uint i(0);i<Nelement;++i){
			uint64_t kmer(originalkmers[i]);
			//~ bool placed=false;
			for(uint ii(0);ii<H;++ii){
				uint64_t h(iterHash64(kmer,0,hashmode)%hashMax);
				buckets[h/256].push_back(kmer);
			}
		}
		for(uint i(0);i<256;++i){
			cout<<buckets[i].size()<<endl;
		}
		cout<<2<<endl;
		
		for(uint i(0); i<256;++i){
			setBool[i]=vector<bool>(buckets[i].size()*gamma,false);
			setSet[i]=vector<uint8_t>(buckets[i].size()*gamma,0);
			for(uint j(0);j<buckets[i].size();++j){
				uint64_t kmer(buckets[i][j]);
				uint64_t h(iterHash64(kmer,0,hashmode)%hashMax);
				uint hpos(iterHash64(kmer,1,hashmode)%setSet[i].size());
				if(setBool[i][hpos]==false){
					setSet[i][hpos]=h%256;
					setBool[i][hpos]=true;
				}else{
					unplaced.push_back(kmer);
				}
			}
		}
		
		cout<<100*(double)unplaced.size()/Nelement<<"%"<<endl;
		
		uint test(0);
		uint FP(0);
		while(test<Nelement){
			seq = randomSeq(k);
			uint64_t kmer(getRepresent(seq2intStranded(seq),k));
			if(set.count(kmer)==0){
				++test;
				uint64_t h(iterHash64(kmer,0,hashmode)%hashMax);
				uint hpos(iterHash64(kmer,1,hashmode)%setSet[h/256].size());
				if(setBool[h/256][hpos]==true and setSet[h/256][hpos]==h%256){
					++FP;
					goto lol2;
				}
			}
			lol2:;
		}
		cout<<FP<<endl;
		//~ cout<<test<<endl;
		cout<<(double)(1000*FP)/(double)test<<endl;
	}
	/*
	srand (time(NULL));
	uint Nelement(1*100*1000);
	vector<uint64_t> originalkmers(Nelement);
	uint k(31);
	string seq;
	for(uint i(0);i<Nelement;++i){
		seq = randomSeq(k);
		originalkmers[i]=(getRepresent(seq2intStranded(seq),k));
	}
	unordered_set<uint64_t> set;
	for(uint i(0);i<Nelement;++i){
		set.insert(originalkmers[i]);
	}
	double gamma(1);
	int hashmode(0);

	uint nbBloom(1),H(1),casesNumber(gamma*Nelement);
	vector<vector<bool>> blooms;
	for(uint i(0);i<nbBloom;++i){
		vector<bool> lol(casesNumber,false);
		blooms.push_back(lol);
	}
	uint offset(0);
	
	//foreach element
	for(uint i(0);i<Nelement;++i){
		uint64_t kmer(originalkmers[i]);
		uint minchangeRequired(H+1);
		uint selectedBloom(nbBloom+1);
		//foreach bloom filter
		for(uint bloomNumber(0);bloomNumber<nbBloom;++bloomNumber){
			//~ cout<<"bloom: "<<bloomNumber<<endl;
			uint changeRequired(0);
			//forech hash function
			for(uint hashnumber(0);hashnumber<H;++hashnumber){
				uint64_t h(iterHash64(kmer,hashnumber+H*(bloomNumber+offset)%nbBloom,hashmode)%casesNumber);=++
				//~ cout<<hashnumber+H*bloomNumber<<endl;
				//~ cout<<h<<endl;
				if(not blooms[(bloomNumber+offset)%nbBloom][h]){
					++changeRequired;
				}
			}
			if(changeRequired<minchangeRequired){
				//~ if(changeRequired<H){
					//~ cout<<"lol"<<endl;
				//~ }
				minchangeRequired=changeRequired;
				selectedBloom=(bloomNumber+offset)%nbBloom;
			}
		}
		//~ cout<<
		for(uint hashnumber(0);hashnumber<H;++hashnumber){
			uint64_t h(iterHash64(kmer,hashnumber+H*selectedBloom,hashmode)%casesNumber);
			blooms[selectedBloom][h]=true;
		}
		offset++;
	}
	
	uint FP(0);
	for(uint test(0);test<Nelement;){
		seq = randomSeq(k);
		uint64_t kmer(getRepresent(seq2intStranded(seq),k));
		if(set.count(kmer)==0){
			++test;
			for(uint ii(0);ii<nbBloom;++ii){
				uint hit(0);
				for(uint iii(0);iii<H;++iii){
					uint64_t h(iterHash64(kmer,iii+H*ii,hashmode)%casesNumber);
					if(blooms[ii][h]){
						++hit;
					}
				}
				if(hit==H){
					++FP;
					goto lol2;
				}
			}
		}
		lol2:;
	}
	for(uint i(0);i<nbBloom;++i){
		cout<<count (blooms[i].begin(),blooms[i].end(), true)<<endl;
	}
	
	cout<<"H: "<<H<<endl;
	cout<<"nb Bloom: "<<nbBloom<<endl;
	cout<<"bit used: "<<casesNumber*nbBloom<<endl;
	cout<<"bit/element: "<<((double)casesNumber*nbBloom/Nelement)<<endl;
	cout<<"FP rate: "<<((double)FP/Nelement)<<endl;
	//~ cout<<FP<<" "<<Nelement<<endl;
	//~ cout<<"FN: "<<FN-kmers.size()<<endl;
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
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
	double gamma(2);
	uint count(2);
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



*/
	return 0;
}
