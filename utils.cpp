#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include "utils.h"


using namespace std;


char randNuc(){
	switch (rand()%4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}


string randomSeq(uint length){
	auto randchar=[]() -> char{
		const char charset[] ="ATCG";
		const uint max_index = (sizeof(charset) - 1);
		return charset[ rand() % max_index ];
	};
	string str(length,0);
	generate_n( str.begin(), length, randchar );
	return str;
}


//TODO:addinsertion deletion
string mutate(string read,int n){
	for(int i(0); i<n; ++i){
		int position(rand()%read.size());
		read[position]=randNuc();
	}
	return read;
}


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


char revcomp (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	return 'A';//error
}


string reversecomplement (const string& s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--){
		rc += revcomp(s[i]);
	}
	return rc;
}


vector<string> getReads(string& ReadFile, uint n){
	vector<string>reads;
	string read,header,inter;
	char c;
	ifstream readFile(ReadFile);
	for(uint i(0);i<n;++i){
		getline(readFile,header);
		getline(readFile,read);
	point:
		c=readFile.peek();
		if(c=='>'){
			if(read.size()>2){
				bool fail(false);
				for(uint j(0);(j)<read.size();++j){
					if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G'){
						fail=true;
						break;
					}
				}
				if(!fail){
					reads.push_back(read);
				}
			}
			read="";
		}else{
			if(!readFile.eof()){
				getline(readFile,inter);
				read+=inter;
				goto point;
			}else{
				if(read.size()>2){
					bool fail(false);
					for(uint j(0);(j)<read.size();++j){
						if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G'){
							fail=true;
							break;
						}
					}
					if(!fail){
						reads.push_back(read);
					}
				}
				return reads;
			}
		}
	}
	return reads;
}


//A sequence and its reverse complement are represented by the same sequence
string getRepresent (const string& str){
	return (min(str,reversecomplement(str)));
}


string getRepresent2(const string& s){
	for (int i = 0; i < (int)s.length(); i++) {
		char c = revcomp(s[s.length() - 1 - i]);
		if (s[i] < c) {
			return s;
		} else if (s[i] > c) {
			return reversecomplement(s);
		}
	}
	return s;
}


minimizer seq2int(const string& seq){
	string str(getRepresent(seq));
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(str[i]);
	}
	return res;
}

char int2nuc(uint n){
switch(n){
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'G';
	case 3:
		return 'T';
}
return 'X';
}


void int2seq(minimizer min, uint n){
	string res;
	uint nuc;
	for(uint i(0);i<n;++i){
		nuc=min%4;
		res+=int2nuc(nuc);
		min>>=2;
	}
	reverse(res.begin(),res.end());
	cout<<res<<endl;
}


minimizer seq2intStranded(const string& seq){
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(seq[i]);
	}
	return res;
}


bool equalStr(const string& seq1, const string& seq2){
	uint size(min(seq1.size(),seq2.size()));
	return (seq1.substr(0,size))==seq2.substr(0,size);
}


//rewrite with ternary operator
bool isCorrect(const string& seq,const string& ref){
	for(uint i(0); i<seq.size(); ++i){
		if(seq[i]!=ref[i]){
			if(seq[i+1]==ref[i]){
				return equalStr(seq.substr(i+2),ref.substr(i+1));
			}
			if(seq[i]==ref[i+1]){
				return equalStr(seq.substr(i+1),ref.substr(i+2));
			}
			return (seq.substr(i+1)==ref.substr(i+1));
		}
	}
	return true;
}


double percentStrandedErrors(uint k, const string& seq, const unordered_multimap<string, string>& genomicKmers, char nuc){
	double inter(0);
	string kmer;
	kmer.reserve(k);
	uint i(0);
	for(; i+k<=seq.size(); ++i){
		kmer=seq.substr(i,k);
		if(kmer.size()!=k){cout<<"wtf"<<endl;}
		auto range(genomicKmers.equal_range(kmer.substr(0,nuc)));
		for(auto it(range.first); it!=range.second; ++it){
			if(isCorrect(kmer.substr(nuc),it->second)){
				++inter;
				break;
			}else{}
		}
	}
	return double(100*inter/(seq.size()-k+1));;
}


uint sketchOrderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2){
	uint res(0);
	unordered_set<minimizer> used;
	for(uint i(0); i<min(sketch1.size(),sketch2.size()); ++i){
		if(used.unordered_set::count(sketch1[i])==0){
			if(sketch1[i]==sketch2[i]){
				++res;
				used.insert(sketch2[i]);
			}
		}
	}
	return res;
}


uint sketchUnorderedComparison(const vector<minimizer>& sketch1, const vector<minimizer>& sketch2){
	uint res(0);
	unordered_set<minimizer> used;
	unordered_set<minimizer> minimizerSet;
	for(uint i(0); i<sketch1.size(); ++i){
		minimizerSet.insert(sketch1[i]);
	}
	for(uint i(0); i<sketch2.size(); ++i){
		if(used.unordered_set::count(sketch2[i])==0){
			if(minimizerSet.unordered_set::count(sketch2[i])!=0){
				++res;
				used.insert(sketch2[i]);
			}
		}
	}
	return res;
}


minimizer cat(minimizer seed, minimizer body, uint n){
	minimizer res(seed<<(2*n));
	return res+=body;
}


minimizer rc(minimizer min,uint n){
	minimizer res(0);
	minimizer offset(1);
	offset<<=(2*n-2);
	for(uint i(0); i<n;++i){
		res+=(3-(min%4))*offset;
		min>>=2;
		offset>>=2;
	}
	return res;
}


minimizer getRepresent(minimizer kmer, uint n){
	return (min(kmer,rc(kmer,n)));
}


minimizer getEnd(minimizer kmer, uint n){
	return kmer%((minimizer)1<<(2*n));
}


minimizer getBegin(minimizer kmer, uint n){
	return (kmer>>(2*n));
}


uint sketchUnorderedComparisonError(const unordered_multimap<string, string>& map1, const unordered_multimap<string, string>& map2){
	uint res(0);
	string beg,end;
	for (auto it=map1.begin(); it!=map1.end(); ++it){
		beg=it->first;
		end=it->second;
		auto ret = map2.equal_range(beg);
		for (auto it2=ret.first; it2!=ret.second; ++it2){
			if(isCorrect(end,it2->second)){
				++res;
			}
		}
	}
	return res;
}


double scoreFromAlignment(const string& seq1,const string& seq2){
	size_t match(0);
	for(size_t i(0);i<seq1.size();++i){
		if(seq1[i]==seq2[i]){
			++match;
		}
	}
	double res((100*match)/(seq1.size()));
	return res;
}
