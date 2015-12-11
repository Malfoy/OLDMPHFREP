#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
// #include <iostream>


uint32_t xor32(uint32_t y){
	y^=(y<<13); y^=(y>>17); return (y^=(y<<15));
}

uint64_t xor64( uint64_t x){
	static uint64_t y = 362436069;
	uint64_t t = (x^(x<<10)); x = y; return y = (y^(y>>10))^(t^(t>>13));
}

uint32_t xor96(uint32_t x){
	static uint32_t  y = 362436069, z = 521288629;
	uint32_t t = (x^(x<<10)); x = y; y = z; return z = (z^(z>>26))^(t^(t>>5));
}

uint32_t xor128(uint32_t x){
	static uint32_t  y = 362436069, z = 521288629,
	                w = 88675123;
	uint32_t t = (x^(x<<11)); x = y; y = z; z = w; return w = (w^(w>>19))^(t^(t>>8));
}

uint32_t xor160(uint32_t x){
	static uint32_t  y = 362436069, z = 521288629,
	                w = 88675123, v = 5783321;
	uint32_t t = (x^(x<<2)); x = y; y = z; z = w; w = v; return v = (v^(v>>4))^(t^(t>>1));
}

uint32_t xorwow(uint32_t x){
	static uint32_t  y = 362436069, z = 521288629,
	                w = 88675123, v = 5783321, d = 6615241;
	uint32_t t = (x^(x>>2)); x = y; y = z; z = w; w = v; v = (v^(v<<4))^(t^(t<<1)); return (d+=362437)+v;
}

#define ROT32(x, y) ((x << y) | (x >> (32 - y))) // avoid effor
uint32_t murmur3_32(uint32_t key, uint32_t seed) {
	static const uint32_t c1 = 0xcc9e2d51;
	static const uint32_t c2 = 0x1b873593;
	static const uint32_t r1 = 15;
	static const uint32_t r2 = 13;
	static const uint32_t m = 5;
	static const uint32_t n = 0xe6546b64;

	uint32_t hash = seed;

	uint32_t k;
		k = key;
		k *= c1;
		k = ROT32(k, r1);
		k *= c2;

		hash ^= k;
		hash = ROT32(hash, r2) * m + n;



	hash ^= (hash >> 16);
	hash *= 0x85ebca6b;
	hash ^= (hash >> 13);
	hash *= 0xc2b2ae35;
	hash ^= (hash >> 16);

	return hash;
}


// uint64_t hash64(uint64_t v){
//     v *= 2654435761;
//     return v >> 32;
// }

//simplest (fastest ?) hash function
uint64_t xorshift64(uint64_t x) {
	// cout<<x<<endl;
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	// cout<<x * UINT64_C(2685821657736338717)<<endl;cin.get();
	return x * UINT64_C(2685821657736338717);
}



uint64_t korenXor(uint64_t x){
	x ^= (x << 21);
	x ^= (x >> 35);
	x ^= (x << 4);
	return x;
}


uint64_t hash64( uint64_t u ){
	return xorshift64(u);
	return korenXor(u);
	// return murmur3_32(u, 69);
	uint64_t v = u * 3935559000370003845 + 2691343689449507681;

	v ^= v >> 21;
	v ^= v << 37;
	v ^= v >>  4;
	v *= 4768777513237032717;

	v ^= v << 20;
	v ^= v >> 41;
	v ^= v <<  5;

	return v;
}


uint64_t iterHash64( uint64_t u , int n){
	if(n==0){
		return hash64(u);
	}
	else{
		return hash64(iterHash64(u, n-1));
	}
}
