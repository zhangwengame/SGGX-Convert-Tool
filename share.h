#if !defined(__SHARE_H_)
#define __SHARE_H_
#include<math.h>
typedef unsigned long long uint64_t;
typedef unsigned int       uint32_t;
struct Point2{
	Point2() :x(0.0), y(0.0){}
	Point2(float x, float y) :x(x), y(y){}
	float x, y;
};
inline float sobol2Single(uint32_t n, uint32_t scramble = 0U) {
	for (uint32_t v = 1U << 31; n != 0; n >>= 1, v ^= v >> 1)
		if (n & 1)
			scramble ^= v;
	return (float)scramble / (float)(1ULL << 32);
}
inline float radicalInverse2Single(uint32_t n, uint32_t scramble = 0U) {
	/* Efficiently reverse the bits in 'n' using binary operations */
#if (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))) || defined(__clang__)
	n = __builtin_bswap32(n);
#else
	n = (n << 16) | (n >> 16);
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
#endif
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);

	// Account for the available precision and scramble
	n = (n >> (32 - 24)) ^ (scramble & ~- (1 << 24));

	return (float)n / (float)(1U << 24);
}
inline Point2 sample02Single(uint32_t n, uint32_t scramble[2]) {
	return Point2(
		radicalInverse2Single(n, scramble[0]),
		sobol2Single(n, scramble[1])
		);
}

struct vec3{
	vec3() :x(0.0), y(0.0), z(0.0){}
	vec3(float x, float y, float z) :x(x), y(y), z(z){}
	float x, y, z;
};
inline float dot(const vec3 &a, const  vec3 &b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline vec3 operator +(const vec3 &a, const  vec3 &b){
	return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline vec3 operator *(const float &a, const  vec3 &b){
	return vec3(a*b.x, a*b.y, a*b.z);
}
inline vec3 normalize(const vec3 &a){
	float length = sqrtf(dot(a, a));
	return vec3(a.x / length, a.y / length, a.z / length);
}

#endif