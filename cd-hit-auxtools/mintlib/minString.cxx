//=================================================================
// This file is a part of the Minimum Template Library.
// By Limin Fu (phoolimin@gmail.com, lmfu@ucsd.edu)
//================================================================= 

#include<ctype.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include"minString.hxx"

BEGIN_NS_MIN

String::String( char ch )
{
	data = NULL;
	size = bufsize = 0;
	SetBytes( & ch, 1 );
}
String::String( const char *s )
{
	data = NULL;
	size = bufsize = 0;
	*this = s;
}
String::String( const char *s, int n )
{
	data = NULL;
	size = bufsize = 0;
	SetBytes( s, n );
}
String::String( const String & other )
{
	data = NULL;
	size = bufsize = 0;
	SetBytes( other.data, other.size );
}

void String::Clear()
{
	if( data ) free( data );
	data = NULL;
	size = bufsize = 0;
}
void String::Reset()
{
	if( size == 0 ) return;
	size = 0;
	data[0] = 0;
}
void String::Reserve( int n )
{
	if( bufsize >= n ) return;
	data = (char*) realloc( data, (n+1)*sizeof(char) );
	bufsize = n;
}
void String::ResetSize( int n )
{
	if( size < n ) Resize( n, 0 );
	size = n;
	data[size] = 0;
}
void String::Resize( int n, char ch )
{
	int i;
	if( bufsize != n ){
		data = (char*) realloc( data, (n+1)*sizeof(char) );
		bufsize = n;
	}
	for(i=size; i<n; i++) data[i] = ch;
	if( data ) data[n] = 0;
	size = bufsize = n;
}
void String::SetBytes( const char *s, int n )
{
	Resize( n );
	if( s ) memcpy( data, s, n*sizeof(char) );
	if( data ) data[n] = 0;
}
void String::AppendBytes( const char *s, int n )
{
	Reserve( size + n );
	memcpy( data + size, s, n*sizeof(char) );
	size += n;
	if( data ) data[size] = 0;
}
void String::Insert( char ch, int at )
{
	if( at < 0 || at > size ) return;
	Reserve( size + 1 );
	if( at < size ) memmove( data + at + 1, data + at, (size - at) * sizeof(char) );
	data[at] = ch;
	data[++size] = '\0';
}
void String::Insert( char *chs, int at )
{
	if( chs == NULL ) return;
	if( at < 0 || at > size ) return;
	int i, N = strlen( chs );
	Reserve( size + N );
	if( at < size ) memmove( data + at + N, data + at, (size - at) * sizeof(char) );
	for(i=0; i<N; i++) data[at+i] = chs[i];
	size += N;
	data[size] = '\0';
}
void String::Erase( int start, int count )
{
	int rest = 0;
	if( count < 0 ) count = size;
	if( start < 0 or start >= size ) return;
	if( start + count > size ) count = size - start;
	rest = size - (start + count);
	if( rest ) memmove( data + start, data + start + count, rest * sizeof(char) );
	size -= count;
	if( data ) data[size] = 0;
}
void String::Chop( int count )
{
	if( size ){
		if( count > size ) count = size;
		Erase( size-count, count );
	}
}
void String::ToUpper()
{
	for(int i=0; i<size; i++) data[i] = toupper( data[i] );
}
void String::ToLower()
{
	for(int i=0; i<size; i++) data[i] = tolower( data[i] );
}
int String::Find( char ch )const
{
	int i;
	for(i=0; i<size; i++) if( data[i] == ch ) return i+1;
	return 0;
}
int String::Find( const String & s )const
{
	int i, j, m = s.Size();
	for(i=0; (i+m)<=size; i++){
		for(j=0; j<m; j++) if( data[i+j] != s[j] ) break;
		if( j == m ) return i+1;
	}
	return 0;
}
String& String::operator=( const char ch )
{
	SetBytes( & ch, 1 );
	return *this;
}
String& String::operator=( const char *s )
{
	if( s == NULL ){
		Clear();
		return *this;
	}
	SetBytes( s, strlen(s) );
	return *this;
}
String& String::operator=( const String & s )
{
	SetBytes( s.data, s.size );
	return *this;
}
void String::operator+=( const char ch )
{
	Reserve( size + 1 );
	data[size] = ch;
	size += 1;
	data[size] = 0;
}
void String::operator+=( const char *s )
{
	if( s == NULL ) return;
	int n = strlen( s );
	Reserve( size + n + 1 );
	memcpy( data + size, s, n*sizeof(char) );
	size += n;
	if( data ) data[size] = 0;
}
void String::operator+=( const String &s )
{
	Reserve( size + s.size );
	memcpy( data + size, s.data, s.size*sizeof(char) );
	size += s.size;
	if( data ) data[size] = 0;
}

String String::FromHex()const
{
	String s;
	char enc[256];
	int i;
	memset( enc, 0, sizeof(char) );
	for(i=0; i<10; i++) enc[i+'0'] = i;
	for(i=0; i<6; i++) enc[i+'a'] = enc[i+'A'] = i+10;
	for(i=0; i<size; i+=2){
		char c1 = enc[ (int)data[i] ];
		char c2 = enc[ (int)data[i+1] ];
		s += (c1<<4)|c2;
	}
	return s;
}
String String::ToHex()const
{
	const char *hex = "0123456789abcdef";
	String shex;
	int i;
	for(i=0; i<size; i++){
		shex += hex[ data[i] >> 4 ];
		shex += hex[ data[i] & 0xf ];
	}
	return shex;
}


static void MD5_Append( String & md5, uint32_t h )
{
	const char *hex = "0123456789abcdef";
	uint32_t k;
	k = (h >> 0) & 0xff;  md5 += hex[k>>4];  md5 += hex[k&0xf];
	k = (h >> 8) & 0xff;  md5 += hex[k>>4];  md5 += hex[k&0xf];
	k = (h >>16) & 0xff;  md5 += hex[k>>4];  md5 += hex[k&0xf];
	k = (h >>24) & 0xff;  md5 += hex[k>>4];  md5 += hex[k&0xf];
}
static void MD5_Update( uint32_t H[4], uint32_t W[16], uint32_t K[64] )
{
	static uint32_t R[64] = {
		7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,
		5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,
		4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,
		6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21
	};
	uint32_t A = H[0];
	uint32_t B = H[1];
	uint32_t C = H[2];
	uint32_t D = H[3];
	uint32_t k;
	for(k=0; k<16; k++){
		uint32_t f = (B & C) | ((~B) & D);
		uint32_t g = k;
		uint32_t t = D;
		uint32_t x = A + f + K[k] + W[g];
		D = C;
		C = B;
		B = B + ((x << R[k]) | (x >> (32-R[k])));
		A = t;
	}
	for(k=16; k<32; k++){
		uint32_t f = (D & B) | ((~D) & C);
		uint32_t g = (k*5 + 1) % 16;
		uint32_t t = D;
		uint32_t x = A + f + K[k] + W[g];
		D = C;
		C = B;
		B = B + ((x << R[k]) | (x >> (32-R[k])));
		A = t;
	}
	for(k=32; k<48; k++){
		uint32_t f = B ^ C ^ D;
		uint32_t g = (k*3 + 5) % 16;
		uint32_t t = D;
		uint32_t x = A + f + K[k] + W[g];
		D = C;
		C = B;
		B = B + ((x << R[k]) | (x >> (32-R[k])));
		A = t;
	}
	for(k=48; k<64; k++){
		uint32_t f = C ^ (B | (~D));
		uint32_t g = (k*7) % 16;
		uint32_t t = D;
		uint32_t x = A + f + K[k] + W[g];
		D = C;
		C = B;
		B = B + ((x << R[k]) | (x >> (32-R[k])));
		A = t;
	}
	H[0] += A;
	H[1] += B;
	H[2] += C;
	H[3] += D;
}

String String::MD5()const
{
	String padding;
	uint64_t n, twop32 = ((uint64_t)1)<<32;
	uint32_t H[4] = { 0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476 };
	uint32_t K[64], W[16];
	int32_t i, k, m, chunks = size / 64;
	uint8_t *data = (uint8_t*) this->data;
	
	for(i=0; i<64; i++) K[i] = (uint32_t) floor( fabs( sin(i+1) ) * twop32 );
	for(i=0; i<chunks; i++){
		for(k=0; k<16; k++){
			uint32_t b = i*64 + k*4;
			uint32_t m = data[b];
			m |= ((uint32_t)data[b+1])<<8;
			m |= ((uint32_t)data[b+2])<<16;
			m |= ((uint32_t)data[b+3])<<24;
			W[k] = m;
		}
		MD5_Update( H, W, K );
	}
	padding.Reserve( 128 );
	padding.size = 64;
	data = (uint8_t*) padding.data;
	m = size - chunks*64;
	if( m ) memcpy( padding.data, this->data + chunks*64, m*sizeof(char) );
	if( m + 8 > 64 ) padding.size = 128;
	chunks = padding.size / 64;
	
	data[m] = 1<<7; // first bit 1 followed by bit 0s;
	for(i=m+1; i<padding.size-8; i++) data[i] = 0;
	n = size * 8;
	// last 64 bits to store the string size in little endian:
	data[i] = n & 0xff;
	data[i+1] = (n >> 8) & 0xff;
	data[i+2] = (n >> 16) & 0xff;
	data[i+3] = (n >> 24) & 0xff;
	data[i+4] = (n >> 32) & 0xff;
	data[i+5] = (n >> 40) & 0xff;
	data[i+6] = (n >> 48) & 0xff;
	data[i+7] = (n >> 56) & 0xff;
	for(i=0; i<chunks; i++){
		for(k=0; k<16; k++){
			uint32_t b = i*64 + k*4;
			uint32_t m = data[b];
			m |= ((uint32_t)data[b+1])<<8;
			m |= ((uint32_t)data[b+2])<<16;
			m |= ((uint32_t)data[b+3])<<24;
			W[k] = m;
		}
		MD5_Update( H, W, K );
	}
	padding.size = 0;
	MD5_Append( padding, H[0] );
	MD5_Append( padding, H[1] );
	MD5_Append( padding, H[2] );
	MD5_Append( padding, H[3] );
	return padding;
}
String String::FromNumber( int i )
{
	char s[30];
	sprintf( s, "%i", i );
	return s;
}

String operator+( const String & left, const String & right )
{
	String res( left );
	res += right;
	return res;
}
bool operator==( const String & left, const String & right )
{
	return Compare( left, right ) ==0;
}
bool operator!=( const String & left, const String & right )
{
	return Compare( left, right ) !=0;
}




#define HASH_SEED  0xda0
unsigned int MurmurHash2( const void * key, int len, unsigned int seed );

unsigned int MakeHash( const String & s, int length )
{
	if( length <= 0 ) length = s.Size();
	if( length > s.Size() ) length = s.Size();
	return MurmurHash2( s.Data(), length, HASH_SEED );
}
unsigned int MakeHash( const String & s, int offset, int length )
{
	if( offset < 0 ) offset = 0;
	if( length <= 0 ) length = s.Size() - offset;
	if( (length + offset) > s.Size() ) length = s.Size() - offset;
	return MurmurHash2( s.Data() + offset, length, HASH_SEED );
}
int Compare( const String & left, const String & right )
{
	const char *L = left.Data();
	const char *R = right.Data();
	int M = left.Size();
	int N = right.Size();
	int I = 0;
	for(; I<M and I<N and *L==*R; I++, L++, R++);
	if( M == N and I == N ) return 0;
	if( I == M ) return -1;
	if( I == N ) return 1;
	if( *L < *R ) return -1;
	return 1;
#if 0
	const char *s1 = left.Data();
	const char *s2 = right.Data();
	int i1=0, n1 = left.Size();
	int i2=0, n2 = right.Size();
	for(; i1<n1 and i2<n2 and *s1==*s2; i1++, i2++, s1++, s2++);
	if( n1 == n2 and i1 == n1 ) return 0;
	if( i1 == n1 ) return -1;
	if( i2 == n2 ) return 1;
	if( *s1 < *s2 ) return -1;
	return 1;
#endif
}
int ComparePrefix( const String & left, const String & right, int maxmm )
{
	const char *L = left.Data();
	const char *R = right.Data();
	int M = left.Size();
	int N = right.Size();
	int K = M < N ? M : N;
	int I = 0, mm = 0;
	for(; (I<K) & ((mm+=(*L!=*R))<=maxmm); I++, L++, R++);
	return I;
}
int ComparePrefix( const String & left, int lstart, const String & right, int rstart, int max, int maxmm )
{
	const char *L = left.Data() + lstart;
	const char *R = right.Data() + rstart;
	int M = left.Size();
	int N = right.Size();
	int I = lstart;
	int J = rstart;
	int mm = 0;
	if( max ){
		if( lstart + max < M ) M = lstart + max;
		if( rstart + max < N ) N = rstart + max;
	}
	for(; I<M and J<N and (mm+=(*L!=*R))<=maxmm; I++, J++, L++, R++);
	return I - lstart;
}
int CompareSuffix( const String & left, int lstart, const String & right, int rstart, int max, int maxmm )
{
	const char *L = left.Data() + lstart;
	const char *R = right.Data() + rstart;
	int M = 0;
	int N = 0;
	int I = lstart;
	int J = rstart;
	int mm = 0;
	if( max ){
		if( lstart + 1 - max >=0 ) M = lstart + 1 - max;
		if( rstart + 1 - max >=0 ) N = rstart + 1 - max;
	}
	for(; I>=M and J>=N and (mm+=(*L!=*R))<=maxmm; I--, J--, L--, R--);
	return lstart - I;
}
int CountMismatch( const String & left, const String & right, int from, int to )
{
	const char *L = left.Data() + from;
	const char *R = right.Data() + from;
	int M = left.Size();
	int N = right.Size();
	int i, mm = 0;
	if( to < 0 ) to = M;
	if( to > M ) to = M;
	if( to > N ) to = N;
	for(i=from; i<to; i++, L++, R++) mm += *L != *R;
	return mm;
}
int CountMismatch2( const String & left, const String & right, int from, int to, int maxmm )
{
	const char *L = left.Data() + from;
	const char *R = right.Data() + from;
	int M = left.Size();
	int N = right.Size();
	int i, mm = 0;
	if( to < 0 ) to = M;
	if( to > M ) to = M;
	if( to > N ) to = N;
	for(i=from; i<to && mm<=maxmm; i++, L++, R++) mm += *L != *R;
	return mm;
}
int CountMismatch3( const String & left, int lstart, const String & right, int rstart, int max, int maxmm )
{
	const char *L = left.Data() + lstart;
	const char *R = right.Data() + rstart;
	int M = left.Size();
	int N = right.Size();
	int I = lstart;
	int J = rstart;
	int mm = 0;
	if( max ){
		if( lstart + max < M ) M = lstart + max;
		if( rstart + max < N ) N = rstart + max;
	}
	for(; I<M and J<N and (mm+=(*L!=*R))<=maxmm; I++, J++, L++, R++);
	return mm;
}

END_NS_MIN
