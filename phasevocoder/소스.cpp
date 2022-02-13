#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <complex.h>
#include <cstdlib>
#include <windows.h>
#include <omp.h>
#include <iostream>
using namespace std;
class fftmath;
class fftcalc;
class doublearray;
template<typename T> class image;
class wavreader;
class wavwriter;
class bmpwriter;
class vocoder;

class fftmath {
protected:
	int  maxn = 0;
	_Dcomplex* cexpb = 0;
	fftmath(int n) : maxn(n) {
		if (cexpb) free(cexpb);
		cexpb = (_Dcomplex*)malloc(sizeof(_Dcomplex) * maxn);
		calculate(-3.14159265358979323846);
	}
	void calculate(double mpii) {
		_Dcomplex temp;
		for (int i = 0; i < maxn; i++) {
			temp._Val[0] = 0;
			temp._Val[1] = mpii * (static_cast<double>(i) / static_cast<double>(maxn));
			cexpb[i] = cexp(temp);
		}
	}
	virtual ~fftmath() {
		if (cexpb) free(cexpb);
	}
public:
	inline void ifft() {
		calculate(3.14159265358979323846);
	}
	inline void fft() {
		calculate(-3.14159265358979323846);
	}
public:
};
class fftcalc :public fftmath {
private:
	int n;
	int n2;
	_Dcomplex* sw;
	_Dcomplex* outt;
public:
	fftcalc(int n) :n(n), fftmath(n) {
		if (n & (n - 1)) throw exception("n must be power of 2");
		n2 = n >> 1;
		outt = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
	}
	void operator()(_Dcomplex* in, _Dcomplex* out) {
		_Dcomplex* in1 = in, * out1 = outt;
		for (int step = (n >>1); step > 1; step >>=1) {
			for (int bias = 0; bias < step; bias++) {
				for (int i = 0; i < n; i +=  (step<<1)) {
					_Dcomplex& ce1 = fftmath::cexpb[i];
					_Dcomplex& tmp2a = in1[i + step + bias];
					_Dcomplex& tmp2b = in1[i + bias];
					_Dcomplex& tmp3a = out1[(i >> 1) + bias];
					_Dcomplex& tmp3b = out1[((i + n) >> 1) + bias];
					tmp3a._Val[0] = tmp2b._Val[0] + (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
					tmp3a._Val[1] = tmp2b._Val[1] + (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
					tmp3b._Val[0] = tmp2b._Val[0] - (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
					tmp3b._Val[1] = tmp2b._Val[1] - (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
				}
			}
			sw = in1;
			in1 = out1;
			out1 = sw;
		}
		for (int i = 0; i < n; i += 2) {
			_Dcomplex& ce1 = fftmath::cexpb[i];
			_Dcomplex& tmp2a = in1[i + 1];
			_Dcomplex& tmp2b = in1[i];
			_Dcomplex& tmp3a = out[(i >> 1)];
			_Dcomplex& tmp3b = out[((i + n) >> 1)];
			tmp3a._Val[0] = tmp2b._Val[0] + (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
			tmp3a._Val[1] = tmp2b._Val[1] + (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
			tmp3b._Val[0] = tmp2b._Val[0] - (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
			tmp3b._Val[1] = tmp2b._Val[1] - (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
		}
		
	}
	int getoutsize() {
		return n2;
	}
	~fftcalc() {
		free(outt);
	}
};
class doublearray {
public:
	double* arr = 0;
	int memo[4];
	int n;
	void create(int n) {
		this->n = n;
		arr = (double*)malloc(sizeof(double) * n);
	}
	int* getmemo() {
		return memo;
	}
	void remove() {
		if (arr) {
			free(arr);
			arr = 0;
		}
	}
};
template<typename T>
class image {
public:
	int x = 0, y = 0;
	T** data = 0;
	void create(int x, int y) {
		this->x = x;
		this->y = y;
		T* baseaddr = (T*)malloc(sizeof(T) * y * x+ sizeof(T*) * y);
		data = (T**)(baseaddr + y * x);
		for (int i = 0; i < y; i++) data[i] = baseaddr+i*x;
	}
	void remove() {
		if (data != 0) {
			free(data[0]);
			data = 0;
		}
	}
	void writebmp(const char* filename, int writephase = 0) {
		bmpwriter bm(filename, x, y);
		for (int i = 0; i < y; i++) bm(data[i], i, writephase);
	}
};

class wavreader {
private:
	void* data;
	int smpl, nbit, count;
	int ndout = 0;
	doublearray da;
public:
	wavreader(const char* filename) {
		int hdr[11];
		FILE* f = fopen(filename, "rb");
		if (f != 0) {
			fread(hdr, 4, 11, f);
			smpl = hdr[6];
			nbit = hdr[8] >> 19;
			if (hdr[9] != 0x61746164) fread(hdr+5, 4, 6, f);
			count = hdr[10] / nbit;
			data = reinterpret_cast<void*>(malloc(hdr[10]));
			da.create(count);
			da.getmemo()[0] = smpl;
			da.getmemo()[1] = nbit;
			fread(data, nbit, count, f);
			fclose(f);
			double* ptr1 = da.arr, * ptr2 = da.arr + count;
			switch (nbit) {
			case 1:
			{
				unsigned char* tmp0 = (unsigned char*)data;
				unsigned char* tmp2 = tmp0 + count;
				while (ptr1 != ptr2) *(ptr1++) = ((double)(*(tmp0++))) / 128.-1.;
			}
			break;
			case 2:
			{
				short* tmp0 = (short*)data;
				short* tmp2 = tmp0 + count;
				while (ptr1 != ptr2) *(ptr1++) = ((double)(*(tmp0++))) / 32768.;
			}
			break;
			case 3:
			{
				char* tmp0 = (char*)data;
				char* tmp2 = tmp0 + 3*count;
				while (ptr1 != ptr2) {
					long tmp = (*(long*)tmp0)&0x00FFFFFF;
					if (tmp & 0x00800000) tmp |= 0xFF000000;
					tmp0 += 3;
					*(ptr1++) = tmp / 8388608.;
				}
				
			}
			break;
			case 4:
			{
				long* tmp0 = (long*)data;
				long* tmp2 = tmp0 + count;
				while (ptr1 != ptr2) *(ptr1++) = ((double)(*(tmp0++))) / 2147483648.;
			}
			break;
			default:
				throw exception("invaild nbit");
				break;
			}
		}
		else throw exception("file not exist");
	}
	operator doublearray() {
		return da;
	}
	~wavreader() {
		free(data);
		da.remove();
	}
};
class wavwriter {
private:
	unsigned long ret[11];
	const char* filename;
	doublearray dat;
	int smplrate, bitdepth;
public:
	wavwriter(const char* filename, doublearray dat) :filename(filename), dat(dat) {
		smplrate = dat.getmemo()[0];
		bitdepth = dat.getmemo()[1];
		ret[0] = 0x46464952;
		ret[2] = 0x45564157;
		ret[3] = 0x20746D66;
		ret[4] = 0x00000010;
		ret[5] = 0x00010001;
		ret[6] = smplrate;
		ret[7] = ret[6] * bitdepth;
		ret[8] = (bitdepth <<19)| bitdepth;
		ret[9] = 0x61746164;
	}
	~wavwriter() {
		ret[10] = dat.n * bitdepth;
		ret[1] = dat.n * bitdepth + 0x24;
		FILE* f = fopen(filename, "wb");
		switch (bitdepth) {
		case 1: {
			unsigned char* temp = (unsigned char*)malloc(sizeof(unsigned char) * dat.n);
#pragma omp parallel for
			for (int i = 0; i < dat.n; i++) {
				if (dat.arr[i] >= .9921875) temp[i] = 255;
				else if (dat.arr[i] < -1.) temp[i] = 0;
				else temp[i] = (unsigned char)(dat.arr[i] * 128.+128.);
			}
			fwrite(ret, 1, 44, f);
			fwrite(temp, 1, dat.n, f);
			free(temp);
			break;
		}
		case 2: {
			short* temp = (short*)malloc(sizeof(short) * dat.n);
#pragma omp parallel for
			for (int i = 0; i < dat.n; i++) {
				if (dat.arr[i] >= .9999695) temp[i] = 32767;
				else if (dat.arr[i] < -1.) temp[i] = -32768;
				else temp[i] = (short)(dat.arr[i] * 32768.);
			}
			fwrite(ret, 1, 44, f);
			fwrite(temp, 2, dat.n, f);
			free(temp);
			break;
		}
		case 3: {
			union {char d[4];long i;} uni;
			char* temp = (char*)malloc(sizeof(char) * 3 * dat.n);
#pragma omp parallel for
			for (int j = 0; j < dat.n; j++) {
				if (dat.arr[j] >= .999999881) uni.i = 8388607;
				else if (dat.arr[j] < -1.) uni.i = -8388608;
				else uni.i= (long)(dat.arr[j] * 8388608.);
				temp[j+j+j] = uni.d[0]; temp[j+j+j + 1] = uni.d[1]; temp[j+j+j + 2] = uni.d[2];
			}
			fwrite(ret, 1, 44, f);
			fwrite(temp, 3, dat.n, f);
			free(temp);
			break;
		}
		case 4: {
			long* temp = (long*)malloc(sizeof(long) * dat.n);
#pragma omp parallel for
			for (int i = 0; i < dat.n; i++) {
				if (dat.arr[i] >= .999999999534339) temp[i] = 2147483647;
				else if (dat.arr[i] < -1.) temp[i] = -2147483648;
				else temp[i] = (long)(dat.arr[i] * 2147483648.);
			}
			fwrite(ret, 1, 44, f);
			fwrite(temp, 4, dat.n, f);
			free(temp);
			break;
		}
		default:
			throw exception("invaild nbit");
			break;
		}
		
		fclose(f);
		
	}
};
class bmpwriter {
private:
	unsigned char* out = 0;
	const unsigned char red[82] = { 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,20,30,40,50,60,70,80,90,100,110,120,130,122,114,106,98,90,82,74,66,58,50,42,34,26,18,10,2,0 };
	const unsigned char grn[82] = { 255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	const unsigned char blu[82] = { 255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0 };
	const unsigned char b1[12] = { 0,0, 0,0, 54,0,0,0 ,40,0,0,0 };
	const unsigned char b2[28] = { 1,0, 24,0 };
	const unsigned char b3[3] = { 0,0,0 };
	FILE* f = 0;
	int w, h;
public:
	bmpwriter(const char* filename, int width, int height) :w(width), h(height) {
		out = (unsigned char*)malloc(3 * w * h);
		f = fopen(filename, "wb");
		if (f) {
			int filesize = 54 + 3 * w * h;
			fwrite("BM", 1, 2, f);
			fwrite(&filesize, 4, 1, f);
			fwrite(b1, 1, 12, f);
			fwrite(&w, 4, 1, f);
			fwrite(&h, 4, 1, f);
			fwrite(b2, 1, 8, f);
			 filesize = 3 * w * h;
			fwrite(&filesize, 4, 1, f);
			fwrite(b2+4, 1, 16, f);
		}
		else throw exception("file write failed");
	}
	template<int T> struct x { typedef int xx; };
	template<typename T,typename U>//SFINAE 1
	void operator()(T* data, U yindex,int writephase) {
		for (int j = 0; j < w; j++) {
			int t = static_cast<int>(-20. * log10(1e-20 + data[j])) + 50;
			if (t < 0) t = 0;
			if (t > 81) t = 81;
			out[(j + yindex * w) * 3 + 2] = red[t]; out[(j + yindex * w) * 3 + 1] = grn[t]; out[(j + yindex * w) * 3] = blu[t];
		}
	}
	template<typename T,typename U=typename x<sizeof(T::_Val)>::xx>//SFINAE 1
	void operator()(T* data, int yindex, int writephase) {
		if (writephase)
			for (int j = 0; j < w; j++) {
				int bri = 31 + static_cast<int>(20. * log10(1e-20 + cabs(data[j])));
				if (bri < 0) bri = 0;
				if (bri > 81) bri = 81;
				double hue = carg(data[j]) / 6.283185307179586476925286766559 + 0.5;
				double r, g, b;
				if (hue < 0.3333333333333) r = 1. - hue * 3.,g = hue * 3.,b = 0;
				else if (hue < 0.6666666666666666) b = hue * 3 - 1.,g = 2 - hue * 3,r = 0;
				else r = hue * 3 - 2,g = 0,b = 3 - hue * 3;
				out[(j + yindex * w) * 3 + 2] = (unsigned char)(r * (double)bri * 1.5);
				out[(j + yindex * w) * 3 + 1] = (unsigned char)(g * (double)bri * 1.5);
				out[(j + yindex * w) * 3] = (unsigned char)(b * (double)bri * 1.5);
			}
		else
			for (int j = 0; j < w; j++) {
				int t = static_cast<int>(-20. * log10(1e-20 + cabs(data[j]))) + 50;
				if (t < 0) t = 0;
				if (t > 81) t = 81;
				out[(j + yindex * w) * 3 + 2] = red[t]; out[(j + yindex * w) * 3 + 1] = grn[t]; out[(j + yindex * w) * 3] = blu[t];
			}
	}
	~bmpwriter() {
		for (int i = 0; i < h; i++) {
			fwrite(out + (w * (h - i - 1) * 3), 3, w, f);
			fwrite(b3, 1, (4 - ((w * 3) & 3)) & 3, f);
		}
		fclose(f);
		free(out);
	}
};

class vocoder {
public:
	int n, overlap, shift, shift2,orig_len,len;
	double dshift, dshift2;
	double ratio;
	fftcalc* ft32[32];
	double* hann2;
	
	double in_var=0;
	image<_Dcomplex> img;
	
	image<double> vol, iffttmp,centerfreq;
	_Dcomplex* dout[32], * dout2[32];
	doublearray da;
	double* oold, * old, * now;
	
	vocoder(int n, int overlap, double ratio, doublearray input, int corr = 1) {
		init(n, overlap, ratio);
		cout << "initialize done" << endl;
		load(input);
		cout << "fft done" << endl;
		//img.writebmp("x.bmp");
		modifyphase();
		cout << "phase correction done" << endl;
		convert();
		cout << "ifft done" << endl;
		for (int i = 0; i < corr; i++) adjustvolume(2, (i%2)*2), cout << "tunning done" << endl;
	}
	vocoder(double ratio, doublearray input, int corr = 1) :
		vocoder(4096, (ratio > 1.2) ? 10 : 10, ratio, input, corr) {}
	vocoder(double speed, double pitch, doublearray input, int corr = 1) :
		vocoder(pitch /  speed, input, corr) {
		da.memo[0] *= pitch;
	}
	void init(int n, int overlap, double ratio) {
		this->n = n;
		this->ratio = ratio;
		this->overlap = overlap;
		shift2 = n / overlap;
		dshift2 = (double)(shift2);
		dshift = dshift2 / ratio;
		shift = (int)dshift;
		for (int i = 0; i < 32; i++) {
			ft32[i] = new fftcalc(n);
			ft32[i]->fft();
			dout[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
			dout2[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		}
		
		hann2 = (double*)malloc(sizeof(double) * n);
		for (int i = 0; i < n; i++) {
			hann2[i] = sin((double)i / n * 3.1415926535898);
			hann2[i] *= hann2[i]; hann2[i] *= hann2[i];
		}
		oold = (double*)malloc(sizeof(double) * n / 2 + 10);
		old = (double*)malloc(sizeof(double) * n / 2 + 10);
		now = (double*)malloc(sizeof(double) * n / 2 + 10);
	}
	void load(doublearray data) {
		int y = (data.n + shift - 1) / shift;
		orig_len = data.n;
		img.create(n, y);
		iffttmp.create(n, y);
		centerfreq.create(n, y);
		da.create((int)((double)y * dshift2) + n);
		
		for (int i = 0; i < orig_len; i++) in_var += data.arr[i] * data.arr[i];
		in_var /= orig_len;
		for (int i = 0; i < 4; i++) da.getmemo()[i] = data.getmemo()[i];
		for (int i = 0; i < 32; i++) ft32[i]->fft();
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < y; i++) {
			int ii,ii2;
			int threadid = omp_get_thread_num();
			ii2 = (int)((double)i * dshift) - n;
			for (int j = 0; j < n; j++) {
				ii = ii2 + j;
				if (ii >= data.n || ii < 0) dout[threadid][j]._Val[0] = 0;
				else dout[threadid][j]._Val[0] = data.arr[ii] *hann2[j];
				dout[threadid][j]._Val[1] = 0;
			}
			(*ft32[threadid])(dout[threadid], img.data[i]);
		}
	}
	image<_Dcomplex>& getimage() {
		return img;
	}
	bool isnear(double upper, double _centerfreq) {
		if (_centerfreq < 100) {
			if (_centerfreq < 50) {
				if (upper - _centerfreq < 4 && _centerfreq - upper < 4) return true;
			}
			else {
				if (upper - _centerfreq < 7 && _centerfreq - upper < 7) return true;
			}
		}
		else {
			if (_centerfreq < 200) {
				if (upper - _centerfreq < 12 && _centerfreq - upper < 12)return true;
			}
			else {
				if (upper - _centerfreq < 24 && _centerfreq - upper < 24) return true;
			}
		}
		return false;
	}
	int phasewidth(int centl) {
		if (centl < 25) {
			if (centl < 10) {
				if (centl < 6) return 1;
				else return 2;
			}
			else {
				if (centl < 15) return 3;
				else return 4;
			}
		}
		else {
			if (centl < 100) {
				if (centl < 50) return 6;
				else return 9;
			}
			else {
				if (centl < 200) return 12;
				else return 16;
			}
		}
	}
	int dist(int a, int b) {
		if (a > b)return a - b;
		else return b - a;
	}
	void modifyphase() {
		image<double> tt;
		tt.create(img.x, img.y);
#pragma omp parallel for num_threads(32)
		for (int a = 1; a < img.y; a++) {
			for (int b = 0; b < n; b++) {
				iffttmp.data[a][b] = cabs(img.data[a][b]);
				tt.data[a][b] = iffttmp.data[a][b];
			}
			int le, cent;
			int ri = 0;
			double _centerfreq, _mag;
			double peak;
			int continuity;
			int cut = 0;
			int ri_backup;
			
			while (1) {
				le = ri;
				ri = le;
				cent = le;
				peak = iffttmp.data[a][le];
				continuity = 1;
				cut = 0;
				while (1) {
					if (ri > n / 2) break;
					if (peak < iffttmp.data[a][ri] && !cut) {
						peak = iffttmp.data[a][ri]; cent = ri;
					}
					if (!ri) ri++;
					else if (iffttmp.data[a][ri] >= iffttmp.data[a][ri + 1]) {
						ri_backup = ri + 1;
						
							if (ri < 50 && iffttmp.data[a][ri] < 0.8 * peak) cut = 1;
							else if (ri < 100 && iffttmp.data[a][ri] < 0.7 * peak) cut = 1;
							else if (ri < 150 && iffttmp.data[a][ri] < 0.7 * peak) cut = 1;
							else if (ri < 200 && iffttmp.data[a][ri] < 0.65 * peak) cut = 1;
							else if (iffttmp.data[a][ri] < 0.65 * peak) cut = 1;
						
					}
					else {
						if (cut && iffttmp.data[a][ri_backup] * 1.1 < iffttmp.data[a][ri]) {
							ri = ri_backup;
							break;
						}
					}
					ri++;
				}
				int centl, centr;
				
				for (centl = le; iffttmp.data[a][centl] < peak * 0.7; centl++);
				for (centr = ri; iffttmp.data[a][centr] < peak * 0.7; centr--);
				cent = (cent+((centl + centr) >> 1))>>1;
				_centerfreq = iffttmp.data[a][cent] * (double)cent;
				_mag = iffttmp.data[a][cent];
				int p = 1;
				while (cent - p >= le && cent + p <= ri && p < 4) {
					_centerfreq += iffttmp.data[a][cent + p] * (double)(cent + p) + iffttmp.data[a][cent - p] * (double)(cent - p);
					_mag += iffttmp.data[a][cent + p] + iffttmp.data[a][cent - p];
					p++;
				}
				if (_mag) _centerfreq /= _mag;
				else _centerfreq = (double)((le + ri) >> 1);
				for (int i = le; i <= ri; i++) centerfreq.data[a][i] = _centerfreq;
				tt.data[a][(int)(_centerfreq+.5)] = 1e9;
				if (ri >= n / 2) break;
			}
		}
	
		double dn = dshift2 * 3.14159265358979323846 /(double)n;
#pragma omp parallel for num_threads(32)
		for (int b = 0; b <= n / 2; b++) {
			double base = carg(img.data[0][b]);
			for (int a = 1; a < img.y; a++) {
				base += dn * (centerfreq.data[a - 1][b] + centerfreq.data[a][b]);
				img.data[a][b] = cexp(_Dcomplex{ 0,base });
				if (!(a & 15)) base = carg(img.data[a][b]);
			}
		}
#pragma omp parallel for num_threads(32)
		for (int a = 0; a < img.y; a++) {
			int b = 0;
			for (; b <= n / 2; b++) {
				img.data[a][b]._Val[0] *= iffttmp.data[a][b];
				img.data[a][b]._Val[1] *= iffttmp.data[a][b];
			}
			for (; b < n; b++) {
				img.data[a][b]._Val[0] = img.data[a][n - b]._Val[0];
				img.data[a][b]._Val[1] = -img.data[a][n - b]._Val[1];
			}
		}
		tt.writebmp("a.bmp");
	}
	
	void modifyphase2() {
		double v1 = 0, v2 = 1;
		if (ratio > 1) {
			//v1 = 0.95 - cos(ratio * 3.1415926535897932384626433832795) * 0.05;
			//v2 = 0.05 + cos(ratio * 3.1415926535897932384626433832795) * 0.05;
		}
#pragma omp parallel for num_threads(32)
		for (int b = 0; b <= n / 2; b++) {
			double op = carg(img.data[1][b]);
			for (int a = 2; a < img.y; a++) {
				double base = carg(img.data[a - 1][b]);
				_Dcomplex& imgdataab = img.data[a][b];
				double arg1 = carg(imgdataab);
				double darg = arg1 - op;
				double tmp2 = dshift2 / ratio * 6.283185307179586476925286766559 * ((double)b) / (double)n;
				double darg2 = (double)floor((tmp2 - darg) / 6.283185307179586476925286766559 + 0.5) * 6.283185307179586476925286766559;
				darg += darg2;
				double tmp = ratio * (darg * v1 + tmp2 * v2);
				op = arg1;
				double mag = cabs(imgdataab);
				imgdataab = cexp(_Dcomplex{ 0,(base + tmp) });
				imgdataab._Val[0] *= mag;
				imgdataab._Val[1] *= mag;
			}
		}
		for (int b = n / 2 + 1; b < n; b++)
			for (int a = 1; a < img.y; a++) {
				img.data[a][b]._Val[0] = img.data[a][n - b]._Val[0];
				img.data[a][b]._Val[1] = -img.data[a][n - b]._Val[1];
			}
	}
	void convert() {
		for (int i = 0; i < 32; i++) ft32[i]->ifft();
		len = (int)((double)orig_len * ratio);
		
		for (int i = 0; i < len; i++) da.arr[i] = 0;
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < img.y; i+=1) {
			int threadid = omp_get_thread_num();
			for (int j = 0; j < n; j++) dout2[threadid][j] = img.data[i][j];
			(*ft32[threadid])(dout2[threadid], dout[threadid]);
			int tmp = (int)(dshift2 * (double)i) - n;
			for (int j = 0; j < n; j++) {
				if (tmp < len) {
					if (tmp >= 0) 
#pragma omp atomic
						da.arr[tmp] += hann2[j] * dout[threadid][j]._Val[0];
				}
				else break;
				tmp++;
			}
		}
		double out_var = 0;
		for (int i = 0; i < len; i++) out_var += da.arr[i] * da.arr[i];
		out_var = sqrt(in_var * len / out_var);
		for (int i = 0; i < len; i++) da.arr[i] *= out_var;
	}
	void convert_old() {
		for (int i = 0; i < 32; i++) ft32[i]->ifft();
		len = (int)((double)orig_len*ratio);
		if (iffttmp.data == 0) iffttmp.create(img.x, img.y);
		for (int i = 0; i < len; i++) da.arr[i] = 0;
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < img.y; i++) {
			int threadid = omp_get_thread_num();
			for (int j = 0; j < n; j++) dout2[threadid][j] = img.data[i][j];
			(*ft32[threadid])(dout2[threadid], dout[threadid]);
			for (int j = 0; j < n; j++) iffttmp.data[i][j] = dout[threadid][j]._Val[0];
		}
		for (int i = 0; i < img.y; i++) {
			int ovlshift = 0;
			
			if (i != 0) {
				int shif = shift2 * i - n;
				double pearson_max = -1;
				int imax = -1;
#pragma omp parallel for num_threads(32)
				for (int k = -(shift2 >> 3); k < (shift2 >> 3); k++) {
					double var1 = 0, cov = 0, var2 = 0;
					for (int m = (n >> 2); m < (n >> 1); m++) {
						double t1, t2;
						if (m + k + shif < len && m + k + shif >= 0)
							t1 = da.arr[m + k + shif];
						else t1 = 0;
						t2 = iffttmp.data[i][m];
						var1 += t1 * t1;
						var2 += t2 * t2;
						cov += t1 * t2;
					}
					double pearson = cov * cov / (1e-12 + var1 * var2);
					if (k > 0) pearson -= (double)k / 6;
					else pearson += (double)k / 6;
					if (pearson > pearson_max) {
#pragma omp critical
						{
							pearson_max = pearson;
							imax = k;
						}
					}
				}
				ovlshift = imax;
			}
			
			int tmp = shift2 * i + ovlshift - n;
			for (int j = 0; j < n; j++) {
				if (tmp < len) {
					if (tmp >= 0) da.arr[tmp] += hann2[j] * iffttmp.data[i][j];
				}
				else break;
				tmp++;
			}
		}
		double out_var=0;
		for (int i = 0; i < len; i++) out_var += da.arr[i]*da.arr[i];
		out_var = sqrt(in_var * len / out_var);
		for (int i = 0; i < len; i++) da.arr[i] *= out_var;
	}

	void adjustvolume(int yres,int yshif) {//unused
		
		int shift3 = 2*n / overlap;
		int yf = 1 >> (yres - 1);
		if (vol.data) {
			if (vol.x - 2 != n  || vol.y - 2 != img.y >> (yres+1)) {
				vol.remove();
				vol.create(n + 2, (img.y >> (yres+1)) + 2);
			}
		}
		else vol.create(n  + 2, (img.y >> (yres+1)) + 2);
		for (int i = 0; i < 32; i++) ft32[i]->fft();
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < (img.y>>1); i++) {
			int ii;
			int threadid = omp_get_thread_num();
			for (int j = 0; j < n; j++) {
				ii = i * shift3 + j - n;
				if (ii >= da.n || ii < 0) dout[threadid][j]._Val[0] = 0;
				else dout[threadid][j]._Val[0] = da.arr[ii] * hann2[j];
				dout[threadid][j]._Val[1] = 0;
			}
			(*ft32[threadid])(dout[threadid], img.data[i]);
		}
#pragma omp parallel for num_threads(32)
		for (int j = 0; j < vol.y - 1; j++) {
			for (int i = 0; i < vol.x - 1; i++) {
				double orig_max = 0, fixed_max = 0;
				for (int m = -yf; m < yf; m++) {
					int b = (j << yres) + m + yshif;
					if (b >= 0 && b < img.y) {
						double temp = cabs(img.data[b][i]);
						if (temp > fixed_max)
							fixed_max = temp;
						if (iffttmp.data[b][i] > orig_max)
							orig_max = iffttmp.data[b<<1][i];
					}
				}
				vol.data[j][i] = ((1e-8 + orig_max) / (1e-8 + fixed_max));
			}
			
		}
		for (int i = 0; i < vol.x - 1; i++)
			vol.data[vol.y - 1][i] = vol.data[vol.y - 2][i];
		for (int i = 0; i < vol.y; i++)
			vol.data[i][vol.x - 1] = vol.data[i][vol.x - 2];
		for (int i = 0; i < 32; i++) ft32[i]->ifft();
		
		for (int i = 0; i < len; i++) da.arr[i] = 0;
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < (img.y>>1); i++) {
			int threadid = omp_get_thread_num();
			int vy = i >> 3;
			double ry = (double)i / 8. - (double)vy;
			for (int j = 0; j < n; j++) {
				int vx = j >> 2;
				double rx = (double)j / 4. - (double)vx;
				double vl = rx * (ry * vol.data[vy + 1][vx + 1] + (1. - ry) * vol.data[vy][vx + 1]) +
					(1. - rx) * (ry * vol.data[vy + 1][vx] + (1. - ry) * vol.data[vy][vx]);
				_Dcomplex& dout2ref = dout2[threadid][j];
				_Dcomplex& imgdataref = img.data[i][j];
				dout2ref._Val[0] = imgdataref._Val[0] * vl;
				dout2ref._Val[1] = imgdataref._Val[1] * vl;
			}
			(*ft32[threadid])(dout2[threadid], dout[threadid]);
#pragma omp critical			
			for (int j = 0; j < n; j++) {
				int tmp = shift3 * i - n + j;
				if (tmp < len && tmp >= 0) da.arr[tmp] += hann2[j] * dout[threadid][j]._Val[0];
			}
		}
		double out_var = 0;
		for (int i = 0; i < len; i++) out_var += da.arr[i] * da.arr[i];
		out_var = sqrt(in_var * len / out_var);
		for (int i = 0; i < len; i++) da.arr[i] *= out_var;
	}
	virtual ~vocoder() {
		
		for (int i = 0; i < 32; i++) {
			delete ft32[i];
			free(dout[i]);
			free(dout2[i]);
		}
		free(oold); free(old); free(now);
		free(hann2);
		img.remove();
		vol.remove();
		iffttmp.remove();
		centerfreq.remove();
		da.remove();
	}
	operator doublearray() const {
		return da;
	}
};

void main(int argc, char* argv[]) {
	unsigned int t = GetTickCount();
	wavwriter(argv[2], vocoder(atof(argv[3]), atof(argv[4]), wavreader(argv[1]), 0));
	printf("%d", GetTickCount() - t);
}
