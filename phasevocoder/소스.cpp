#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <complex.h>
#include <cstdlib>
#include <windows.h>
#include <omp.h>
#include <iostream>
#define debug
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
	double lev(image<T> ref) {
		//writebmp("a.bmp");
		//ref.writebmp("b.bmp");
		//exit(0);
		double mid = 0;
		double cent = 0;
		double s1 = 0,s2=0;
#pragma omp parallel for reduction(+:mid)
		for (int i = 2; i < y - 40; i++) {
			for (int j = 0; j < 205; j++) {
				s1 += data[i][j];
				s2+= ref.data[i+9][j];
			}
		}
		cent = s1 / s2;
		for (int i = 2; i < y - 40; i++) {
			for (int j = 0; j < 250; j++) {
				double temp = (1e-9+data[i][j]) / (1e-9+ref.data[i+9][j]*cent);
				if (temp < 1) temp = 1 / temp;
				//if (temp < 1.2) temp = 1.2;
				mid += (temp - 1)/(10+j)*pow(data[i][j],0.2);
			}
		}
		return mid;
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
int _a[13] = { 0, };
double _d[12] = { -0.167782,0.165182,0.421332,0.184772,0.489644,0.933740,0.666953,0.476673,0.832487,0.803771,0.708964 ,0.1};
class vocoder {
public:
	int n, overlap, shift, shift2, orig_len, len;
	double dshift, dshift2;
	double ratio;
	fftcalc* ft32[32];
	double* hann2;

	double in_var = 0;
	image<_Dcomplex> img;

	image<double> vol, iffttmp, centerfreq;
	_Dcomplex* dout[32], * dout2[32];
	doublearray da;
	double* oold, * old, * now;

	vocoder(int n, int overlap, double ratio, doublearray input) {
		unsigned int t = GetTickCount();
		init(n, overlap, ratio);
		//cout << "initialize done" << endl;
		load(input);
		//cout << "fft done" << endl;
		//img.writebmp("x.bmp");
		modifyphase();
		//cout << "phase correction done" << endl;
		convert();
		//cout << "ifft done" << endl;
		cout << GetTickCount() - t << endl;
		
	}
	vocoder(double ratio, doublearray input) :
		vocoder(4096, 20, ratio, input) {}
	vocoder(double speed, double pitch, doublearray input) :
		vocoder(pitch / speed, input) {
		da.memo[0] *= pitch;
	}
	void init(int n, int overlap, double ratio) {
		this->n = n;
		this->ratio = ratio;
		this->overlap = overlap;
		shift2 = (int)((double)n *ratio/ (double)overlap);
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
			hann2[i] *= hann2[i]*hann2[i]; //hann2[i] *= hann2[i];
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
			int ii, ii2;
			int threadid = omp_get_thread_num();
			ii2 = (int)((double)i * dshift) - n;
			for (int j = 0; j < n; j++) {
				ii = ii2 + j;
				if (ii >= data.n || ii < 0) dout[threadid][j]._Val[0] = 0;
				else dout[threadid][j]._Val[0] = data.arr[ii] * hann2[j];
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
				if (upper - _centerfreq < 4 + _a[1] && _centerfreq - upper < 4 + _a[1]) return true;
			}
			else {
				if (upper - _centerfreq < 7 + _a[2] && _centerfreq - upper < 7 + _a[2]) return true;
			}
		}
		else {
			if (_centerfreq < 200) {
				if (upper - _centerfreq < 12 + _a[3] && _centerfreq - upper < 12 + _a[3])return true;
			}
			else {
				if (upper - _centerfreq < 24 + _a[4] && _centerfreq - upper < 24 + _a[4]) return true;
			}
		}
		return false;
	}
	void modifyphase() {
		//image<double> tt;
		//tt.create(img.x, img.y);
#pragma omp parallel for num_threads(32)
		for (int a = 1; a < img.y; a++) {
			for (int b = 0; b < n; b++) {
				iffttmp.data[a][b] = cabs(img.data[a][b]);
				//tt.data[a][b] = iffttmp.data[a][b];
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

						if (ri < 50 && iffttmp.data[a][ri] <_d[5] * peak) cut = 1;
						else if (ri < 100 && iffttmp.data[a][ri] < _d[6] * peak) cut = 1;
						else if (ri < 150 && iffttmp.data[a][ri] < _d[7]* peak) cut = 1;
						else if (ri < 200 && iffttmp.data[a][ri] < _d[8] * peak) cut = 1;
						else if (iffttmp.data[a][ri] < _d[9] * peak) cut = 1;

					}
					else {
						if (cut && iffttmp.data[a][ri_backup] * (1+_d[11]) < iffttmp.data[a][ri]) {
							ri = ri_backup;
							break;
						}
					}
					ri++;
				}
				int centl, centr;

				for (centl = le; iffttmp.data[a][centl] < peak * _d[10]; centl++);
				for (centr = ri; iffttmp.data[a][centr] < peak * _d[10]; centr--);
				cent = (cent + ((centl + centr) >> 1)) >> 1;
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
				
				if (ri >= n / 2) break;
			}
		}
#pragma omp parallel for
		for (int b = 0; b <= n / 2; b++) {
			old[b] = centerfreq.data[0][b];
			now[b] = centerfreq.data[1][b];
		}
		for (int a = 2; a < img.y - 2; a++) {
			double* tmp;
			tmp = oold;
			oold = old; old = now; now = tmp;
			for (int b = 0; b <= n / 2; b++) {
				now[b] = centerfreq.data[a][b];
			}
			for (int b = 0; b < n / 2;) {
				double centerfreq2o, centerfreq0, centerfreq1, centerfreq2, centerfreq3, centerfreq4;
				int centerfreq_int;
				centerfreq2o = now[b];
				centerfreq2 = centerfreq2o;
				centerfreq_int = (int)(centerfreq2 + .5);
			
				centerfreq0 = oold[centerfreq_int]; centerfreq1 = old[centerfreq_int]; centerfreq3 = centerfreq.data[a + 1][centerfreq_int]; centerfreq4 = centerfreq.data[a + 2][centerfreq_int];
				
				if (isnear(centerfreq1, centerfreq2)) {
					if (isnear(centerfreq2, centerfreq3)) {
						centerfreq2 = _d[0] * centerfreq0 + _d[1] * centerfreq1 + (1-_d[0]-_d[1]) * centerfreq2 + _d[1] * centerfreq3 + _d[0] * centerfreq4;

						centerfreq_int = (int)(centerfreq2 + .5);// tt.data[a][centerfreq_int] = 1e9;
					
						
						for (; centerfreq.data[a][b] == centerfreq2o; b++) centerfreq.data[a][b] = centerfreq2;
					}
					else {
						centerfreq2 = _d[3] * centerfreq0 + _d[4] * centerfreq1 + (1-_d[3]-_d[4]) * centerfreq2;
						centerfreq_int = (int)(centerfreq2 + .5); //tt.data[a][centerfreq_int] = 1e9;
						for (; centerfreq.data[a][b] == centerfreq2o; b++) centerfreq.data[a][b] = centerfreq2;
					}
				}
				else {
					if (isnear(centerfreq2, centerfreq3)) {
						centerfreq2 = (1 - _d[3] - _d[4]) * centerfreq2 + _d[4] * centerfreq3 + _d[3] * centerfreq4;
						centerfreq_int = (int)(centerfreq2 + .5);//tt.data[a][centerfreq_int] = 1e9;
						for (; centerfreq.data[a][b] == centerfreq2o; b++) centerfreq.data[a][b] = centerfreq2;
					}
					else {
						//tt.data[a][centerfreq_int] = 1e9;
						for (; centerfreq.data[a][b] == centerfreq2o; b++);

					}
				}
			}
		}
		double dn = (dshift - dshift2) * 2*3.14159265358979323846 / (double)n;
#pragma omp parallel for num_threads(32)
		for (int b = 0; b <= n / 2; b++) {
			double op = carg(img.data[1][b]);
			double base = op;
			for (int a = 2; a < img.y; a++) {

				if (!(a & 15))
					base = carg(img.data[a - 1][b]);

				double arg1 = carg(img.data[a][b]);
				double darg = arg1 - op;
				darg -= dn * centerfreq.data[a][b];

				op = arg1;
				base += darg;
				img.data[a][b] = cexp(_Dcomplex{ 0,base });


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
		//tt.writebmp("a.bmp");
		//tt.remove();
	}

	void modifyphase2() {
		double v1 = 0, v2 = 1;
		
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
		for (int i = 0; i < img.y; i += 1) {
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
		len = (int)((double)orig_len * ratio);
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
		double out_var = 0;
		for (int i = 0; i < len; i++) out_var += da.arr[i] * da.arr[i];
		out_var = sqrt(in_var * len / out_var);
		for (int i = 0; i < len; i++) da.arr[i] *= out_var;
	}
	double getcurr(image<double> ref) {
		int y = (da.n + shift2 - 1) / shift2;
		img.remove();
		img.create(n, y);
		iffttmp.remove();
		iffttmp.create(n, y);
		for (int i = 0; i < 32; i++) ft32[i]->fft();
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < y; i++) {
			int ii;
			int threadid = omp_get_thread_num();
			for (int j = 0; j < n; j++) {
				ii = i * shift2 + j - n;
				if (ii >= da.n || ii < 0) dout[threadid][j]._Val[0] = 0;
				else dout[threadid][j]._Val[0] = da.arr[ii] * hann2[j];
				dout[threadid][j]._Val[1] = 0;
			}
			(*ft32[threadid])(dout[threadid], img.data[i]);
		}
#ifdef debug
		//img.writebmp("k.bmp");
#endif
#pragma omp parallel for num_threads(32)
		for (int a = 0; a < y; a++) {
			for (int b = 0; b < n; b++) {
				iffttmp.data[a][b] = cabs(img.data[a][b]);
			}
		}
		return iffttmp.lev(ref);

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
class vocoder_mtlb {
public:
	int n, overlap, shift2, orig_len, len;
	double dshift2;
	double ratio;
	fftcalc* ft32[32];
	double* hann2;

	double in_var = 0;
	image<_Dcomplex> img1, img2;

	image<double> iffttmp1;
	_Dcomplex* dout[32], * dout2[32];
	doublearray da;
	double* oold, * old, * now;

	vocoder_mtlb(int n, int overlap, double ratio, doublearray input) {
		unsigned int t = GetTickCount();
		init(n, overlap, ratio);
		//cout << "initialize done" << endl;
		load(input);
		//cout << "fft done" << endl;
		//img.writebmp("x.bmp");
		modifyphase();
		//cout << "phase correction done" << endl;
		convert();
		//cout << "ifft done" << endl;
		cout << GetTickCount() - t << endl;
	}
	vocoder_mtlb(double ratio, doublearray input) :
		vocoder_mtlb(1024, 4, ratio, input) {}
	vocoder_mtlb(double speed, double pitch, doublearray input) :
		vocoder_mtlb(pitch / speed, input) {
		da.memo[0] *= pitch;
	}
	void init(int n, int overlap, double ratio) {
		this->n = n;
		this->ratio = ratio;
		this->overlap = overlap;
		shift2 = n / overlap;
		dshift2 = (double)(shift2);

		for (int i = 0; i < 32; i++) {
			ft32[i] = new fftcalc(n);
			ft32[i]->fft();
			dout[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
			dout2[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		}

		hann2 = (double*)malloc(sizeof(double) * n);
		for (int i = 0; i < n; i++) {
			hann2[i] = sin((double)i / n * 3.1415926535898);
			hann2[i] *= hann2[i];// *hann2[i]; //hann2[i] *= hann2[i];
		}

	}
	void load(doublearray data) {
		int y = (data.n + shift2 - 1) / shift2;
		orig_len = data.n;
		img1.create(n, y);
		iffttmp1.create(n, y);

		da.create((int)((double)y * dshift2 * ratio) + n);

		for (int i = 0; i < orig_len; i++) in_var += data.arr[i] * data.arr[i];
		in_var /= orig_len;
		for (int i = 0; i < 4; i++) da.getmemo()[i] = data.getmemo()[i];
		for (int i = 0; i < 32; i++) ft32[i]->fft();
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < y; i++) {
			int ii, ii2;
			int threadid = omp_get_thread_num();
			ii2 = (int)((double)i * dshift2) - n;
			for (int j = 0; j < n; j++) {
				ii = ii2 + j;
				if (ii >= data.n || ii < 0) dout[threadid][j]._Val[0] = 0;
				else dout[threadid][j]._Val[0] = data.arr[ii] * hann2[j];
				dout[threadid][j]._Val[1] = 0;
			}
			(*ft32[threadid])(dout[threadid], img1.data[i]);
		}
	}
	image<_Dcomplex>& getimage() {
		return img1;
	}
	void modifyphase() {
		int x = img1.x;
		int y1 = img1.y;
		int y2 = (int)((double)y1 * ratio);
		img2.create(x, y2);
#pragma omp parallel for num_threads(32)
		for (int a = 0; a < img1.y; a++) {
			for (int b = 0; b < n; b++) {
				iffttmp1.data[a][b] = cabs(img1.data[a][b]);
			}
		}
		double dn = (dshift2) * 2 * 3.14159265358979323846 / (double)n;

		double dratio = 1 / ratio;
#pragma omp parallel for num_threads(32)
		for (int b = 0; b <= n / 2; b++) {
			int neut = 0;
			double prec = 0;
			double phas = 0;
			neut = 0; prec = 0;
			double ph1 = carg(img1.data[0][b]);
			double ph2 = carg(img1.data[1][b]);
			phas = ph1;
			for (int a = 0; a < img2.y; a++) {

				while (prec > 1) {
					prec -= 1; if (neut < img1.y - 2) neut++;
					ph1 = ph2; ph2 = carg(img1.data[neut + 1][b]);
				}
				double adv = ph2 - ph1 - dn * (double)b;
				while (adv > 3.14159265358979323846) adv -= 2 * 3.14159265358979323846;

				double mag = ((1. - prec) * iffttmp1.data[neut][b] + prec * iffttmp1.data[neut + 1][b]);
				img2.data[a][b] = cexp(_Dcomplex{ 0,phas });
				img2.data[a][b]._Val[0] *= mag;
				img2.data[a][b]._Val[1] *= mag;
				phas += adv + dn * (double)b;
				prec += dratio;
			}
		}
#pragma omp parallel for num_threads(32)
		for (int a = 0; a < img2.y; a++) {
			int b = n / 2 + 1;

			for (; b < n; b++) {
				img2.data[a][b]._Val[0] = img2.data[a][n - b]._Val[0];
				img2.data[a][b]._Val[1] = -img2.data[a][n - b]._Val[1];
			}
		}
		//img1.writebmp("2.bmp");
		//img2.writebmp("a.bmp");
		//tt.remove();
	}


	void convert() {

		for (int i = 0; i < 32; i++) ft32[i]->ifft();
		len = (int)((double)orig_len * ratio);

		for (int i = 0; i < len; i++) da.arr[i] = 0;
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < img2.y; i += 1) {
			int threadid = omp_get_thread_num();
			for (int j = 0; j < n; j++) dout2[threadid][j] = img2.data[i][j];
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


	virtual ~vocoder_mtlb() {

		for (int i = 0; i < 32; i++) {
			delete ft32[i];
			free(dout[i]);
			free(dout2[i]);
		}
		free(hann2);
		img1.remove();
		img2.remove();
		iffttmp1.remove();

		da.remove();
	}
	operator doublearray() const {
		return da;
	}
};
class vocoder_mtlb_sin3 {
public:
	int n, overlap, shift2, orig_len, len;
	double dshift2;
	double ratio;
	fftcalc* ft32[32];
	double* hann2;

	double in_var = 0;
	image<_Dcomplex> img1, img2;

	image<double> iffttmp1;
	_Dcomplex* dout[32], * dout2[32];
	doublearray da;
	double* oold, * old, * now;

	vocoder_mtlb_sin3(int n, int overlap, double ratio, doublearray input) {
		unsigned int t = GetTickCount();
		init(n, overlap, ratio);
		//cout << "initialize done" << endl;
		load(input);
		//cout << "fft done" << endl;
		//img.writebmp("x.bmp");
		modifyphase();
		//cout << "phase correction done" << endl;
		convert();
		//cout << "ifft done" << endl;
		cout << GetTickCount() - t << endl;
	}
	vocoder_mtlb_sin3(double ratio, doublearray input) :
		vocoder_mtlb_sin3(4096, 20, ratio, input) {}
	vocoder_mtlb_sin3(double speed, double pitch, doublearray input) :
		vocoder_mtlb_sin3(pitch / speed, input) {
		da.memo[0] *= pitch;
	}
	void init(int n, int overlap, double ratio) {
		this->n = n;
		this->ratio = ratio;
		this->overlap = overlap;
		shift2 = n / overlap;
		dshift2 = (double)(shift2);

		for (int i = 0; i < 32; i++) {
			ft32[i] = new fftcalc(n);
			ft32[i]->fft();
			dout[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
			dout2[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		}

		hann2 = (double*)malloc(sizeof(double) * n);
		for (int i = 0; i < n; i++) {
			hann2[i] = sin((double)i / n * 3.1415926535898);
			hann2[i] *= hann2[i] *hann2[i]; //hann2[i] *= hann2[i];
		}

	}
	void load(doublearray data) {
		int y = (data.n + shift2 - 1) / shift2;
		orig_len = data.n;
		img1.create(n, y);
		iffttmp1.create(n, y);

		da.create((int)((double)y * dshift2 * ratio) + n);

		for (int i = 0; i < orig_len; i++) in_var += data.arr[i] * data.arr[i];
		in_var /= orig_len;
		for (int i = 0; i < 4; i++) da.getmemo()[i] = data.getmemo()[i];
		for (int i = 0; i < 32; i++) ft32[i]->fft();
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < y; i++) {
			int ii, ii2;
			int threadid = omp_get_thread_num();
			ii2 = (int)((double)i * dshift2) - n;
			for (int j = 0; j < n; j++) {
				ii = ii2 + j;
				if (ii >= data.n || ii < 0) dout[threadid][j]._Val[0] = 0;
				else dout[threadid][j]._Val[0] = data.arr[ii] * hann2[j];
				dout[threadid][j]._Val[1] = 0;
			}
			(*ft32[threadid])(dout[threadid], img1.data[i]);
		}
	}
	image<_Dcomplex>& getimage() {
		return img1;
	}
	void modifyphase() {
		int x = img1.x;
		int y1 = img1.y;
		int y2 = (int)((double)y1 * ratio);
		img2.create(x, y2);
#pragma omp parallel for num_threads(32)
		for (int a = 0; a < img1.y; a++) {
			for (int b = 0; b < n; b++) {
				iffttmp1.data[a][b] = cabs(img1.data[a][b]);
			}
		}
		double dn = (dshift2) * 2 * 3.14159265358979323846 / (double)n;

		double dratio = 1 / ratio;
#pragma omp parallel for num_threads(32)
		for (int b = 0; b <= n / 2; b++) {
			int neut = 0;
			double prec = 0;
			double phas = 0;
			neut = 0; prec = 0;
			double ph1 = carg(img1.data[0][b]);
			double ph2 = carg(img1.data[1][b]);
			phas = ph1;
			for (int a = 0; a < img2.y; a++) {

				while (prec > 1) {
					prec -= 1; if (neut < img1.y - 2) neut++;
					ph1 = ph2; ph2 = carg(img1.data[neut + 1][b]);
				}
				double adv = ph2 - ph1 - dn * (double)b;
				while (adv > 3.14159265358979323846) adv -= 2 * 3.14159265358979323846;

				double mag = ((1. - prec) * iffttmp1.data[neut][b] + prec * iffttmp1.data[neut + 1][b]);
				img2.data[a][b] = cexp(_Dcomplex{ 0,phas });
				img2.data[a][b]._Val[0] *= mag;
				img2.data[a][b]._Val[1] *= mag;
				phas += adv + dn * (double)b;
				prec += dratio;
			}
		}
#pragma omp parallel for num_threads(32)
		for (int a = 0; a < img2.y; a++) {
			int b = n / 2 + 1;

			for (; b < n; b++) {
				img2.data[a][b]._Val[0] = img2.data[a][n - b]._Val[0];
				img2.data[a][b]._Val[1] = -img2.data[a][n - b]._Val[1];
			}
		}
		//img1.writebmp("2.bmp");
		//img2.writebmp("a.bmp");
		//tt.remove();
	}


	void convert() {

		for (int i = 0; i < 32; i++) ft32[i]->ifft();
		len = (int)((double)orig_len * ratio);

		for (int i = 0; i < len; i++) da.arr[i] = 0;
#pragma omp parallel for num_threads(32)
		for (int i = 0; i < img2.y; i += 1) {
			int threadid = omp_get_thread_num();
			for (int j = 0; j < n; j++) dout2[threadid][j] = img2.data[i][j];
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


	virtual ~vocoder_mtlb_sin3() {

		for (int i = 0; i < 32; i++) {
			delete ft32[i];
			free(dout[i]);
			free(dout2[i]);
		}
		free(hann2);
		img1.remove();
		img2.remove();
		iffttmp1.remove();

		da.remove();
	}
	operator doublearray() const {
		return da;
	}
};
void mai2n(int argc, char* argv[]) {
	double hi = 0;
	double fir, dec;
	int _aa[13];
	double _dd[12];
	int j = 0, k = 0;
	wavreader wr = wavreader("src.wav");
	wavreader wr2 = wavreader("targ.wav");
	vocoder v = vocoder(.5, 1, wr);
	vocoder ans = vocoder(1, 1, wr2);
	hi = v.getcurr(ans.iffttmp);
	//ans.img.writebmp("a.bmp");
	fir = hi;
	dec = 1;
	//return ;
	for (int i = 0; i < 13; i++) _aa[i] = _a[i];
	for (int i = 0; i < 12; i++) _dd[i] = _d[i];
	while (1) {
		if (!(rand() % 4)) {
			_a[j] += rand() % 3 - 1;
			if (_a[j] < -4) _a[j] = -4;
			else if (_a[j] > 4) _a[j] = 4;
			j++;
			if (j == 13) j = 0;
		}
		else 
		 {
			_d[k] += dec * ((double)(rand() % 32768) / 163840. - 0.1);
			if (k < 5) {
				if (_d[k] < -1) _d[k] = -1;
				else if (_d[k] > 1) _d[k] = 1;
			}
			else {
				if (_d[k] < 0) _d[k] = 0;
				else if (_d[k] > 1) _d[k] = 1;
			}
			k++;
			if (k == 12) k = 0;
		}

		vocoder t = vocoder(.5, 1, wr);
		double xn = t.getcurr(ans.iffttmp);
		if(xn < hi)
		{
			hi = xn;
			dec = hi / fir; dec *= dec;
			printf("%lf\n", hi);
			wavwriter("aa.wav", t);
			FILE* f = fopen("a.txt", "w");
			for (int i = 0; i < 13; i++) fprintf(f, "%i,", _a[i]);
			fprintf(f, "\n");
			for (int i = 0; i < 12; i++) fprintf(f, "%lf,", _d[i]);
			fclose(f);
			for (int i = 0; i < 13; i++) _aa[i] = _a[i];
			for (int i = 0; i < 12; i++) _dd[i] = _d[i];
		}
		else {
		printf("x ", hi, xn);
		for (int i = 0; i < 13; i++) _a[i] = _aa[i];
		for (int i = 0; i < 12; i++) _d[i] = _dd[i];
		}
	}
}
void main(int argc, char* argv[]) {
	
	wavwriter(argv[2], vocoder(atof(argv[3]), atof(argv[4]), wavreader(argv[1])));
	wavwriter(argv[2], vocoder_mtlb(atof(argv[3]), atof(argv[4]), wavreader(argv[1])));
	wavwriter(argv[2], vocoder_mtlb_sin3(atof(argv[3]), atof(argv[4]), wavreader(argv[1])));
	
}
