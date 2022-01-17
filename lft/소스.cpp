#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <windows.h>
#include <omp.h>
#include <iostream>
using namespace std;
_Dcomplex operator+(_Dcomplex a, _Dcomplex b) {
	_Dcomplex tmp;
	tmp._Val[0] = a._Val[0] + b._Val[0];
	tmp._Val[1] = a._Val[1] + b._Val[1];
	return tmp;
}
_Dcomplex operator-(_Dcomplex a, _Dcomplex b) {
	_Dcomplex tmp;
	tmp._Val[0] = a._Val[0] - b._Val[0];
	tmp._Val[1] = a._Val[1] - b._Val[1];
	return tmp;
}
_Dcomplex operator*(_Dcomplex a, _Dcomplex b) {
	_Dcomplex tmp;
	tmp._Val[0] = a._Val[0] * b._Val[0] - a._Val[1] * b._Val[1];
	tmp._Val[1] = a._Val[1] * b._Val[0] + a._Val[0] * b._Val[1];
	return tmp;
}
_Dcomplex operator*(_Dcomplex a, double b) {
	_Dcomplex tmp;
	tmp._Val[0] = a._Val[0] * b;
	tmp._Val[1] = a._Val[1] * b;
	return tmp;
}
_Dcomplex operator*(double b, _Dcomplex a) {
	_Dcomplex tmp;
	tmp._Val[0] = a._Val[0] * b;
	tmp._Val[1] = a._Val[1] * b;
	return tmp;
}
class fftmath {
protected:
	_Dcomplex mpii = { 0,-3.14159265358979323846 };
	
	int  maxn = 0;
	int overlap;
	_Dcomplex* cexpb = 0;
	fftmath(int overlap, int n) :overlap(overlap), maxn(n) {
		if (cexpb) free(cexpb);
		cexpb = (_Dcomplex*)malloc(sizeof(_Dcomplex) * maxn);
		for (int i = 0; i < maxn; i++) cexpb[i] = cexp(mpii * (static_cast<double>(i) / static_cast<double>(maxn)));
	}

	virtual ~fftmath() {
		if (cexpb) free(cexpb);
	}
public:
	void ifft() {
		mpii._Val[1] = 3.14159265358979323846;
		for (int i = 0; i < maxn; i++) cexpb[i] = cexp(mpii * (static_cast<double>(i) / static_cast<double>(maxn)));
	}
	void fft(){
		mpii._Val[1] = -3.14159265358979323846;
		for (int i = 0; i < maxn; i++) cexpb[i] = cexp(mpii * (static_cast<double>(i) / static_cast<double>(maxn)));
	}
public:
};
class fftcalc :public fftmath {
private:
	int n;
	int n2;
	_Dcomplex* sw, * sw2;
public:
	fftcalc(int overlap, int n) :n(n), fftmath(overlap, n) {
		if (overlap < 1) throw exception("overlap must be positive integer");
		if (n & (n - 1)) throw exception("n must be power of 2");
		n2 = n >> 1;
	}
	
	void operator()(_Dcomplex* in,_Dcomplex* out) {
		sw2 = out;
		for (int step = n / 2; step > 0; step /= 2) {
			for (int bias = 0; bias < step; bias++) {
				for (int i = 0; i < n; i += 2 * step) {
					_Dcomplex& ce1 = fftmath::cexpb[i];
					_Dcomplex& tmp2a = in[i + step + bias];
					_Dcomplex& tmp2b = in[i + bias];
					_Dcomplex& tmp3a = sw2[(i >> 1) + bias];
					_Dcomplex& tmp3b = sw2[((i + n) >> 1) + bias];
					
					tmp3a._Val[0] = tmp2b._Val[0] + (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
					tmp3a._Val[1] = tmp2b._Val[1] + (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
					tmp3b._Val[0] = tmp2b._Val[0] - (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
					tmp3b._Val[1] = tmp2b._Val[1] - (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
					//_Dcomplex t = cexpb[(maxn / n) * i] * in[i + step + bias];
					//sw2[(i >> 1) + bias] = in[i + bias] + t;
					//sw2[((i + n) >> 1) + bias] = in[i + bias] - t;
				}
			}
			sw = in;
			in = sw2;
			sw2 = sw;
		}
		for (int i = 0; i < n; i++) out[i] = in[i];
	}
	int getoutsize() {
		return n2;
	}
};

class doublearray {
public:
	double* arr = 0;
	int n;
	void create(int n) {
		this->n = n;
		arr = (double*)malloc(sizeof(double) * n);
	}
	void copy(doublearray da) {
		n = da.n;
		arr = (double*)malloc(sizeof(double) * n);
		for (int i = 0; i < n; i++) arr[i] = da.arr[i];
	}
	void remove() {
		if (arr) {
			free(arr);
			arr = 0;
		}
	}
};
class wavreader {
private:
	void* data;
	int smpl, nbit, count;
	_Dcomplex* dout = 0;
	int ndout = 0;
	double noiselevel = 0;
public:
	wavreader(const char* filename) {
		int hdr[11];
		FILE* f = fopen(filename, "rb");
		if (f != 0) {
			fread(hdr, 4, 11, f);
			smpl = hdr[6];
			nbit = hdr[8] >> 19;
			count = hdr[10] / nbit;
			data = reinterpret_cast<void*>(malloc(hdr[10]));
			fread(data, nbit, count, f);
			fclose(f);
		}
		else throw exception("file not exist");

	}
	_Dcomplex* operator()(int offset, int n) {
		if (n > ndout) {
			ndout = n;
			dout = (_Dcomplex*)realloc((void*)dout, n * sizeof(_Dcomplex));
		}
		_Dcomplex* ptr1 = dout, * ptr2 = dout + n;
		switch (nbit) {
		case 1:
		{
			char* tmp0 = (char*)data;
			char* tmp = tmp0 + offset;
			char* tmp2 = tmp0 + count;
			while (ptr1 != ptr2) {
				if (tmp < tmp2 && tmp >= tmp0)
					ptr1->_Val[0] = ((double)(*(tmp++) + 128)) / 128.;
				else ptr1->_Val[0] = 0., tmp++;
				ptr1->_Val[1] = 0.;
				ptr1++;
			}
		}
		break;
		case 2:
		{
			short* tmp0 = (short*)data;
			short* tmp = tmp0 + offset;
			short* tmp2 = tmp0 + count;
			while (ptr1 != ptr2) {
				if (tmp < tmp2 && tmp >= tmp0)
					ptr1->_Val[0] = ((double)(*(tmp++))) / 32768.;
				else ptr1->_Val[0] = 0., tmp++;
				ptr1->_Val[1] = 0.;
				ptr1++;
			}
		}
		break;
		case 4:
		{
			long* tmp0 = (long*)data;
			long* tmp = tmp0 + offset;
			long* tmp2 = tmp0 + count;
			while (ptr1 != ptr2) {
				if (tmp < tmp2 && tmp >= tmp0) ptr1->_Val[0] = ((double)(*(tmp++))) / 2147483648.;
				else ptr1->_Val[0] = 0., tmp++;
				ptr1->_Val[1] = 0.;
				ptr1++;
			}
		}
		break;
		default:
			throw exception("invaild nbit");
			break;
		}
		return dout;
	}
	int getsize() {
		return count;
	}
	int getsmplrate() {
		return smpl;
	}
	void* getptr() {
		return data;
	}
	int getnbit() {
		return nbit;
	}
};
class wavwriter {
	unsigned long ret[11];
	const char* filename;
public:
	wavwriter(int samprate, const char* filename) :filename(filename) {
		
		
		ret[0] = 0x46464952;
		ret[2] = 0x45564157;
		ret[3] = 0x20746D66;
		ret[4] = 0x00000010;
		ret[5] = 0x00010001;
		ret[6] = samprate;
		ret[7] = ret[6] * 2;
		ret[8] = 0x00100002;
		ret[9] = 0x61746164;
	}
	void operator()(doublearray dat) {
		short* temp = (short*)malloc(sizeof(short) * dat.n);
		ret[10] = dat.n * 2;
		ret[1] = dat.n * 2 + 0x24;
		FILE* f = fopen(filename, "wb");
		for (int i = 0; i < dat.n; i++) {
			if (dat.arr[i] >= 32768.) temp[i] = 32767;
			else if (dat.arr[i] < -32768.) temp[i] = -32768.;
			else temp[i] = (short)(dat.arr[i] * 32768.);
		}
		fwrite(ret, 1, 44, f);
		fwrite(temp, 2, dat.n, f);
		fclose(f);
		free(temp);
	}
	~wavwriter() {
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
			fwrite(b2, 1, 28, f);
		}
		else throw exception("file write failed");
	}
	void operator()(double* data, int yindex) {
		for (int j = 0; j < w; j++) {
			int t = static_cast<int>(-20. * log10(1e-20 + data[j])) + 50;
			if (t < 0) t = 0;
			if (t > 81) t = 81;
			out[(j + yindex * w) * 3 + 2] = red[t]; out[(j + yindex * w) * 3 + 1] = grn[t]; out[(j + yindex * w) * 3] = blu[t];
		}
	}
	void operator()(_Dcomplex* data, int yindex) {
		for (int j = 0; j < w; j++) {
			int t = static_cast<int>(-20. * log10(1e-20 + cabs(data[j]))) + 50;
			if (t < 0) t = 0;
			if (t > 81) t = 81;
			out[(j + yindex * w) * 3 + 2] = red[t]; out[(j + yindex * w) * 3 + 1] = grn[t]; out[(j + yindex * w) * 3] = blu[t];
		}
	}
	void p(_Dcomplex* data, int yindex) {
		for (int j = 0; j < w; j++) {
			double bri = static_cast<int>(-20. * log10(1e-20 + cabs(data[j]))) + 50;
			if (bri < 0) bri = 0;
			if (bri > 81) bri = 81;
			bri = 81 - bri;
			double hue = carg(data[j])/ 6.283185307179586476925286766559+0.5;
			double r, g, b;
			if (hue < 0.3333333333333) {
				r = 1.-hue * 3.;
				g = hue * 3.;
				b = 0;
			}
			else if (hue < 0.6666666666666666) {
				b = hue * 3 - 1.;
				g = 2 - hue * 3;
				r = 0;
			}
			else {
				r = hue * 3 - 2;
				g = 0;
				b = 3 - hue * 3;
			}
			out[(j + yindex * w) * 3 + 2] = (unsigned char)(r*(double)bri*1.5); 
			out[(j + yindex * w) * 3 + 1] = (unsigned char)(g * (double)bri * 1.5);
			out[(j + yindex * w) * 3] = (unsigned char)(b * (double)bri * 1.5);
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
class image {
public:
	int x, y;
	double* temp;
	_Dcomplex** data;
	void createimage(int x, int y) {
		this->x = x;
		this->y = y;
		temp = (double*)malloc(sizeof(double) * x);
		data = (_Dcomplex**)malloc(sizeof(_Dcomplex*) * (y));
		for (int i = 0; i < y; i++) 
			data[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * x);
	}
	void removeimage() {
		if (data != 0) {
			for (int i = 0; i < y; i++)free(data[i]);
			free(data);
			data = 0;
			free(temp);
		}
	}
	void writebmp(const char* filename) {
		bmpwriter bm(filename, x, y);
		for (int i = 0; i < y; i++) {
			bm.p(data[i], i);
		}
	}
};

class driver {
public:

	int n, overlap;
	int shift;
	int shift2;
	double ratio;
	_Dcomplex* dctmp;
	fftcalc* ft;
	fftcalc* ft32[32];
	_Dcomplex* dctmp32[32];
	int x;
	double vol;
	int smpl = 44100;
	_Dcomplex im = { 0,1 };
	double* hann;
	double sinetable[65536];

	driver(int n, int overlap,double ratio) :n(n), ratio(ratio), overlap(overlap) {
		ft = new fftcalc(overlap, n);
		ft->fft();
		for (int i = 0; i < 32; i++) {
			ft32[i] = new fftcalc(overlap, n);
			ft32[i]->fft();
		}
		for (int i = 0; i < 32; i++) ft32[i] = new fftcalc(overlap, n);
		shift2 = n / overlap;
		shift = (int)((double)shift2/ratio);
		x = n;
		vol = 2. / (double)(n*overlap);
		dctmp = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		hann = (double*)malloc(sizeof(double) * n);
		for (int i = 0; i < n; i++) {
			hann[i] = sin((double)i / n * 3.1415926535898);
			hann[i] *= hann[i]; hann[i] *= hann[i];
		}
		for (int i = 0; i < 65536; i++) sinetable[i] = sin(6.2831853071795864 / 65536. * (double)i);

	}

	image wfile2img(char* filename) {
		wavreader wr(filename);
		smpl = wr.getsmplrate();
		int datasize = wr.getsize();
		int nbit = wr.getnbit();
		int y = (datasize + shift - 1) / shift;
		_Dcomplex* dout[32];
		for(int i=0;i<32;i++)
			dout[i]= (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		image fft1;
		fft1.createimage(x, y);
		void* data = wr.getptr();
		switch (nbit) {
		case 1:
#pragma omp parallel for num_threads(32)
			for (int i = 0; i < y; i++) {
				int numthread = omp_get_thread_num();
				int offset = i * shift - n;
				_Dcomplex* ptr1 = dout[numthread], * ptr2 = dout[numthread] + n;
				char* tmp0 = (char*)data;
				char* tmp = tmp0 + offset;
				char* tmp2 = tmp0 + datasize;
				while (ptr1 != ptr2) {
					if (tmp < tmp2 && tmp >= tmp0)
						ptr1->_Val[0] = ((double)(*(tmp++) + 128)) / 128.;
					else ptr1->_Val[0] = 0., tmp++;
					ptr1->_Val[1] = 0.;
					ptr1++;
				}
				for (int i = 0; i < n; i++) dout[numthread][i]._Val[0] *= hann[i];
				(*ft32[numthread])(dout[numthread], fft1.data[i]);
			}
			break;
		case 2:
#pragma omp parallel for num_threads(32)
			for (int i = 0; i < y; i++) {
				int numthread = omp_get_thread_num();
				int offset = i * shift - n;
				_Dcomplex* ptr1 = dout[numthread], * ptr2 = dout[numthread] + n;
				short* tmp0 = (short*)data;
				short* tmp = tmp0 + offset;
				short* tmp2 = tmp0 + datasize;
				while (ptr1 != ptr2) {
					if (tmp < tmp2 && tmp >= tmp0)
						ptr1->_Val[0] = ((double)(*(tmp++))) / 32768.;
					else ptr1->_Val[0] = 0., tmp++;
					ptr1->_Val[1] = 0.;
					ptr1++;
				}
				for (int i = 0; i < n; i++) dout[numthread][i]._Val[0] *= hann[i];
				(*ft32[numthread])(dout[numthread], fft1.data[i]);
			}
			break;
		
		case 4:
#pragma omp parallel for num_threads(32)
			for (int i = 0; i < y; i++) {
				int numthread = omp_get_thread_num();
				int offset = i * shift - n;
				_Dcomplex* ptr1 = dout[numthread], * ptr2 = dout[numthread] + n;
				long* tmp0 = (long*)data;
				long* tmp = tmp0 + offset;
				long* tmp2 = tmp0 + datasize;
				while (ptr1 != ptr2) {
					if (tmp < tmp2 && tmp >= tmp0) ptr1->_Val[0] = ((double)(*(tmp++))) / 2147483648.;
					else ptr1->_Val[0] = 0., tmp++;
					ptr1->_Val[1] = 0.;
					ptr1++;
				}
				for (int i = 0; i < n; i++) dout[numthread][i]._Val[0] *= hann[i];
				(*ft32[numthread])(dout[numthread], fft1.data[i]);
			}
			break;
		
		default:
			throw exception("invaild nbit");
			break;
		}
		for (int i = 0; i < 32; i++) free(dout[i]);
		return fft1;
	}
	image wav2img(doublearray data) {
		int y = (data.n + shift - 1) / shift;
		image fft1;
		fft1.createimage(x, y);
		int ii;
		for (int i = 0; i < y; i++) {
			for (int j = 0; j < n; j++) {
				ii = i * shift + j - n;
				if (ii >= data.n || ii < 0) dctmp[j]._Val[0] = 0;
				else dctmp[j]._Val[0] = data.arr[ii];
				dctmp[j]._Val[1] = 0;
			}
			for (int i = 0; i < n; i++) dctmp[i]._Val[0] *= hann[i];
			(*ft)(dctmp, fft1.data[i]);
		}
		return fft1;
	}

	void modifyphase(image img) {
		double v1=1, v2=0;
		if (ratio > 1) {
			 v1 = 0.95 - cos(ratio * 3.1415926535897932384626433832795) * 0.05;
			v2 = 0.05 + cos(ratio * 3.1415926535897932384626433832795) * 0.05;
		}
#pragma omp parallel for
		for (int b = 0; b <= x / 2; b++) {
			double op = carg(img.data[1][b]);
			for (int a = 2; a < img.y; a++) {
				double base = carg(img.data[a - 1][b]);
				double arg1 = carg(img.data[a][b]);
				double darg = arg1 - op;
				double tmp2 = (double)shift2/ratio * 6.283185307179586476925286766559 * ((double)b) / (double)x;
				double darg2 = (double)floor((tmp2 - darg)/ 6.283185307179586476925286766559+0.5)* 6.283185307179586476925286766559;
				darg += darg2;
				double tmp = ratio * (darg*v1+tmp2*v2);
				op = carg(img.data[a][b]);
				img.data[a][b] = cexp(im * (base +tmp)) * cabs(img.data[a][b]);
			}
		}
		for (int b = x / 2 + 1; b < x; b++)
			for (int a = 1; a < img.y; a++) {
				img.data[a][b]._Val[0] = img.data[a][x - b]._Val[0];
				img.data[a][b]._Val[1] = -img.data[a][x - b]._Val[1];
			}
	}
	void overlapf(double* targ, int shif,int y, _Dcomplex *src) {
		for (int i = 0; i < n; i++) src[i]._Val[0] *= vol;
		if (shif == 0) {
			for (int j = 0; j < n; j++)
				targ[j] += hann[j] *src[j]._Val[0];
			return;
		}
		double pearson_max = -1;
		int imax=-1;
#pragma omp parallel for
		for (int i = -(shift2>>3); i < (shift2>>3); i++) {
			double var1 = 0, cov = 0, var2 = 0;
			for (int j = n>>2; j < n>>1; j++) {
				double t1, t2;
				if (j + i + shif < y)
					t1 = targ[j + i + shif];
				else t1 = 0;
				t2 = src[j]._Val[0];
				var1 += t1 * t1;
				var2 += t2 * t2;
				cov += t1 * t2;
			}
			double pearson = cov * cov / (1e-12 + var1 * var2);
			if (i > 0) pearson -= (double)i / 6;
			else pearson += (double)i / 6;
			if (pearson > pearson_max) {
#pragma omp critical
				{
					pearson_max = pearson;
					imax = i;
				}
			}
		}
		
		for (int j = 0; j < n; j++)
			if(imax+j+shif<y)
			targ[imax + j+shif] += hann[j]*src[j]._Val[0];
	}
	doublearray img2wav(image img) {
		for(int i=0;i<32;i++) ft32[i]->ifft();
		_Dcomplex* dout[32];
		for (int i = 0; i < 32; i++)
			dout[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		int len = img.y * shift2+img.x;
		doublearray da;
		da.create(len);
		for (int i = 0; i < len; i++) {
			da.arr[i] = 0;
		}
		for (int i = 0; i < img.y; i++) {
			(*ft32[0])(img.data[i], dout[0]);
			overlapf(da.arr,shift2 * i,len, dout[0]);
		}
		return da;
	}
	void img2wfile(image img, const char* filename) {
		wavwriter wwr(smpl,filename);
		wwr(img2wav(img));
	}
	
	virtual ~driver() {
		free(dctmp);
	}
};



void main(int argc, char* argv[]) {
	double ratio = atof(argv[3]);
	int n = 4096;
	int overlap;
	if (ratio > 1.2) overlap = 4;
	else if (ratio > 0.9) overlap = 8;
	else overlap = 16;
	driver drv(n, overlap,ratio);
	image img1 = drv.wfile2img(argv[1]);
	drv.modifyphase(img1);
	drv.img2wfile(img1, argv[2]);
}
