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
class fftmath {
protected:
	const _Dcomplex mpii = { 0,-3.14159265358979323846 };
	int outsize = 0, n = 0;
	int overlap;
	_Dcomplex* out = 0, * cexpb = 0;
	fftmath(int overlap, int n) :overlap(overlap), n(n) {
		if (cexpb) free(cexpb);
		cexpb = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		for (int i = 0; i < n; i++) cexpb[i] = cexp(mpii * (static_cast<double>(i) / static_cast<double>(n)));
	}
	int log2n(int n) {
		int ret = -1;
		while (n) ret++, n >>= 1;
		return ret;
	}
	virtual ~fftmath() {
		if (cexpb) free(cexpb);
	}
public:
	int getoutsize() {
		return outsize;
	}
	int getinsize() {
		return n;
	}
};
class lftcalc :virtual public fftmath {
private:
	int* ovldata = 0;
	int** LeftEndIndex;
	int** RightEndIndex;
	int* StageSize;
	int nnly;
	int stepsize;
	int nnly2;

	_Dcomplex* tmp1, * tmp2, * tmp3, * sw;
	_Dcomplex* mem1, * mem2;
	int  r, s, q;
public:
	double* centerfreq;
	double* vol;
	lftcalc(int overlap, int n, int r, int s, int q) : r(r), s(s), q(q), fftmath(overlap, n) {
		if (overlap < 1) throw exception("overlap must be positive integer");
		if (n & (n - 1)) throw exception("n must be power of 2");
		if (n < 2) throw exception("n must be larger than 1");
		if (r & (r - 1)) throw exception("r must be power of 2");
		if (r < 2) throw exception("r must be larger than 1");
		if (s & (s - 1)) throw exception("s must be power of 2");
		if (s < 1) throw exception("s must be larger than 0");
		if (s >= r) throw exception("s must be smaller than r/2");
		if (q < 2) throw exception("q must be larger than 1");
		if (q & (q - 1)) throw exception("q must be power of 2");
		nnly = log2n(n);
		stepsize = log2n(r / s);
		LeftEndIndex = (int**)malloc(sizeof(int*) * nnly);
		RightEndIndex = (int**)malloc(sizeof(int*) * nnly);
		StageSize = (int*)malloc(sizeof(int) * nnly);

		mem1 = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		mem2 = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		int t = 0;
		int nn = nnly / stepsize;

		for (int i = 0; i < nn; i++)
			StageSize[i] = nnly - i * stepsize;
		for (int i = 0; i < nnly; i++) {
			LeftEndIndex[i] = (int*)malloc(sizeof(int) * StageSize[i]);
			RightEndIndex[i] = (int*)malloc(sizeof(int) * StageSize[i]);
		}
		LeftEndIndex[0][StageSize[0] - 1] = 0;
		RightEndIndex[0][StageSize[0] - 1] = r;
		outsize = r;
		int p;
		t++;
		int ss = s;
		int rr = r;
		while (1) {
			rr *= q;
			LeftEndIndex[t][StageSize[t] - 1] = ss;
			p = rr << (t * stepsize);
			if (p >= (n >> 1)) {
				RightEndIndex[t][StageSize[t] - 1] = rr * (n >> 1) / p;
				outsize += RightEndIndex[t][StageSize[t] - 1] - ss;
				nnly2 = t;
				break;
			}
			else RightEndIndex[t][StageSize[t] - 1] = rr;
			t++;
			outsize += rr - ss;
			ss *= q;
		}
		for (int i = 0; i <= nnly2; i++) {
			for (int j = StageSize[i] - 2; j >= 0; j--) {
				if (RightEndIndex[i][j + 1] > (n >> (i * stepsize + 1))) {
					LeftEndIndex[i][j] = 0;
					RightEndIndex[i][j] = n >> (i * stepsize);
				}
				else {
					LeftEndIndex[i][j] = LeftEndIndex[i][j + 1] << 1;
					RightEndIndex[i][j] = RightEndIndex[i][j + 1] << 1;
				}
			}
		}
		out = (_Dcomplex*)malloc(sizeof(_Dcomplex) * outsize);
		centerfreq = (double*)malloc(sizeof(double) * outsize);
		vol = (double*)malloc(sizeof(double) * outsize);
		int j = 0, k = 1, l = 0, m = 0;
		for (int i = 0; i < outsize; i++) {
			centerfreq[i] = 65536. * ((double)j) / (double)n;
			vol[i] = (double)k / (double)n;
			j += k; l++;
			if (l >= RightEndIndex[m][StageSize[m] - 1]) {
				m++;
				l = LeftEndIndex[m][StageSize[m] - 1];
				k *= r / s;
			}
		}
	}
	_Dcomplex* ft(_Dcomplex* in) {
		int step, tstep;
		int _n;
		int cnt = 0;
		_Dcomplex* inptr = 0;
		for (int p = 0; p <= nnly2; p++) {
			_n = n >> (p * stepsize);
			inptr = in + (n - _n);
			tmp1 = mem1;
			tmp2 = inptr;
			tmp3 = mem2;

			for (int step2 = 0; step2 < StageSize[p]; step2++) {
				step = _n >> (step2 + 1);
				tstep = step << 1;
				if (RightEndIndex[p][step2] == _n) {
					int i = 0;
					int end = LeftEndIndex[p][step2] << 1;
					for (; i < end; i += tstep) {
						for (int offset = 0; offset < step; offset++) {

							_Dcomplex& ce1 = fftmath::cexpb[i << (p * stepsize)];
							_Dcomplex& tmp2a = tmp2[i + step + offset];
							_Dcomplex& tmp2b = tmp2[i + offset];
							_Dcomplex& tmp3b = tmp1[((i + _n) >> 1) + offset];
							tmp3b._Val[0] = tmp2b._Val[0] - (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
							tmp3b._Val[1] = tmp2b._Val[1] - (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
							//_Dcomplex t = fftmath::cexpb[i << (p*stepsize)] * tmp2[i + step + offset];
							//tmp1[((i + _n) >> 1) + offset] = tmp2[i + offset] - t;
						}
					}
					for (; i < _n; i += tstep) {
						for (int offset = 0; offset < step; offset++) {
							_Dcomplex& ce1 = fftmath::cexpb[i << (p * stepsize)];
							_Dcomplex& tmp2a = tmp2[i + step + offset];
							_Dcomplex& tmp2b = tmp2[i + offset];
							_Dcomplex& tmp3a = tmp1[(i >> 1) + offset];
							_Dcomplex& tmp3b = tmp1[((i + _n) >> 1) + offset];
							tmp3a._Val[0] = tmp2b._Val[0] + (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
							tmp3a._Val[1] = tmp2b._Val[1] + (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
							tmp3b._Val[0] = tmp2b._Val[0] - (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
							tmp3b._Val[1] = tmp2b._Val[1] - (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
							//_Dcomplex t = fftmath::cexpb[i << (p*stepsize)] * tmp2[i + step + offset];
							//tmp1[(i >> 1) + offset] = tmp2[i + offset] + t;
							//tmp1[((i + _n) >> 1) + offset] = tmp2[i + offset] - t;
						}
					}
				}
				else {
					int end = RightEndIndex[p][step2] << 1;
					for (int i = (LeftEndIndex[p][step2] << 1); i < end; i += tstep) {
						for (int offset = 0; offset < step; offset++) {
							_Dcomplex& ce1 = fftmath::cexpb[i << (p * stepsize)];
							_Dcomplex& tmp2a = tmp2[i + step + offset];
							_Dcomplex& tmp2b = tmp2[i + offset];
							_Dcomplex& tmp3a = tmp1[(i >> 1) + offset];
							tmp3a._Val[0] = tmp2b._Val[0] + (ce1._Val[0] * tmp2a._Val[0] - ce1._Val[1] * tmp2a._Val[1]);
							tmp3a._Val[1] = tmp2b._Val[1] + (ce1._Val[1] * tmp2a._Val[0] + ce1._Val[0] * tmp2a._Val[1]);
							//_Dcomplex t = fftmath::cexpb[i << (p*stepsize)] * tmp2[i + step + offset];
							//tmp1[(i >> 1) + offset] = tmp2[i + offset] + t;

						}
					}
				}
				if (tmp2 == inptr) {
					tmp2 = tmp1;
					tmp1 = tmp3;
				}
				else {
					sw = tmp1;
					tmp1 = tmp2;
					tmp2 = sw;
				}
			}
			int end = RightEndIndex[p][StageSize[p] - 1];
			for (int i = LeftEndIndex[p][StageSize[p] - 1]; i < end; i++) out[cnt++] = tmp2[i];
		}
		return out;
	}
	double* operator()(_Dcomplex* in, double* dout) {
		out = ft(in);
		for (int b = 0; b < outsize; b++)
			dout[b] = cabs(out[b]);
		return dout;
	}

	virtual ~lftcalc() {
		for (int i = 0; i < nnly; i++) {
			free(LeftEndIndex[i]);
			free(RightEndIndex[i]);
		}
		free(LeftEndIndex);
		free(RightEndIndex);
		free(StageSize);
		free(mem1);
		free(mem2);
		free(out);
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
};
class wavwriter {
	unsigned long* ret;
	const char* filename;
public:
	wavwriter(int samprate, const char* filename) :filename(filename) {
		ret = reinterpret_cast<unsigned long*>(malloc(44));
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
		free(ret);
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
	double** data;

	void createimage(int x, int y) {
		this->x = x;
		this->y = y;
		data = (double**)malloc(sizeof(double*) * (y + 1));
		for (int i = 0; i < y + 1; i++) data[i] = (double*)malloc(sizeof(double) * x);
	}
	void removeimage() {
		if (data != 0) {
			for (int i = 0; i < y; i++)free(data[i]);
			free(data);
			data = 0;
		}
	}
	void writebmp(char* filename) {
		bmpwriter bm(filename, x, y);
		for (int i = 0; i < y; i++) {
			bm(data[i], i);
		}
	}
};

class driver {
public:
	static int srandv;
	int n, r, overlap;
	int s, t;
	int shift;
	_Dcomplex* dctmp;
	lftcalc* ft;
	lftcalc* ft32[32];
	_Dcomplex* dctmp32[32];
	int x;
	double sinetable[65536];
	doublearray syncp32[32];
	int syncplen;
	image syncpi32[32];
	double** centerfreql;
	double* vol;
	driver(int n, int r, int overlap) :n(n), r(r), overlap(overlap) {
		s = r / 2;
		t = 2;
		shift = n / overlap;
		ft = new lftcalc(overlap, n, r, s, t);
		for (int i = 0; i < 32; i++) ft32[i] = new lftcalc(overlap, n, r, s, t);
		syncplen = 2 * n;// -shift;
		for (int i = 0; i < 32; i++) syncp32[i].create(syncplen);
		x = ft->getoutsize();
		vol = ft->vol;
		centerfreql = (double**)malloc(sizeof(double*) * 2);
		centerfreql[0] = ft->centerfreq;
		centerfreql[1] = (double*)malloc(sizeof(double) * x);
		for (int i = 0; i < x - 1; i++) {
			centerfreql[1][i] = .5 * (centerfreql[0][i] + centerfreql[0][i + 1]);
		}
		centerfreql[1][x - 1] = 16384. + .5 * centerfreql[0][x - 1];

		dctmp = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
		for (int i = 0; i < 32; i++) {
			dctmp32[i] = (_Dcomplex*)malloc(sizeof(_Dcomplex) * n);
			syncpi32[i].createimage(x, overlap);
		}
		for (int i = 0; i < 65536; i++) sinetable[i] = sin(6.2831853071795864 / 65536. * (double)i);
	}
	inline double fastsin(unsigned short tpiv) const {
		return sinetable[tpiv];
	}
	double fastrand() {
		srandv = srandv * 1103515245 + 12345;
		double temp = (double)((srandv >> 16) & 65535);
		srandv = srandv * 1103515245 + 12345;
		return temp + (double)((srandv >> 16) & 32767) / 32768.;
	}
	image wfile2img(char* filename) {
		wavreader wr(filename);
		int datasize = wr.getsize();
		int y = (datasize + shift - 1) / shift;
		_Dcomplex* tmp;
		image fft1;
		fft1.createimage(x, y);
		for (int i = 0; i < y; i++) {
			tmp = wr(i * shift - n, n);
			(*ft)(tmp, fft1.data[i]);
		}
		for (int i = 0; i < x; i++)
			fft1.data[y][i] = 0;
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
			(*ft)(dctmp, fft1.data[i]);
		}
		for (int i = 0; i < x; i++)
			fft1.data[y][i] = 0;
		return fft1;
	}
	image stratchimage(image img, double ratio) {
		int yy = (int)((double)img.y / ratio);
		image fft2;
		fft2.createimage(x, yy);
		double idx = 0;
		for (int i = 0; i < yy; i++, idx += ratio) {
			int lo = (int)idx;
			int hi = lo + 1;
			double rl = idx - (double)lo;
			double rh = (double)hi - idx;
			for (int j = 0; j < x; j++)
				fft2.data[i][j] = img.data[lo][j] * rl + img.data[hi][j] * rh;
		}
		for (int i = 0; i < x; i++)
			fft2.data[yy][i] = 0;
		return fft2;
	}
	image diffimage(image a, image b) {
		image c;
		c.createimage(x, a.y);
		for (int i = 0; i < a.y; i++)
			for (int j = 0; j < x; j++) {
				if (a.data[i][j] < b.data[i][j]) c.data[i][j] = b.data[i][j] - a.data[i][j];
				else c.data[i][j] = a.data[i][j] - b.data[i][j];
			}
		return c;
	}
	float fastsqrt(float z) {
		union { float f; uint32_t i; } val = { z };
		val.i -= 1 << 23;
		val.i >>= 1;
		val.i += 1 << 29;
		return val.f;
	}
	double getmag(image a) {
		double v = 0;
		for (int i = 0; i < a.y; i++)
			for (int j = 0; j < x; j++)
				v += fastsqrt(a.data[i][j]);
		return v;
	}
	doublearray sync(image a) {
		doublearray da;
		int arrsize = shift * a.y;
		da.create(arrsize);
		double* phase = (double*)malloc(sizeof(double) * x);

		for (int i = 0; i < x; i++) phase[i] = fastrand();
		for (int i = 0; i < da.n; i++) da.arr[i] = 0;
#pragma omp parallel for
		for (int k = 0; k < x; k++) {
			printf("%d/%d\n", k, x);
			for (int i = 0; i < a.y; i++) {
				int idx = i * shift;
				for (int j = 0; j < shift; j++) {
					phase[k] += centerfreql[0][k];
					while (phase[k] > 262144.) phase[k] -= 262144.;
					double tmp = 1.4 * vol[k] * a.data[i][k] * fastsin((unsigned short)(phase[k]));
#pragma omp atomic
					da.arr[idx] += tmp;
					idx++;
				}

			}
		}
		return da;
	}
	void sync_part(image& targ, doublearray& srcd, image& srci, image& diff) {
		int y = targ.y;
		//for (int xsh = 0; xsh < 2; xsh ++) 
		int xsh = 0;
		int ysh = 0;
		{
#pragma omp parallel for num_threads(32)
			for (int ys = ysh; ys < y; ys += 2 * overlap) {
				printf("!%d/%d\n", ys, y);
				int threadid =  omp_get_thread_num();
				for (int xs = 0; xs < x; xs += 20) {
					//printf(".%d/%d\n", xs, x);
					for (int xc = xs; (xc < xs + 20) && (xc < x); xc++) {
							for (int yc = ys; (yc < ys + 2 * overlap) && (yc < y); yc+=2) {
								int ycl = yc + 1;
								int ych = yc + overlap+1;
								if (ych >= y) ych = y;
								int ycm = (ycl + ych) >> 1;
								
								int p = (yc + 1) * shift - n;
								double phase = fastrand();
								for (int q = 0; q < syncplen; q++) {
									if (p + q < 0 || p + q >= srcd.n) syncp32[threadid].arr[q] = 0;
									else syncp32[threadid].arr[q] = srcd.arr[p + q];
								}
								for (int q = n - shift; q < n+shift; q++) {
									phase += centerfreql[xsh][xc];
									while (phase > 262144.) phase -= 262144.;
									syncp32[threadid].arr[q] += 0.3 * vol[xc] * diff.data[ycm][xc] * fastsin((unsigned short)(phase));
								}
								for (int i = 0; i < overlap+1; i++) {
									for (int j = 0; j < n; j++) {//x
										int ii = i * shift + j;
										dctmp32[threadid][j]._Val[0] = syncp32[threadid].arr[ii];
										dctmp32[threadid][j]._Val[1] = 0;
									}
									(*ft32[threadid])(dctmp32[threadid], syncpi32[threadid].data[i]);
								}
								double orig = 0; double ddiff;
								for (int i = 0; i < overlap+1; i++) {
									if (ycl + i < y) {
										for (int j = 0; j < x; j++) {
											ddiff = srci.data[ycl + i][j] - targ.data[ycl + i][j];
											if (ddiff < 0) orig += fastsqrt(-ddiff);
											else orig += fastsqrt(ddiff);
											ddiff = syncpi32[threadid].data[i][j] - targ.data[ycl + i][j];
											if (ddiff < 0) orig -= fastsqrt(-ddiff);
											else orig -= fastsqrt(ddiff);
										}
									}
								}
								if (orig>0) {//¼öÁ¤
									p = yc * shift;
									for (int i = 0; i < shift*2; i++)
										srcd.arr[p + i] = syncp32[threadid].arr[n - shift + i];
									for (int i = 0; i < overlap+1; i++) {
										if (ycl + i < y)
											for (int j = 0; j < x; j++) 
												srci.data[ycl + i][j] = syncpi32[threadid].data[i][j];
									}
								}

							}
						
					}
				}
			}
		}
	}
	virtual ~driver() {
		for(int i=0;i<32;i++)syncp32[i].remove();
	}
};

int driver::srandv = 0;

void main(int argc, char* argv[]) {
	driver::srandv = time(0);
	double ratio = .8;
	int n = 16384;
	int r = 32;
	int overlap = 64;

	driver drv(n, r, overlap);
	image img1 = drv.wfile2img((char*)"a.wav");
	image img2 = drv.stratchimage(img1, ratio);
	doublearray da = drv.sync(img2);
	image img3 = drv.wav2img(da);
	image img4 = drv.diffimage(img2, img3);
	wavwriter ww(44100, "b.wav");
	ww(da);
	img2.writebmp((char*)"a.bmp");
	img3.writebmp((char*)"b.bmp");
	img4.writebmp((char*)"c.bmp");
	wavwriter ww2(44100, "b1.wav");
	for (int i = 0; i < 20; i++) {
		cout << drv.getmag(img2) << " " << drv.getmag(img3) << " " << drv.getmag(img4) << endl;
		drv.sync_part(img2, da, img3, img4);
		img4 = drv.diffimage(img2, img3);
		ww2(da);
		
	}
	

	cout << drv.getmag(img2) << " " << drv.getmag(img3) << " " << drv.getmag(img4) << endl;
}
