#include <iostream>
#include <iomanip>
#include <windows.h>
#include <gdiplus.h>
#include <fstream>
#include <functional>
using namespace Gdiplus;
using namespace std;
#pragma comment (lib, "Gdiplus.lib")
void Log(string tag, string s)
{
#ifdef _DEBUG
	cout << "[" << tag << "]" << s << endl;
#endif
}
template <class T>
void fileLog(string filename, T& obj)
{
#ifdef _DEBUG
	ofstream file;
	file.open(filename);
	file << obj;
	file.close;
#endif
}
template <class T>
class MATRIX
{
private:
	int W, H;
	T** A;
public:
	MATRIX() { A = 0; W = H = 0; }
	MATRIX(int h, int w) { init(h, w); }
	MATRIX(MATRIX& other) { init(other.H, other.W); *this = other; }
	MATRIX(MATRIX&& other) { init(other.H, other.W); *this = other; }
	~MATRIX() { clear(); }
	void init(int h, int w)
	{
		W = w;
		H = h;
		A = new T * [H];
		for (int i = 0; i < H; i++)
		{
			A[i] = new T[W];
			for (int j = 0; j < W; j++)
				A[i][j] = 0;
		}
	}
	void clear()
	{
		if (A)
		{
			for (int i = 0; i < H; i++)
				delete[] A[i];
			delete[] A;
		}
		A = 0;
	}
	void set(int i, int j, const T& v)
	{
		A[i][j] = v;
	}
	template <int h, int w>
	void set(T(&m)[h][w]) {
		init(h, w);
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] = *(*(m + j) + i);
	}
	void set(MATRIX& other) {
		clear();
		init(other.H, other.W);
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] = other.A[i][j];
	}
	int getH() { return H; }
	int getW() { return W; }
	int countValue(T v) {
		int n = 0;
		for (int y = 0; y < H; y++)
			for (int x = 0; x < W; x++)
				n += (A[x][y] == v);
		return n;
	}
	void bound(T min, T max) {
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
			{
				if (A[i][j] >= max) A[i][j] = max;
				else if (A[i][j] <= min) A[i][j] = min;
			}
	}
	void interbound(T min, T max) {
		int i, j;
		T a = 0;
		T b = 0;
		a = b = A[0][0];
		for (i = 0; i < H; i++)
			for (j = 0; j < W; j++)
			{
				if (A[i][j] >= a) a = A[i][j];
				else if (A[i][j] <= b) b = A[i][j];
			}
		for (i = 0; i < H; i++)
			for (j = 0; j < W; j++)
				A[i][j] = min + (max - min) * (A[i][j] - b) / (a - b);
	}
	void transpose()
	{
		MATRIX t(W, H);
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				t.A[j][i] = A[i][j];
		clear();
		init(t.H, t.W);
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[j][i] = t.A[i][j];
	}
	T dot(MATRIX& m)
	{
		T d = 0;
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				d += m.A[i][j] * A[i][j];
		return d;
	}
	void setAll(T c)
	{
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] = c;
	}
	T& operator()(int i, int j) { return A[i - 1][j - 1]; }
	T* operator[](int i) { return A[i]; }
	MATRIX& operator=(const MATRIX& m)
	{
		clear();
		init(m.H, m.W);
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] = m.A[i][j];
		return *this;
	}
	void operator*=(const T& s)
	{
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] *= s;
	}
	void operator*=(const double& s)
	{
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] *= s;
	}
	void operator/=(const T& s)
	{
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] /= s;
	}
	void operator+=(const MATRIX<T>& m)
	{
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] += m.A[i][j];
	}
	void operator-=(const MATRIX<T>& m)
	{
		for (int i = 0; i < H; i++)
			for (int j = 0; j < W; j++)
				A[i][j] -= m.A[i][j];
	}
	MATRIX<T> operator*(const MATRIX<T>& m)
	{
		if (W != m.H)throw exception("can not cross product");
		MATRIX<T> C(H, m.W);
		for (int i = 0; i < C.H; i++)
			for (int j = 0; j < C.W; j++)
				for (int k = 0; k < W; k++)
					C[i][j] += A[i][k] * m.A[k][j];
		return C;
	}
	MATRIX<T> conv(MATRIX<T> m)
	{
		MATRIX<T> C(H, W);
		int dx = m.W / 2;
		int dy = m.H / 2;
		int n = m.W * m.H;
		for (int y = dy; y < H - dy; y++)
			for (int x = dx; x < W - dx; x++)
			{
				C[y][x] = 0;
				for (int i = 0; i < m.H; i++)
					for (int j = 0; j < m.W; j++)
						C[y][x] += (A[y + i - dy][x + j - dx] * m(i + 1, j + 1));
			}
		return C;
	}
	friend ostream& operator<<(ostream& os, const MATRIX<T>& m)
	{
		for (int i = 0; i < m.H; i++)
		{
			for (int j = 0; j < m.W; j++)
				os << setw(5) << m.A[i][j];
			os << endl;
		}
		return os;
	}
	friend istream& operator>>(istream& is, MATRIX<T>& m)
	{
		for (int i = 0; i < m.H; i++)
			for (int j = 0; j < m.W; j++)
				is >> m.A[i][j];
		return is;
	}
};
class APP
{
private:
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	Gdiplus::Status status;
	MATRIX<int> A[3];
	MATRIX<int> G;
	int W = 0, H = 0;
	std::wstring wfilename;
	Bitmap* bmp = 0;
	int GetEncoderClsid(const WCHAR* format, CLSID* pClsid) {
		UINT num = 0, size = 0;
		GetImageEncodersSize(&num, &size);
		if (size == 0) return -1;

		ImageCodecInfo* pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
		if (!pImageCodecInfo) return -1;

		GetImageEncoders(num, size, pImageCodecInfo);
		for (UINT j = 0; j < num; j++) {
			if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0) {
				*pClsid = pImageCodecInfo[j].Clsid;
				free(pImageCodecInfo);
				return j;
			}
		}
		free(pImageCodecInfo);
		return -1;
	}
public:
	APP()
	{
		status = Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
	}
	~APP()
	{
		if (bmp)delete bmp;
		Gdiplus::GdiplusShutdown(gdiplusToken);
	}
	const int getH() { return H; }
	const int getW() { return W; }
	bool read(string filename)
	{
		Log("Bitmap", "read");
		wfilename = std::wstring(filename.begin(), filename.end());
		bmp = Bitmap::FromFile(wfilename.c_str());
		W = bmp->GetWidth();
		H = bmp->GetHeight();
		Color pixelColor;
		for (int i = 0; i < 3; i++)
			A[i].init(H, W);
		for (int y = 0; y < H; y++)
			for (int x = 0; x < W; x++) {
				bmp->GetPixel(x, y, &pixelColor);
				A[2](y + 1, x + 1) = pixelColor.GetRed();
				A[1](y + 1, x + 1) = pixelColor.GetGreen();
				A[0](y + 1, x + 1) = pixelColor.GetBlue();
			}
		return true;
	}
	bool write(string filename)
	{
		Log("Bitmap", "write");
		wstring wfilename = std::wstring(filename.begin(), filename.end());
		CreateAndSaveBitmap(wfilename, PixelFormat24bppRGB, W, H, L"image/bmp");
		return true;
	}
	CLSID GetEncoderClsid(const WCHAR* format) {
		UINT num = 0, size = 0;
		GetImageEncodersSize(&num, &size);
		if (size == 0) return CLSID{};
		ImageCodecInfo* pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
		if (pImageCodecInfo == nullptr) return CLSID{};
		GetImageEncoders(num, size, pImageCodecInfo);
		for (UINT i = 0; i < num; i++) {
			if (wcscmp(pImageCodecInfo[i].MimeType, format) == 0) {
				CLSID clsid = pImageCodecInfo[i].Clsid;
				free(pImageCodecInfo);
				return clsid;
			}
		}
		free(pImageCodecInfo);
		return CLSID{};
	}
	void CreateAndSaveBitmap(const std::wstring& filename, PixelFormat format, int width, int height, const std::wstring& mimeType)
	{
		Bitmap bmp(width, height, format);
		Rect rect(0, 0, width, height);
		BitmapData bmpData;
		bmp.LockBits(&rect, ImageLockModeWrite, format, &bmpData);
		BYTE* scan0 = static_cast<BYTE*>(bmpData.Scan0);
		int stride = bmpData.Stride;
		int bytesPerPixel = (format == PixelFormat24bppRGB) ? 3 : 4;
		for (int y = 0; y < height; y++) {
			BYTE* row = scan0 + y * stride;
			for (int x = 0; x < width; x++) {
				BYTE* pixel = row + x * bytesPerPixel;
				pixel[0] = A[0](y + 1, x + 1);
				pixel[1] = A[1](y + 1, x + 1);
				pixel[2] = A[2](y + 1, x + 1);
			}
		}
		bmp.UnlockBits(&bmpData);
		CLSID clsid = GetEncoderClsid(mimeType.c_str());
		if (clsid == CLSID{}) {
			std::wcerr << L"Failed to get encoder CLSID for " << mimeType << std::endl;
			return;
		}
		Status stat = bmp.Save(filename.c_str(), &clsid, nullptr);
		if (stat == Ok)
			std::wcout << L"Saved: " << filename << std::endl;
		else
			std::wcerr << L"Save failed: " << stat << std::endl;
	}
	//============= 期末處理 =============
	void fisheye()
	{
		Log("finalproject", "fisheye");
		MATRIX<int> out[3];
		for (int i = 0; i < 3; i++)
		{
			out[i].init(H, W);
			out[i].setAll(255);
		}
		double cx = W / 2.0, cy = H / 2.0;
		double radius = min(W, H) / 2.0;
		for (int y = 0; y < H; y++)
			for (int x = 0; x < W; x++) {
				double dx = x - cx;
				double dy = y - cy;
				double r = sqrt(dx * dx + dy * dy);
				if (r < radius) {
					double theta = atan2(dy, dx);

					double strength = 0.6; // 扭曲強度 0.0 ~ 1.0,
					double nr = r + (r * r / radius - r) * strength;

					int sx = int(cx + nr * cos(theta));
					int sy = int(cy + nr * sin(theta));
					if (sx >= 0 && sx < W && sy >= 0 && sy < H)
						for (int c = 0; c < 3; c++)
							out[c][y][x] = A[c][sy][sx];
				}
			}
		for (int i = 0; i < 3; i++) A[i] = out[i];
	}
	void mosaicCircle(int blockSize = 20, int dotRadius = 10)
	{
		Log("finalproject", "mosaicCircle");
		MATRIX<int> output[3];
		for (int c = 0; c < 3; c++)
		{
			output[c].init(H, W);
			output[c].setAll(255);
		}
		for (int by = 0; by < H; by += blockSize)
			for (int bx = 0; bx < W; bx += blockSize) {
				int r = 0, g = 0, b = 0, count = 0;
				for (int y = by; y < min(by + blockSize, H); y++)
					for (int x = bx; x < min(bx + blockSize, W); x++) {
						r += A[0][y][x];
						g += A[1][y][x];
						b += A[2][y][x];
						count++;
					}
				r /= count; g /= count; b /= count;
				drawCircle(output, bx + blockSize / 2, by + blockSize / 2, dotRadius, r, g, b);
			}
		for (int c = 0; c < 3; c++) A[c] = output[c];
	}
	void drawCircle(MATRIX<int> output[3], int cx, int cy, int radius, int r, int g, int b)
	{
		for (int y = max(cy - radius, 0); y <= min(cy + radius, H - 1); y++)
			for (int x = max(cx - radius, 0); x <= min(cx + radius, W - 1); x++) {
				int dx = x - cx, dy = y - cy;
				if (dx * dx + dy * dy <= radius * radius) {
					output[0][y][x] = r;
					output[1][y][x] = g;
					output[2][y][x] = b;
				}
			}
	}
	//-------------------------------------------
	void runAllPic(const string& effectName, const function<void()>& effectFunc, int count) // effect 函式指標
	{
		Log("finalproject", "runAllPic");
		char filename[100], output[100];
		for (int i = 0; i < count; i++)
		{
			sprintf_s(filename, "testpic%d.bmp", i);
			read(filename);
			effectFunc();
			sprintf_s(output, "testpic%d_%s.bmp", i, effectName.c_str());
			write(output);
		}
	}
	void finalproject() {
		Log("finalproject", "------------------------");
		runAllPic("fisheye", [this]() { fisheye(); }, 7);
		runAllPic("mosaicCircle", [this]() { mosaicCircle(20, 10); }, 7);
	}
	bool run()
	{
		if (status == Gdiplus::Ok) {
			finalproject();
			return true;
		}
		else
			return false;
	}
};
int main() {
	APP app;
	app.run();
	return 0;
}