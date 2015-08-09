#include "share.h"
#include "brent.h"
#include "random.h"
#include "sampler.h"
#include "GaussianFiber.h"
#include "sggx.h"
#include <fstream>
class LowDiscrepancySampler : public Sampler {
public:
	//LowDiscrepancySampler() : Sampler(Properties()) { }

	LowDiscrepancySampler(/*const Properties &props*/) /*: Sampler(props) */{
		/* Sample count (will be rounded up to the next power of two) */
		m_sampleCount = 4;

		/* Dimension, up to which which low discrepancy samples are guaranteed to be available. */
		m_maxDimension = 4;

		/*	if (!math::isPowerOfTwo(m_sampleCount)) {
		m_sampleCount = math::roundToPowerOfTwo(m_sampleCount);
		Log(EWarn, "Sample count should be a power of two -- rounding to "
		SIZE_T_FMT, m_sampleCount);
		}*/

		m_samples1D = new float*[m_maxDimension];
		m_samples2D = new Point2*[m_maxDimension];

		for (size_t i = 0; i<m_maxDimension; i++) {
			m_samples1D[i] = new float[m_sampleCount];
			m_samples2D[i] = new Point2[m_sampleCount];
		}

		m_random = new Random();
	}

	/*LowDiscrepancySampler(Stream *stream, InstanceManager *manager)
	: Sampler(stream, manager) {
	m_random = static_cast<Random *>(manager->getInstance(stream));
	m_maxDimension = stream->readSize();

	m_samples1D = new Float*[m_maxDimension];
	m_samples2D = new Point2*[m_maxDimension];
	for (size_t i = 0; i<m_maxDimension; i++) {
	m_samples1D[i] = new Float[(size_t)m_sampleCount];
	m_samples2D[i] = new Point2[(size_t)m_sampleCount];
	}
	}*/

	virtual ~LowDiscrepancySampler() {
		for (size_t i = 0; i<m_maxDimension; i++) {
			delete[] m_samples1D[i];
			delete[] m_samples2D[i];
		}
		delete[] m_samples1D;
		delete[] m_samples2D;
		delete m_random;
	}
	inline void generate1D(float *samples, size_t sampleCount) {

		uint32_t scramble = m_random->nextULong() & 0xFFFFFFFF;
		for (size_t i = 0; i < sampleCount; ++i)
			samples[i] = radicalInverse2Single((uint32_t)i, scramble);
		m_random->shuffle(samples, samples + sampleCount);
	}

	inline void generate2D(Point2 *samples, size_t sampleCount) {
		union {
			uint64_t qword;
			uint32_t dword[2];
		} scramble;

		scramble.qword = m_random->nextULong();

		for (size_t i = 0; i < sampleCount; ++i)
			samples[i] = sample02Single((uint32_t)i, scramble.dword);
		m_random->shuffle(samples, samples + sampleCount);
	}

	void generate(const Point2 &) {
		for (size_t i = 0; i<m_maxDimension; ++i) {
			generate1D(m_samples1D[i], m_sampleCount);
			generate2D(m_samples2D[i], m_sampleCount);
		}

		for (size_t i = 0; i<m_req1D.size(); i++)
			generate1D(m_sampleArrays1D[i], m_sampleCount * m_req1D[i]);

		for (size_t i = 0; i<m_req2D.size(); i++)
			generate2D(m_sampleArrays2D[i], m_sampleCount * m_req2D[i]);

		m_sampleIndex = 0;
		m_dimension1D = m_dimension2D = 0;
		m_dimension1DArray = m_dimension2DArray = 0;
	}

	void advance() {
		m_sampleIndex++;
		m_dimension1D = m_dimension2D = 0;
		m_dimension1DArray = m_dimension2DArray = 0;
	}

	void setSampleIndex(size_t sampleIndex) {
		m_sampleIndex = sampleIndex;
		m_dimension1D = m_dimension2D = 0;
		m_dimension1DArray = m_dimension2DArray = 0;
	}

	float next1D() {
		//	Assert(m_sampleIndex < m_sampleCount);
		if (m_dimension1D < m_maxDimension)
			return m_samples1D[m_dimension1D++][m_sampleIndex];
		else
			return m_random->nextFloat();
	}

	Point2 next2D() {
		//Assert(m_sampleIndex < m_sampleCount);
		if (m_dimension2D < m_maxDimension)
			return m_samples2D[m_dimension2D++][m_sampleIndex];
		else
			return Point2(m_random->nextFloat(), m_random->nextFloat());
	}

	/*std::string toString() const {
	std::ostringstream oss;
	oss << "LowDiscrepancySampler[" << endl
	<< "  sampleCount = " << m_sampleCount << "," << endl
	<< "  dimension = " << m_maxDimension << endl
	<< "]";
	return oss.str();
	}*/

	//MTS_DECLARE_CLASS()
public:
	Random* m_random;
private:

	size_t m_maxDimension;
	size_t m_dimension1D;
	size_t m_dimension2D;
	float **m_samples1D;
	Point2 **m_samples2D;
};
/**
* \brief Flake distribution for simulating rough fibers
*
* This class implements the Gaussian flake distribution proposed in
*
* "Building Volumetric Appearance Models of Fabric using
* Micro CT Imaging" by Shuang Zhao, Wenzel Jakob, Steve Marschner,
* and Kavita Bala, ACM SIGGRAPH 2011
*
* \author Wenzel Jakob
*/

float m_cosPhi[256], m_sinPhi[256];
float m_cosTheta[256], m_sinTheta[256];
int Sum = 0;
int readInt(std::fstream &f){
	char tmpint[4];
	f.read(tmpint, 4);
	int ret = *((int*)tmpint);
	return ret;
}
float readFloat(std::fstream &f){
	char tmpfloat[4];
	f.read(tmpfloat, 4);
	float ret = *((float*)tmpfloat);
	return ret;
}
inline void writeInt(std::fstream &f, int &a){
	f.write((char *)(&a), 4);
}
inline void writeFloat(std::fstream &f, float &a){
	f.write((char *)(&a), 4);
}
inline void writeVec(std::fstream &f, vec3 &v){
	writeFloat(f, v.x);
	writeFloat(f, v.y);
	writeFloat(f, v.z);
}
float TMPS[300000][6];
void readVolumn(char* filename, char* nfilename){
	std::fstream f;
	std::fstream fo;
	f.open(filename, std::ios::binary | std::ios::in);
	fo.open(nfilename, std::ios::binary | std::ios::out);
	//fo2.open(nfilename2, std::ios::binary | std::ios::out);
	printf("%s\n", filename);
	if (f.fail())
		printf("Fail!\n");
	Sum++;
	char s[3];
	unsigned char v;
	int encode, xcell, ycell, zcell, channel;
	float xmin, ymin, zmin, xmax, ymax, zmax;
	f.read(s, 3);
	fo.write(s, 3);
	f.read((char *)&v, 1);
	fo.write((char *)&v, 1);
	encode = readInt(f);
	encode = 5;
	writeInt(fo, encode);
	xcell = readInt(f);
	ycell = readInt(f);
	zcell = readInt(f);

	int lod = 2;
	int len = 1 << lod;
	int nxcell = xcell >> lod, nycell = ycell >> lod, nzcell = zcell >> lod;
	writeInt(fo, nxcell); writeInt(fo, nycell); writeInt(fo, nzcell);

	channel = readInt(f);
	channel = 6;
	writeInt(fo, channel);
	xmin = readFloat(f);
	ymin = readFloat(f);
	zmin = readFloat(f);
	xmax = readFloat(f);
	ymax = readFloat(f);
	zmax = readFloat(f);
	writeFloat(fo, xmin); writeFloat(fo, ymin); writeFloat(fo, zmin);
	writeFloat(fo, xmax); writeFloat(fo, ymax); writeFloat(fo, zmax);
	SggxPhaseFunction s1(0.3);

	for (int k = 0; k < zcell; k++)
		for (int j = 0; j < ycell; j++)
			for (int i = 0; i < xcell; i++)
			{
				vec3 ori;
				unsigned char theta, phi;
				int index = (k*ycell + j)*xcell + i;
				f.read((char*)&theta, 1);
				f.read((char*)&phi, 1);
				ori.x = m_cosPhi[phi] * m_sinTheta[theta];
				ori.y = m_sinPhi[phi] * m_sinTheta[theta];
				ori.z = m_cosTheta[theta];
				assert(!f.eof());

				//Sxx:0 Syy:1 Szz:2 Sxy:3 Sxz:4 Syz:5
				for (int l = 0; l<6; l++) TMPS[index][l] = 0.0;
				float length = dot(ori, ori);
				if (fabs(length) > 1e-5)
					s1.GetS(ori, TMPS[index]);
				for (int l = 0; l < 6; l++)
					assert(TMPS[index][l] < 1.001);
	
						
				/*float S[6];
				for (int l = 0; l < 6; l++)
					S[l] = TMPS[index][l];
				if (S[0] * S[0] + S[1] * S[1] + S[2] * S[2] + S[3] * S[3] + S[4] * S[4] + S[5] * S[5] > 1e-5)
				{
					float max = 0.0;
					for (int j = 0; j <= 360; j++)
					{
						for (int i = 0; i <= 180; i++)
						{
							float radi = 2 * M_PI*(1.0*i / 360);
							float radj = 2 * M_PI*(1.0*j / 360);
							vec3 ori = vec3(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(radi));
							float value = sigma(ori, S[0], S[1], S[2], S[3], S[4], S[5]);
							if (value > max)
								max = value;
						}

					}
					if (max > 1.0001)
					{
						printf("guale %7f\n", max);
						system("pause");
					}						
				}*/

			}
	for (int k = 0; k < nzcell; k++)
		for (int j = 0; j < nycell; j++)
			for (int i = 0; i < nxcell; i++)
			{
				if (k == 18 && i == 0 && j == 18)
					k = k*1;
				float S[6];
				memset(S, 0, sizeof(S));
				int bx = i << lod, by = j << lod, bz = k << lod;
				for (int nk = bz; nk < bz + len; nk++)
					for (int nj = by; nj<by + len; nj++)
						for (int ni = bx; ni<bx + len; ni++)
						{
							int index = (nk*ycell + nj)*xcell + ni;
							for (int l = 0; l < 6; l++)
								S[l] += TMPS[index][l];
						}
				for (int l = 0; l < 6; l++)
					S[l] = S[l] / (len*len*len);
				for (int l = 0; l < 6; l++)
					assert(S[l] < 1.001);
				vec3 Sigma, r;
				if (S[0] * S[0] + S[1] * S[1] + S[2] * S[2] + S[3] * S[3] + S[4] * S[4] + S[5] * S[5]> 1e-5)
				{
					Sigma.x = sqrtf(S[0]); Sigma.y = sqrtf(S[1]); Sigma.z = sqrtf(S[2]);
					r.x = S[3] / sqrtf(S[0] * S[1]);
					r.y = S[4] / sqrtf(S[0] * S[2]);
					r.z = S[5] / sqrtf(S[1] * S[2]);
				/*	float max = 0.0;
					for (int jj = 0; jj <= 360; jj++)
					{
						for (int ii = 0; ii <= 180; ii++)
						{
							float radi = 2 * M_PI*(1.0*ii / 360);
							float radj = 2 * M_PI*(1.0*jj / 360);
							vec3 ori = vec3(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(radi));
							float value = sigma(ori, S[0], S[1], S[2], S[3], S[4], S[5]);
							if (value > max)
								max = value;
						}

					}
					printf("%f\n", max);
					system("pause");*/
				}
	
				float rate = 255.0;
				assert(Sigma.x*rate >= 0 && Sigma.x*rate <= 255);
				assert(Sigma.y*rate >= 0 && Sigma.y*rate <= 255);
				assert(Sigma.z*rate >= 0 && Sigma.z*rate <= 255);
				assert((r.x + 1)*0.5*rate >= 0 && (r.x + 1)*0.5*rate <= 255);
				assert((r.y + 1)*0.5*rate >= 0 && (r.y + 1)*0.5*rate <= 255);
				assert((r.z + 1)*0.5*rate >= 0 && (r.z + 1)*0.5*rate <= 255);
				unsigned char sx = Sigma.x*rate, sy = Sigma.y*rate, sz = Sigma.z*rate;
				unsigned char rx = (r.x + 1)*0.5*rate, ry = (r.y + 1)*0.5*rate, rz = (r.z + 1)*0.5*rate;
				fo.write((char*)&sx, 1); fo.write((char*)&sy, 1); fo.write((char*)&sz, 1);
				fo.write((char*)&rx, 1); fo.write((char*)&ry, 1); fo.write((char*)&rz, 1);
			}
	printf("x:%f->%f %d\n", xmin, xmax, xcell);
	printf("y:%f->%f %d\n", ymin, ymax, ycell);
	printf("z:%f->%f %d\n", zmin, zmax, zcell);
	assert((xcell == 64) && (ycell == 64) && (zcell == 64));
	f.close();

	//system("pause");
}
unsigned char TMPD[300000];
void downsampleDensity(char* filename, char* nfilename){
	std::fstream f;
	std::fstream fo;
	f.open(filename, std::ios::binary | std::ios::in);
	fo.open(nfilename, std::ios::binary | std::ios::out);
	//fo2.open(nfilename2, std::ios::binary | std::ios::out);
	printf("%s\n", filename);
	if (f.fail())
		printf("Fail!\n");
	Sum++;
	char s[3];
	unsigned char v;
	int encode, xcell, ycell, zcell, channel;
	float xmin, ymin, zmin, xmax, ymax, zmax;
	f.read(s, 3);
	fo.write(s, 3);
	f.read((char *)&v, 1);
	fo.write((char *)&v, 1);
	encode = readInt(f);
	encode = 3;
	writeInt(fo, encode);
	xcell = readInt(f);
	ycell = readInt(f);
	zcell = readInt(f);

	int lod = 2;
	int len = 1 << lod;
	int nxcell = xcell >> lod, nycell = ycell >> lod, nzcell = zcell >> lod;
	writeInt(fo, nxcell); writeInt(fo, nycell); writeInt(fo, nzcell);

	channel = readInt(f);
	channel = 1;
	writeInt(fo, channel);
	xmin = readFloat(f);
	ymin = readFloat(f);
	zmin = readFloat(f);
	xmax = readFloat(f);
	ymax = readFloat(f);
	zmax = readFloat(f);
	writeFloat(fo, xmin); writeFloat(fo, ymin); writeFloat(fo, zmin);
	writeFloat(fo, xmax); writeFloat(fo, ymax); writeFloat(fo, zmax);
	

	for (int k = 0; k < zcell; k++)
		for (int j = 0; j < ycell; j++)
			for (int i = 0; i < xcell; i++)
			{
				vec3 ori;
				int index = (k*ycell + j)*xcell + i;
				f.read((char*)&TMPD[index], 1);
			}
	for (int k = 0; k < nzcell; k++)
		for (int j = 0; j < nycell; j++)
			for (int i = 0; i < nxcell; i++)
			{
			
				int sumd=0;
				int bx = i << lod, by = j << lod, bz = k << lod;
				for (int nk = bz; nk < bz + len; nk++)
					for (int nj = by; nj<by + len; nj++)
						for (int ni = bx; ni<bx + len; ni++)
						{
							int index = (nk*ycell + nj)*xcell + ni;
							sumd += TMPD[index];
						}
				if (sumd != 0)
					sumd *= 1;
				sumd = int(sumd*1.0 / (len*len*len)+0.50);
				assert(sumd >= 0 && sumd <= 255);
				unsigned char o = sumd;
				fo.write((char *)&o, 1);
			}
	printf("x:%f->%f %d\n", xmin, xmax, xcell);
	printf("y:%f->%f %d\n", ymin, ymax, ycell);
	printf("z:%f->%f %d\n", zmin, zmax, zcell);
	assert((xcell == 64) && (ycell == 64) && (zcell == 64));
	f.close();

	//system("pause");
}
int main(){
	//GaussianFiberDistribution g1(0.3);
	/*	Random* m_random = new Random();
	Sampler* m_sampler = new LowDiscrepancySampler();
	SggxPhaseFunction s1(0.1);
	m_sampler->generate(Point2(0,0));*/
	/*	float S[5][6];
	SggxPhaseFunction s1(0.1);
	s1.GetS(vec3(1, 0, 1), S[0]);
	s1.GetS(vec3(0.7, 0, 1), S[1]);
	s1.GetS(vec3(0, 0.3, 1), S[2]);
	s1.GetS(vec3(0.8, -0.3, 1), S[3]);
	for (int i = 0; i < 6; i++)
	S[4][i] = 0.25*(S[0][i] + S[1][i] + S[2][i] + S[3][i]);
	for (int k = 0; k < 5; k++)
	{
	float max = 0.0;
	for (int j = 0; j <= 360; j++)
	{
	for (int i = 0; i <= 180; i++)
	{
	float radi = 2 * M_PI*(1.0*i / 360);
	float radj = 2 * M_PI*(1.0*j / 360);
	vec3 ori = vec3(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(radi));
	float value = sigma(ori, S[k][0], S[k][1], S[k][2], S[k][3], S[k][4], S[k][5]);
	if (value > max)
	max = value;
	}

	}
	printf("%d: %f\n", k, max);
	}*/
	//s1.OutputDistribution(S[0], "s0.txt");
	//s1.OutputDistribution(S[1], "s1.txt");
	//s1.OutputDistribution(S[2], "s2.txt");
	//s1.OutputDistribution(S[3], "s3.txt");

	//s1.OutputDistribution(S[4], "s4.txt");

//	system("pause");

	for (int i = 0; i<255; i++) {
		float angle = (float)i * ((float)M_PI / 255.0f);
		m_cosPhi[i] = cosf(2.0f * angle);
		m_sinPhi[i] = sinf(2.0f * angle);
		m_cosTheta[i] = cosf(angle);
		m_sinTheta[i] = sinf(angle);
	}
	m_cosPhi[255] = m_sinPhi[255] = 0;
	m_cosTheta[255] = m_sinTheta[255] = 0;


	std::fstream f;
	f.open("data/volume_description.vol", std::ios::binary | std::ios::in);
	float xmin, ymin, zmin, xmax, ymax, zmax;
	int xgrid, ygrid, zgrid;
	xmin = readFloat(f);
	ymin = readFloat(f);
	zmin = readFloat(f);
	xmax = readFloat(f);
	ymax = readFloat(f);
	zmax = readFloat(f);
	xgrid = readInt(f);
	ygrid = readInt(f);
	zgrid = readInt(f);
	char filenameo[100],filenamed[100], nfilenames[100], nfilenamed[100];
	char prefix[20] = "data/volume_";
	char postfixo[20] = "-orientation.vol";
	char postfixd[20] = "-density.vol";
	char nprefix[20] = "data/volume_";
	char npostfixs[20] = "-sggx_lod2.vol";
	char npostfixd[20] = "-density_lod2.vol";
	while (1){
		int x, y, z;
		x = readInt(f);
		if (f.eof()) break;
		y = readInt(f);
		z = readInt(f);
		sprintf_s(filenameo, "%s%03i_%03i_%03i%s", prefix, x, y, z, postfixo);
		sprintf_s(nfilenames, "%s%03i_%03i_%03i%s", nprefix, x, y, z, npostfixs);
		sprintf_s(filenamed, "%s%03i_%03i_%03i%s", prefix, x, y, z, postfixd);
		sprintf_s(nfilenamed, "%s%03i_%03i_%03i%s", nprefix, x, y, z, npostfixd);
		readVolumn(filenameo, nfilenames);
		//downsampleDensity(filenamed, nfilenamed);
		//	system("pause");
	}
	f.close();
	printf("x:%d y:%d z:%d sum:%d\n", xgrid, ygrid, zgrid, Sum);

	system("pause");
}