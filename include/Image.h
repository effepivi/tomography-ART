#ifndef IMAGE_H
#define IMAGE_H

#include <string>
#include <vector>

using namespace std;


class Image;


Image operator*(float k, const Image& anImage);


class Image
{
public:
    Image();

    Image(const Image& anImage);

    Image(int cols, int rows, int slices, float aDefaultValue = 0.0);

    Image(const char* aFileName);
    Image(const string& aFileName);

    void loadTIFF(const char* aFileName);
    void loadTIFF(const string& aFileName);

    void saveTIFF(const char* aFileName) const;
    void saveTIFF(const string& aFileName) const;

    void saveTEXT(const char* aFileName) const;
    void saveTEXT(const string& aFileName) const;

    void allocate(int cols, int rows, int slices, float aDefaultValue = 0.0);
    int size() const;
    void clear();

    Image& operator=(const Image& anImage);

    Image& operator+=(const Image& anImage);
    Image& operator*=(const Image& anImage);

    Image& operator*=(float k);
    Image& operator/=(float k);
    Image& operator*=(double k);
    Image& operator/=(double k);

    Image operator-(const Image& anImage) const;
    Image operator*(const Image& anImage) const;
    Image operator/(const Image& anImage) const;

    Image operator-(float k) const;
    Image operator*(float k) const;
    Image operator/(float k) const;
    Image operator*(double k) const;
    Image operator/(double k) const;

    Image radon(const vector<double>& aSetOfAnglesInDegrees, int aWidth = 0) const;
    Image iradon(const vector<double>& aSetOfAnglesInDegrees) const;

    Image iradonART(const vector<double>& aSetOfAnglesInDegrees,
        double lambda, int n) const;

    Image rotate(double anAngleInDegrees) const;

    Image setCanvasSize(int cols, int rows, int slices, float aDefaultValue = 0.0) const;

    float& getPixel(int x, int y, int z = 0);
    const float& getPixel(int x, int y, int z = 0) const;

    float& operator()(int x, int y, int z = 0);
    const float& operator()(int x, int y, int z = 0) const;

    int getCols() const;
    int getRows() const;
    int getSlices() const;

    Image getROI(unsigned int i,
        unsigned int j,
        unsigned int k,
        unsigned int aWidth,
        unsigned int aHeight,
        unsigned int aDepth) const;

    double getMean() const;
    double getStd() const;

    double rmse(const Image& anImage) const;
    double zncc(const Image& anImage) const;

    Image square() const;

private:
    int m_cols;
    int m_rows;
    int m_slices;
    vector<float> m_p_voxel_set;
};


#include "Image.inl"


#endif
