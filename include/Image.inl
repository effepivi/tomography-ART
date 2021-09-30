#include <iostream>
#include <cmath>


//---------------------------------------------------
inline Image operator*(float k, const Image& anImage)
//---------------------------------------------------
{
    return anImage * k;
}


//--------------------
inline Image::Image():
//--------------------
    m_cols(0),
    m_rows(0),
    m_slices(0)
//--------------------
{}


//----------------------------------------
inline Image::Image(const Image& anImage):
//----------------------------------------
    m_cols(anImage.m_cols),
    m_rows(anImage.m_rows),
    m_slices(anImage.m_slices),
    m_p_voxel_set(anImage.m_p_voxel_set)
//----------------------------------------
{}


//-----------------------------------------------------------------------
inline Image::Image(int cols, int rows, int slices, float aDefaultValue):
//-----------------------------------------------------------------------
    m_cols(cols),
    m_rows(rows),
    m_slices(slices),
    m_p_voxel_set(vector<float>(m_cols * m_rows * m_slices, aDefaultValue))
//-----------------------------------------------------------------------
{}


//----------------------------------------
inline Image::Image(const char* aFileName)
//----------------------------------------
{
    loadTIFF(aFileName);
}


//------------------------------------------
inline Image::Image(const string& aFileName)
//------------------------------------------
{
    loadTIFF(aFileName);
}


//--------------------------------------------------
inline void Image::loadTIFF(const string& aFileName)
//--------------------------------------------------
{
    loadTIFF(aFileName.c_str());
}

//--------------------------------------------------------
inline void Image::saveTIFF(const string& aFileName) const
//--------------------------------------------------------
{
    saveTIFF(aFileName.c_str());
}

//--------------------------------------------------------
inline void Image::saveTEXT(const string& aFileName) const
//--------------------------------------------------------
{
    saveTEXT(aFileName.c_str());
}


//------------------------------------------------------------------------------
inline void Image::allocate(int cols, int rows, int slices, float aDefaultValue)
//------------------------------------------------------------------------------
{
    m_cols = cols;
    m_rows = rows;
    m_slices = slices;
    m_p_voxel_set = vector<float>(m_cols * m_rows * m_slices, aDefaultValue);
}


//------------------------
inline void Image::clear()
//------------------------
{
    m_cols = 0;
    m_rows = 0;
    m_slices = 0;
    m_p_voxel_set.clear();
}


//----------------------------
inline int Image::size() const
//----------------------------
{
    return m_p_voxel_set.size();

    // same as
    // return m_cols * m_rows * m_slices;
}


//--------------------------------------------------
inline Image& Image::operator=(const Image& anImage)
//--------------------------------------------------
{
    m_cols = anImage.m_cols;
    m_rows = anImage.m_rows;
    m_slices = anImage.m_slices;
    m_p_voxel_set = anImage.m_p_voxel_set;

    return *this;
}


//---------------------------------------------------
inline Image& Image::operator+=(const Image& anImage)
//---------------------------------------------------
{
    // Deal with images of different sizes
    unsigned int min_width( std::min(m_cols,  anImage.m_cols));
    unsigned int min_height(std::min(m_rows, anImage.m_rows));
    unsigned int min_depth( std::min(m_slices,  anImage.m_slices));

    for (int z = 0; z < min_depth; ++z)
        for (int y = 0; y < min_height; ++y)
            for (int x = 0; x < min_width; ++x)
            {
                getPixel(x, y, z) += anImage(x, y, z);
            }

    return *this;
}


//---------------------------------------------------
inline Image& Image::operator*=(const Image& anImage)
//---------------------------------------------------
{
    // Deal with images of different sizes
    unsigned int min_width( std::min(m_cols,  anImage.m_cols));
    unsigned int min_height(std::min(m_rows, anImage.m_rows));
    unsigned int min_depth( std::min(m_slices,  anImage.m_slices));

    for (int z = 0; z < min_depth; ++z)
        for (int y = 0; y < min_height; ++y)
            for (int x = 0; x < min_width; ++x)
            {
                getPixel(x, y, z) *= anImage(x, y, z);
            }

    return *this;
}


//--------------------------------------
inline Image& Image::operator*=(float k)
//--------------------------------------
{
    for (auto& ite : m_p_voxel_set)
    {
        ite *= k;
    }

    return *this;
}


//--------------------------------------
inline Image& Image::operator/=(float k)
//--------------------------------------
{
    for (auto& ite : m_p_voxel_set)
    {
        ite /= k;
    }

    return *this;
}


//---------------------------------------
inline Image& Image::operator*=(double k)
//---------------------------------------
{
    for (auto& ite : m_p_voxel_set)
    {
        ite *= k;
    }

    return *this;
}


//---------------------------------------
inline Image& Image::operator/=(double k)
//---------------------------------------
{
    for (auto& ite : m_p_voxel_set)
    {
        ite /= k;
    }

    return *this;
}


//-------------------------------------------------------
inline Image Image::operator-(const Image& anImage) const
//-------------------------------------------------------
{
    // Deal with images of different sizes
    unsigned int min_width( std::min(m_cols,  anImage.m_cols));
    unsigned int min_height(std::min(m_rows, anImage.m_rows));
    unsigned int min_depth( std::min(m_slices,  anImage.m_slices));

    Image temp(min_width, min_height, min_depth, 0.0);

    for (int z = 0; z < min_depth; ++z)
        for (int y = 0; y < min_height; ++y)
            for (int x = 0; x < min_width; ++x)
            {
                temp(x, y, z) = getPixel(x, y, z) - anImage(x, y, z);
            }

    return temp;
}


//-------------------------------------------------------
inline Image Image::operator*(const Image& anImage) const
//-------------------------------------------------------
{
    // Deal with images of different sizes
    unsigned int min_width( std::min(m_cols,  anImage.m_cols));
    unsigned int min_height(std::min(m_rows, anImage.m_rows));
    unsigned int min_depth( std::min(m_slices,  anImage.m_slices));

    Image temp(min_width, min_height, min_depth, 0.0);

    for (int z = 0; z < min_depth; ++z)
        for (int y = 0; y < min_height; ++y)
            for (int x = 0; x < min_width; ++x)
            {
                temp(x, y, z) = getPixel(x, y, z) * anImage(x, y, z);
            }

    return temp;
}


//-------------------------------------------------------
inline Image Image::operator/(const Image& anImage) const
//-------------------------------------------------------
{
    // Deal with images of different sizes
    unsigned int min_width( std::min(m_cols,  anImage.m_cols));
    unsigned int min_height(std::min(m_rows, anImage.m_rows));
    unsigned int min_depth( std::min(m_slices,  anImage.m_slices));

    Image temp(min_width, min_height, min_depth, 0.0);

    for (int z = 0; z < min_depth; ++z)
        for (int y = 0; y < min_height; ++y)
            for (int x = 0; x < min_width; ++x)
            {
                temp(x, y, z) = getPixel(x, y, z) / anImage(x, y, z);
            }

    return temp;
}


//------------------------------------------
inline Image Image::operator-(float k) const
//------------------------------------------
{
    Image temp = *this;

    for (auto& ite : temp.m_p_voxel_set)
    {
        ite -= k;
    }

    return temp;
}


//------------------------------------------
inline Image Image::operator*(float k) const
//------------------------------------------
{
    Image temp = *this;

    for (auto& ite : temp.m_p_voxel_set)
    {
        ite *= k;
    }

    return temp;
}


//------------------------------------------
inline Image Image::operator/(float k) const
//------------------------------------------
{
    Image temp = *this;

    for (auto& ite : temp.m_p_voxel_set)
    {
        ite /= k;
    }

    return temp;
}


//-------------------------------------------
inline Image Image::operator*(double k) const
//-------------------------------------------
{
    Image temp = *this;

    for (auto& ite : temp.m_p_voxel_set)
    {
        ite *= k;
    }

    return temp;
}


//-------------------------------------------
inline Image Image::operator/(double k) const
//-------------------------------------------
{
    Image temp = *this;

    for (auto& ite : temp.m_p_voxel_set)
    {
        ite /= k;
    }

    return temp;
}


//------------------------------------------------
inline float& Image::getPixel(int x, int y, int z)
//------------------------------------------------
{
    if (x < 0 || x > m_cols || y < 0 || y > m_rows || z < 0 || z > m_slices)
    {
        cerr << "Invalid pixel index." << endl;
        exit(EXIT_FAILURE);
    }

    return m_p_voxel_set[z * m_cols * m_rows + y * m_cols + x];
}


//------------------------------------------------------------
inline const float& Image::getPixel(int x, int y, int z) const
//------------------------------------------------------------
{
    if (x < 0 || x > m_cols || y < 0 || y > m_rows || z < 0 || z > m_slices)
    {
        cerr << "Invalid pixel index." << endl;
        exit(EXIT_FAILURE);
    }

    return m_p_voxel_set[z * m_cols * m_rows + y * m_cols + x];
}


//--------------------------------------------------
inline float& Image::operator()(int x, int y, int z)
//--------------------------------------------------
{
    return getPixel(x, y, z);
}


//--------------------------------------------------------------
inline const float& Image::operator()(int x, int y, int z) const
//--------------------------------------------------------------
{
    return getPixel(x, y, z);
}


//-------------------------------
inline int Image::getCols() const
//-------------------------------
{
    return m_cols;
}


//-------------------------------
inline int Image::getRows() const
//-------------------------------
{
    return m_rows;
}


//---------------------------------
inline int Image::getSlices() const
//---------------------------------
{
    return m_slices;
}


//---------------------------------------------------
inline Image Image::getROI(unsigned int i,
                           unsigned int j,
                           unsigned int k,
                           unsigned int aWidth,
                           unsigned int aHeight,
                           unsigned int aDepth) const
//---------------------------------------------------
{
    // Create a black image
    Image roi(aWidth, aHeight, aDepth);


    // Process every slice of the ROI
	for (int z = 0; z < aDepth; ++z)
    {
        // Process every row of the ROI
        for (int y = 0; y < aHeight; ++y)
        {
            // Process every column of the ROI
            for (int x = 0; x < aWidth; ++x)
            {
                unsigned int index_i(x + i);
                unsigned int index_j(y + j);
                unsigned int index_k(z + k);

                // The pixel index is not valid
                if ((index_i >= m_cols) ||
                        (index_j >= m_rows) ||
                        (index_k >= m_slices))
                {
                    cerr << "Invalid pixel index" << endl;
                    exit(EXIT_FAILURE);
                }

                // Get the pixel intensity from the current instance
                float intensity(getPixel(index_i, index_j, index_k));

                // Set the pixel of the ROI
                roi(x, y, z) = intensity;
            }
        }
    }

    return (roi);
}


//----------------------------------
inline double Image::getMean() const
//----------------------------------
{
    double accu = 0.0;

    int npixel = m_p_voxel_set.size();
    for (int i = 0; i < npixel; ++i)
    {
        accu += m_p_voxel_set[i];
    }

    return accu / npixel;
}


//---------------------------------
inline double Image::getStd() const
//---------------------------------
{
    double variance = 0.0;
    double avg = getMean();

    int npixel = m_p_voxel_set.size();
    for (int i = 0; i < npixel; ++i)
    {
        variance += std::pow(m_p_voxel_set[i] - avg, 2);
    }

    variance /= npixel;

    return std::sqrt(variance);
}


//---------------------------------------------------
inline double Image::rmse(const Image& anImage) const
//---------------------------------------------------
{
    return std::sqrt((*this - anImage).square().getMean());
}


//---------------------------------------------------
inline double Image::zncc(const Image& anImage) const
//---------------------------------------------------
{
    Image temp1 = *this - getMean();
    Image temp2 = anImage - anImage.getMean();

    temp1 /= getStd();
    temp2 /= anImage.getStd();

    return (temp1 * temp2).getMean();
}


//--------------------------------
inline Image Image::square() const
//--------------------------------
{
    Image temp = *this;
    temp *= temp;

    return temp;
}
