#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "Image.h"

#ifdef __APPLE__
typedef unsigned int TIFF_UINT32_T;
typedef unsigned short TIFF_UINT16_T;
#endif

#define DEBUG


//-------------------------------------------------------------------------------
Image Image::radon(const vector<double>& aSetOfAnglesInDegrees, int aWidth) const
//-------------------------------------------------------------------------------
{
    unsigned int sinogram_cols = aWidth;

    if (aWidth == 0)
    {
        unsigned int square_width;
        square_width = 4 + std::max(m_cols, m_rows);
        square_width *= square_width;

        sinogram_cols = floor(std::sqrt(2.0 * square_width));
    }

    double scaling_factor(double(sinogram_cols) / double(m_cols));

    Image sinogram(sinogram_cols, aSetOfAnglesInDegrees.size(), m_slices);

    Image padded_image(setCanvasSize(sinogram_cols, scaling_factor * m_rows, m_slices));

    for (unsigned int angle_id = 0; angle_id < aSetOfAnglesInDegrees.size(); ++angle_id)
    {
        Image rotated_image(padded_image.rotate(aSetOfAnglesInDegrees[angle_id]));

// #ifdef DEBUG
//         stringstream filename;
//         filename.fill('0');
//         filename << "rotated_image_" << std::setw(3) << angle_id << ".tif";
//         rotated_image.saveTIFF(filename.str());
// #endif

		for (int j = 0; j < rotated_image.m_rows; ++j)
        {
#pragma omp parallel for collapse(2)
            for (int i = 0; i < rotated_image.m_cols; ++i)
            {
                for (int k = 0; k < m_slices; ++k)
                {
                    sinogram(i, angle_id, k) += rotated_image(i, j, k);
                }
            }
        }
    }

    return (sinogram);
}


//--------------------------------------------------------------------
Image Image::iradon(const vector<double>& aSetOfAnglesInDegrees) const
//--------------------------------------------------------------------
{
    Image reconstruction(getCols(), getCols(), getSlices());

    double angle(0.0);


    for (unsigned int angle_id(0);
            angle_id < aSetOfAnglesInDegrees.size();
            ++angle_id)
    {
        double angle = aSetOfAnglesInDegrees[angle_id];

        Image back_projection(getCols(), getCols(), getSlices(), 0.0);

        int nslices = reconstruction.getSlices();
        int nrows = reconstruction.getRows();
        int ncols = reconstruction.getCols();

#pragma omp parallel for collapse(2)
        for (int k = 0; k < nslices; ++k)
        {
            for (int j = 0; j < nrows; ++j)
            {
                for (int i = 0; i < ncols; ++i)
                {
                    back_projection(i, j, k) += getPixel(i, angle_id, k);
                }
            }
        }

        reconstruction += back_projection.rotate(-angle);
    }

    //reconstruction = reconstruction.setCanvasSize(image_width, image_width, this->m_slices);
    reconstruction *= M_PI / (2.0 * aSetOfAnglesInDegrees.size());

    return reconstruction;
}


//-----------------------------------------------------------------
Image Image::iradonART(const vector<double>& aSetOfAnglesInDegrees,
                       double lambda, int n) const
//-----------------------------------------------------------------
{
    // Create an homogeneous image
    Image reconstruction(m_cols, m_cols, m_slices, 0.0);

    Image ones(m_cols, m_cols, m_slices, 1.0);

    Image AT = ones.radon(aSetOfAnglesInDegrees, getCols());
    AT /= M_PI / (2.0 * aSetOfAnglesInDegrees.size());

    Image ATA = AT.iradon(aSetOfAnglesInDegrees);

    // Process every iterations
    for (int i = 0; i < n; ++i)
    {
        cout << "Iteration " << i + 1 << " of " << n << endl;
        Image estimated_projections = reconstruction.radon(aSetOfAnglesInDegrees, getCols());
        Image error = *this - estimated_projections;
        Image reconstructed_error = error.iradon(aSetOfAnglesInDegrees);
        reconstructed_error /= M_PI / (2.0 * aSetOfAnglesInDegrees.size());

        reconstruction += (lambda * reconstructed_error) / ATA;

        // Make sure no voxel is negative as it is not possible
        int npixel = reconstruction.size();
#pragma omp parallel for
        for (int pixel_id = 0; pixel_id < npixel; ++pixel_id)
        {
            reconstruction.m_p_voxel_set[pixel_id] = max(0.0f, reconstruction.m_p_voxel_set[pixel_id]); // Avoid using branching
        }

#ifdef DEBUG

        stringstream filename1;
        stringstream filename2;
        stringstream filename3;
        stringstream filename4;

        filename1.fill('0');
        filename2.fill('0');
        filename3.fill('0');
        filename4.fill('0');

        filename1 << "estimated_projections_" << std::setw(3) << i << ".tif";
        filename2 << "error_iteration_" << std::setw(3) << i << ".tif";
        filename3 << "reconstructed_error_iteration_" << std::setw(3) << i << ".tif";
        filename4 << "art_iteration_" << std::setw(3) << i << ".tif";

        estimated_projections.saveTIFF(filename1.str());
        error.saveTIFF(filename2.str());
        reconstructed_error.saveTIFF(filename3.str());
        reconstruction.saveTIFF(filename4.str());
#endif
    }

    return reconstruction;
}





//-----------------------------------------
void Image::loadTIFF(const char* aFileName)
//-----------------------------------------
{
    // Empty the image
    clear();

    // Open the TIFF file
    TIFF* tif = TIFFOpen(aFileName, "r");

    // The TIFF file is open
    if (tif)
    {
	    uint32 width = 0;
	    uint32 height = 0;

        // Read dimensions of image
        if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) != 1)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Failed to read width" << endl;

            throw error_message.str();
        }

        if (TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height) != 1)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Failed to read height" << endl;

            throw error_message.str();
        }

        // Get the number of bits per pixel
        TIFF_UINT16_T number_of_bits_per_pixel = 0;
        if (TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &number_of_bits_per_pixel) != 1)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Failed to read number of bits per pixel" << endl;

            throw error_message.str();
        }

        // Get the number of samples per pixel
        TIFF_UINT16_T number_of_samples_per_pixel = 0;
        if (TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &number_of_samples_per_pixel) != 1)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Failed to read orientation" << endl;

            throw error_message.str();
        }

        // Get the format of a sample
        TIFF_UINT16_T format_of_a_sample = 0;
        if (TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &format_of_a_sample) != 1)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Failed to read data type" << endl;

            throw error_message.str();
        }

        if (number_of_bits_per_pixel != 32)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Not a valid TIFF format for this application" << endl;

            throw error_message.str();
        }

        if (number_of_samples_per_pixel != 1)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Not a valid TIFF format for this application" << endl;

            throw error_message.str();
        }

        if (format_of_a_sample != SAMPLEFORMAT_IEEEFP)
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Not a valid TIFF format for this application" << endl;

            throw error_message.str();
        }

	    uint32 npixels = width * height;
        tsize_t scanline_size = TIFFScanlineSize(tif);
        uint32*	raster = (uint32*) _TIFFmalloc(height * scanline_size);

        m_p_voxel_set = vector<float>(npixels);
        if (scanline_size / width == 4 && (scanline_size % width) == 0)
        {
            for (uint32 y = 0; y < height; ++y)
            {
                TIFFReadScanline(tif, &m_p_voxel_set[y * width], y);
            }
        }
        else
        {
            stringstream error_message;
            error_message << "EXCEPTION caught: " << endl <<
                "File: " << __FILE__ << endl <<
                "Function: " << __FUNCTION__ << endl <<
                "Line: " << __LINE__ << endl <<
                "Message: " << "Not a valid TIFF format for this application" << endl;

            throw error_message.str();
        }

        _TIFFfree(raster);
        TIFFClose(tif);

        // Save the image size
	    m_cols = width;
	    m_rows = height;
        m_slices = 1;
    }
    // The TIFF file is not open
    else
    {
        stringstream error_message;
        error_message << "EXCEPTION caught: " << endl <<
            "File: " << __FILE__ << endl <<
            "Function: " << __FUNCTION__ << endl <<
            "Line: " << __LINE__ << endl <<
            "Message: " << "Failed to open image (" << aFileName << ")" << endl;

        throw error_message.str();
    }
}


//-----------------------------------------------
void Image::saveTIFF(const char* aFileName) const
//-----------------------------------------------
{
    TIFF* tif = TIFFOpen(aFileName, "w");

    // The file is open
    if (tif)
    {
        unsigned int number_of_pixels = m_cols * m_rows;
        TIFF_UINT16_T number_of_bits_per_pixel = 32;

        // Set the image width
        TIFF_UINT32_T image_width = m_cols;


        // Write all the slices
        for (unsigned int k(0); k < m_slices; k++)
        {
            // Create a new directory
            if (m_slices > 1)
            {
                if (TIFFWriteDirectory(tif) != 1)
                {
                    cerr << "Could not write the TIFF file" << endl;
                    exit(1);
                }
            }

            if (TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, image_width) != 1)
            {
                cerr << "Could not set the image width in the TIFF file" << endl;
                exit(1);
            }

            // Set the image height
            TIFF_UINT32_T image_height = m_rows;
            if (TIFFSetField(tif, TIFFTAG_IMAGELENGTH, image_height) != 1)
            {
                cerr << "Could not set the image height in the TIFF file" << endl;
                exit(1);
            }

            // Set the number of bits per pixel
            if (TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, number_of_bits_per_pixel) != 1)
            {
                cerr << "Could not set the number of bits per pixel in the TIFF file" << endl;
                exit(1);
            }

            // Set the compression
            TIFF_UINT16_T compression = COMPRESSION_LZW;
            if (TIFFSetField(tif, TIFFTAG_COMPRESSION, compression) != 1)
            {
                cerr << "Could not set the compression in the TIFF file" << endl;
                exit(1);
            }

            // Set the photometric interpretation
            TIFF_UINT16_T photometric_interpretation = PHOTOMETRIC_MINISBLACK;
            if (TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, photometric_interpretation) != 1)
            {
                cerr << "Could not set the photometric interpretation in the TIFF file" << endl;
                exit(1);
            }

            // Set the orientation
            TIFF_UINT16_T orientation = ORIENTATION_TOPLEFT;
            if (TIFFSetField(tif, TIFFTAG_ORIENTATION, orientation) != 1)
            {
                cerr << "Could not set the orientation in the TIFF file" << endl;
                exit(1);
            }

            // Set the number of samples per pixel
            TIFF_UINT16_T number_of_samples_per_pixel = 1;
            if (TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, number_of_samples_per_pixel) != 1)
            {
                cerr << "Could not set the orientation in the TIFF file" << endl;
                exit(1);
            }
    /*
            // Set the min sample value
            if (TIFFSetField(tif, TIFFTAG_MINSAMPLEVALUE, min_value) != 1)
            {
                cerr << "Could not set the min sample value in the TIFF file" << endl;
                exit(1);
            }

            // Set the max sample value
            if (TIFFSetField(tif, TIFFTAG_MAXSAMPLEVALUE, max_value) != 1)
            {
                cerr << "Could not set the max sample value in the TIFF file" << endl;
                exit(1);
            }*/

            // Set the planar configuration
            TIFF_UINT16_T planar_configuration = PLANARCONFIG_CONTIG;
            if (TIFFSetField(tif, TIFFTAG_PLANARCONFIG, planar_configuration) != 1)
            {
                cerr << "Could not set the planar configuration in the TIFF file" << endl;
                exit(1);
            }

            // Set the format of a sample
            TIFF_UINT16_T format_of_a_sample = SAMPLEFORMAT_IEEEFP;
            if (TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, format_of_a_sample) != 1)
            {
                cerr << "Could not set the data type in the TIFF file" << endl;
                exit(1);
            }

            // Length in memory of one row of pixel in the image
            tsize_t number_of_bytes_per_line = number_of_samples_per_pixel * m_cols;

            // We set the strip size of the file to be size of one row of pixels
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, number_of_bytes_per_line));

            // Write the image data, line by line
            int error_code = 0;
            for (unsigned int j = 0; j < m_rows; ++j)
            {
                if (TIFFWriteScanline(tif, const_cast<float*>(&m_p_voxel_set[0]) + m_cols * j, j, 0) != 1)
                {
                    // Close the file
                    TIFFClose(tif);

                    // Generate an error
                    cerr << "Could not write the TIFF file" << endl;
                    exit(1);
                }
            }

            // Flush the data in the TIFF file
            if (m_slices > 1)
            {
                if (TIFFFlush(tif) != 1)
                {
                    cerr << "Could not write the TIFF file" << endl;
                    exit(1);
                }
            }
        }

        TIFFClose(tif);
    }
    // The TIFF file is not open
    else
    {
        cerr << "Failed to open image" << endl;
        exit(EXIT_FAILURE);
    }
}


//-----------------------------------------------
void Image::saveTEXT(const char* aFileName) const
//-----------------------------------------------
{
    ofstream output(aFileName);

    if (output.is_open())
    {
        for (int y = 0; y < m_rows; ++y)
        {
            for (int x = 0; x < m_cols; ++x)
            {
                output << m_p_voxel_set[x + y * m_cols];

                if (x < m_cols - 1) output << " ";
            }

            if (y < m_rows - 1) output << endl;
        }
    }
    else
    {
        cerr << "Failed to write image in " << aFileName << endl;
        exit(EXIT_FAILURE);
    }
}


//-----------------------------------------------
Image Image::rotate(double anAngleInDegrees) const
//-----------------------------------------------
{
    Image temp(m_cols, m_rows, m_slices);

    double angle_in_radian(anAngleInDegrees * M_PI / 180);

    double sinf = std::sin(angle_in_radian);
    double cosf = std::cos(angle_in_radian);

    double centre_x = (m_cols - 1.0) / 2.0;
    double centre_y = (m_rows - 1.0) / 2.0;

#pragma omp parallel for collapse(3)
    for (int k = 0; k < m_slices; ++k)
    {
        for (unsigned j = 0; j < m_rows; ++j)
        {
            for (unsigned i = 0; i < m_cols; ++i)
            {
                double temp_x(i - centre_x);
                double temp_y(j - centre_y);

                // temp_x -= centre_x;
                // temp_y -= centre_y;

                // double x = temp_x * cosf - temp_y * sinf + centre_x;
                // double y = temp_x * sinf + temp_y * cosf + centre_y;

                double x = temp_x * cosf + temp_y * sinf + centre_x;
                double y = temp_y * cosf - temp_x * sinf + centre_y;


                // x += 0.5;
                // y += 0.5;


                // x += (m_cols - 1.0) / 2.0;
                // y += (m_rows - 1.0) / 2.0;
                //
                int x1 = floor(x);
                int y1 = floor(y);
                //
                int x2 = x1 + 1;
                int y2 = y1 + 1;

                if (x1 >= 0 && x1 < m_cols &&
                    y1 >= 0 && y1 < m_rows)
                {
                    temp(i, j, k) = getPixel(x1, y1, k);

                    if (x2 >= 0 && x2 < m_cols &&
                            y2 >= 0 && y2 < m_rows)
                    {
                        double c11 = getPixel(x1, y1, k);
                        double c12 = getPixel(x1, y2, k);
                        double c22 = getPixel(x2, y2, k);
                        double c21 = getPixel(x2, y1, k);

                        double cx1 = c11 + (c21 - c11) * (x - double(x1)) / double(x2 - x1);
                        double cx2 = c12 + (c22 - c12) * (x - double(x1)) / double(x2 - x1);

                        double cy  = cx1 + (cx2 - cx1) * double(y - y1) / double(y2 - y1);

                        temp(i, j, k) = cy;
                    }
                }
            }
        }
    }

    return (temp);
}



//-----------------------------------------------------------------------------------
Image Image::setCanvasSize(int cols, int rows, int slices, float aDefaultValue) const
//-----------------------------------------------------------------------------------
{
    Image temp(cols, rows, slices, aDefaultValue);

    int x_offset((int(cols)  - int(m_cols))  / 2);
    int y_offset((int(rows) - int(m_rows)) / 2);
    int z_offset((int(slices) - int(m_slices)) / 2);

	for (int k = 0; k < m_slices; ++k)
    {
        for (int j = 0; j < m_rows; ++j)
        {
            for (int i = 0; i < m_cols; ++i)
            {
                if (i + x_offset >= 0 && i + x_offset < cols &&
                    j + y_offset >= 0 && j + y_offset < rows &&
                    k + z_offset >= 0 && k + z_offset < slices)
                {
                    temp(i + x_offset, j + y_offset, k + z_offset) = getPixel(i, j, k);
                }
            }
        }
    }

    return (temp);
}
