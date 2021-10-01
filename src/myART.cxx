#include <vector>

#include "Image.h"

using namespace std;

int main(int argc, char** argv)
{
    int error_code = EXIT_SUCCESS;

    try
    {
        //Image reference_CT("../reference_CT.tif");

        Image reference_CT("../data/shepp_logan/reference_CT.tif");
        Image reference_sinogram("../data/shepp_logan/sinogram.tif");

        vector<double> angle_set;
        double angle_step = 180.0 / reference_sinogram.getRows();
        for (int i = 0; i < reference_sinogram.getRows(); ++i)
        {
            angle_set.push_back(i * angle_step);
        }

        Image my_sinogram = reference_CT.radon(angle_set);
        my_sinogram.saveTIFF("my_sinogram.tif");

        Image my_reconstruction_sbp = my_sinogram.iradon(angle_set);
        my_reconstruction_sbp.saveTIFF("my_reconstruction.tif");

        int x = (my_reconstruction_sbp.getCols() - reference_CT.getCols()) / 2;
        int y = (my_reconstruction_sbp.getRows() - reference_CT.getRows()) / 2;
        int z = (my_reconstruction_sbp.getSlices() - reference_CT.getSlices()) / 2;

        my_reconstruction_sbp = my_reconstruction_sbp.getROI(x, y, z,
            reference_CT.getCols(),
            reference_CT.getRows(),
            reference_CT.getSlices());

        my_reconstruction_sbp.saveTIFF("my_cropped_reconstruction_sbp.tif");

        Image my_reconstruction_art = my_sinogram.iradonART(angle_set, 1.0, 200);
        my_reconstruction_art.saveTIFF("my_reconstruction_ART.tif");

        x = (my_reconstruction_art.getCols() - reference_CT.getCols()) / 2;
        y = (my_reconstruction_art.getRows() - reference_CT.getRows()) / 2;
        z = (my_reconstruction_art.getSlices() - reference_CT.getSlices()) / 2;

        my_reconstruction_art = my_reconstruction_art.getROI(x, y, z,
            reference_CT.getCols(),
            reference_CT.getRows(),
            reference_CT.getSlices());

        my_reconstruction_art.saveTIFF("my_cropped_reconstruction_art.tif");


        x = (my_sinogram.getCols() - reference_sinogram.getCols()) / 2;
        y = (my_sinogram.getRows() - reference_sinogram.getRows()) / 2;
        z = (my_sinogram.getSlices() - reference_sinogram.getSlices()) / 2;

        my_sinogram = my_sinogram.getROI(x, y, z,
            reference_sinogram.getCols(),
            reference_sinogram.getRows(),
            reference_sinogram.getSlices());

        my_sinogram.saveTIFF("my_cropped_sinogram.tif");

        cout << "RMSE(ref, SBP): " << reference_CT.rmse(my_reconstruction_sbp) << endl;
        cout << "ZNCC(ref, SBP): " << 100.0 * reference_CT.zncc(my_reconstruction_sbp) << "%" << endl;
        cout << "RMSE(ref, ART): " << reference_CT.rmse(my_reconstruction_art) << endl;
        cout << "ZNCC(ref, ART): " << 100.0 * reference_CT.zncc(my_reconstruction_art) << "%" << endl;
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << endl;
        error_code = EXIT_FAILURE;
    }
    catch (const std::string& e)
    {
        cerr << e << endl;
        error_code = EXIT_FAILURE;
    }
    catch (const char* e)
    {
        cerr << e << endl;
        error_code = EXIT_FAILURE;
    }
    catch (...)
    {
        cerr << "EXCEPTION caught: unknown error" << endl;
        error_code = EXIT_FAILURE;
    }

    return error_code;
}
