#include <vector>
#include <string>

#include "Image.h"
#include "cmdline.h"

using namespace std;

int main(int argc, char** argv)
{
    int error_code = EXIT_SUCCESS;

    try
    {
        struct gengetopt_args_info args_info;
        if (cmdline_parser(argc, argv, &args_info) != 0)
        {
            throw "Wrong command line arguments";
        }

        // Optional argument
        Image reference_CT;
        if (args_info.reference_CT_given == 1)
        {
            cout << "Load reference CT volume:\t" << args_info.reference_CT_arg << endl << endl;
            reference_CT.loadTIFF(args_info.reference_CT_arg);
        }

        cout << "Load reference sinogram:\t" << args_info.reference_sinogram_arg << endl;
        Image reference_sinogram(args_info.reference_sinogram_arg);
        cout << "Number of colums in the sinogram:\t" << reference_sinogram.getCols() << endl;
        cout << "Number of angles in the sinogram:\t" << reference_sinogram.getRows() << endl;
        cout << "Number of slices in the sinogram:\t" << reference_sinogram.getSlices() << endl;

        vector<double> angle_set;
        double angle_step = 180.0 / reference_sinogram.getRows();
        for (int i = 0; i < reference_sinogram.getRows(); ++i)
        {
            angle_set.push_back(i * angle_step);
        }
        cout << "First angle (in degrees):\t" << angle_set[0] << endl;
        cout << "Last angle (in degrees):\t" << angle_set[angle_set.size() - 1] << endl;

        if (angle_set.size() >= 2)
        {
            cout << "Angular step (in degrees):\t" << angle_set[1] << endl;
        }
        cout << endl;

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

        Image my_reconstruction_art = my_sinogram.iradonART(angle_set, 1.0, args_info.iterations_arg);
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

        cout << endl << "Save reconstructed volume:\t" << args_info.reconstruction_arg << endl;

        // The reference CT volume is loaded
        if (reference_CT.size())
        {
            cout << endl << "Compare the reference CT with the reconstructed CT:" << endl;
            cout << "\tRMSE(ref, SBP): " << reference_CT.rmse(my_reconstruction_sbp) << endl;
            cout << "\tZNCC(ref, SBP): " << 100.0 * reference_CT.zncc(my_reconstruction_sbp) << "%" << endl << endl;
            cout << "\tRMSE(ref, ART): " << reference_CT.rmse(my_reconstruction_art) << endl;
            cout << "\tZNCC(ref, ART): " << 100.0 * reference_CT.zncc(my_reconstruction_art) << "%" << endl;
        }
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
