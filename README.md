# 3D Geometric Nonlinear Diffusion Filter

This [paper](https://arafathm.github.io/assets/pdf/mah2015b.pdf) presented two near real-time nonlinear anisotropic diffusion filtering (NADF) methods for the 2D and 3D X-ray computed tomography (CT) and magnetic resonance (MR) image denoising. Typically, NADFs are preferred for the medical image denoising due to its edge preserving feature though they are computationally expensive. Recently, a computation-time efficient 2D NADF has been proposed which uses local pixel intensity-based geometric parameters for diffusion. But it has limitations resulting from (i) its assumption that the neighboring pixels are non-noisy while deciding on an interrogated pixel being noisy or not, and (ii) its confinement of working only on a 2D image. Motivated from this, we propose an improved 2D NADF method that uses additional neighboring pixels in an effective way to lower the noise impact on the estimated geometric parameters. We also extend our 2D method into 3D that considers all the three directions for information diffusion. The performance of the proposed methods is evaluated using a 3D synthetic phantom, and in vivo CT and MR data.

Please note that we shared a prelminary version of the proposed 3D approach (i.e. gnldf3D.m function), where only North-South, East-West, and Up-down pixels are used during diffusion. This code can be easily modified for the proposed 2D and 3D approaches. 

## Citations
If you find this work useful, please cite the first or all of the following papers:
```
@inproceedings{hussain2015towards,
  title={Towards real-time 3D geometric nonlinear diffusion filter and its application to CT and MR imaging},
  author={Hussain, Mohammad Arafat and Shourov, Riad Mashrub and Khan, Shamima Nasrin},
  booktitle={2015 18th International Conference on Computer and Information Technology (ICCIT)},
  pages={462--467},
  year={2015},
  organization={IEEE}
}
```
