# [sba](https://users.ics.forth.gr/~lourakis/sba/) - *Generic Sparse Bundle Adjustment Based on the Levenberg-Marquardt Algorithm*

This is *sba*, a library for generic sparse [bundle adjustment](http://en.wikipedia.org/wiki/Bundle_adjustment) in 3D reconstruction that is distributed under the GPLv2.
Bundle adjustment simultaneously refines an initial motion and structure estimate by minimizing the reprojection error between the observed and predicted image points.
Solving this problem with a general-purpose implementation of a non-linear least squares algorithm (e.g.,  [Levenberg-Marquardt](http://en.wikipedia.org/wiki/Levenberg-Marquardt_algorithm)) incurs high computational costs due to the large number of unknowns involved.

sba exploits the sparse block structure of the underlying [normal equations](https://en.wikipedia.org/wiki/Ordinary_least_squares#Normal_equations) with a tailored sparse variant of the LM algorithm that leads to considerable computational gains.
sba is generic in the sense that it grants the user full control over the definition of the parameters describing cameras and 3D structure. For example, see this application to [SfM for light field images](https://openaccess.thecvf.com/content_CVPR_2019/html/Nousias_Large-Scale_Metric_Structure_From_Motion_for_Unordered_Light_Fields_CVPR_2019_paper.html).

More details can be found in the corresponding [publication](http://www.ics.forth.gr/~lourakis/sba/sba-toms.pdf) and [sba's web site](https://users.ics.forth.gr/~lourakis/sba/).

## Required libraries
sba requires [LAPACK](https://github.com/Reference-LAPACK/lapack) or an equivalent numerical linear algebra library.

## Build
After installing the required dependencies, create a ``build`` directory in the root of the cloned repository and run ``cmake`` .

## Matlab interface
sba has a matlab mex interface in the ``matlab`` subdirectory. See the included ``README.txt`` for more information.

## Cite as
If you use this code in your published work, please cite the following [paper](http://www.ics.forth.gr/~lourakis/sba/sba-toms.pdf):<br><br>
<pre>
@article{Lourakis09sba,
    author={M. I. A. Lourakis and A. A. Argyros},
    title={{SBA}: A Software Package for Generic Sparse Bundle Adjustment},
    journal={ACM Transactions on Mathematical Software},
    volume={36},
    number={1},
    year={2009},
    pages={1-30},
    publisher={ACM},
    address={New York, NY, USA},
    doi={http://doi.acm.org/10.1145/1486525.1486527}
}
<pre>
