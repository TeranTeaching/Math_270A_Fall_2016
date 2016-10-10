This is a test.

Copyright (c) 2016 Theodore Gast, Chuyuan Fu, Chenfanfu Jiang, Joseph Teran

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

If the code is used in an article, the following paper shall be cited:
@techreport{qrsvd:2016,
  title={Implicit-shifted Symmetric QR Singular Value Decomposition of 3x3 Matrices},
  author={Gast, Theodore and Fu, Chuyuan and Jiang, Chenfanfu and Teran, Joseph},
  year={2016},
  institution={University of California Los Angeles}
}

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

################################################################################
ImplicitQRSVD.h implements 2D and 3D polar decompositions and SVDs.
Tools.h provides a random number generator and a timer.
The code is tested with g++ 5.3. It uses std=c++14. It relies on Eigen3 and is header-only.
################################################################################
To run the benchmark:
    make;
    ./svd
################################################################################
To use the SVD code: (T may be float or double)
2D Polar:
    Eigen::Matrix<T, 2, 2> A,R,S;
    A<<1,2,3,4;
    JIXIE::polarDecomposition(A, R, S);
    // R will be the closest rotation to A
    // S will be symmetric
2D SVD:
    Eigen::Matrix<T, 2, 2> A;
    A<<1,2,3,4;
    Eigen::Matrix<T, 2, 1> S;
    Eigen::Matrix<T, 2, 2> U;
    Eigen::Matrix<T, 2, 2> V;
    JIXIE::singularValueDecomposition(A,U,S,V);
    // A = U S V'
    // U and V will be rotations
    // S will be singular values sorted by decreasing magnitude. Only the last one may be negative.
3D Polar:
    Eigen::Matrix<T, 3, 3> A,R,S;
    A<<1,2,3,4,5,6;
    JIXIE::polarDecomposition(A, R, S);
    // R will be the closest rotation to A
    // S will be symmetric
3D SVD:
    Eigen::Matrix<T, 3, 3> A;
    A<<1,2,3,4,5,6;
    Eigen::Matrix<T, 3, 1> S;
    Eigen::Matrix<T, 3, 3> U;
    Eigen::Matrix<T, 3, 3> V;
    JIXIE::singularValueDecomposition(A,U,S,V);
    // A = U S V'
    // U and V will be rotations
    // S will be singular values sorted by decreasing magnitude. Only the last one may be negative.
################################################################################
