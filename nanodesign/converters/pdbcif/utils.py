
# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# 2021.02.25.: moved matrix math to seperate file @Elija, Feigl

"""
This module is contains various functions.
"""
import numpy as np
from math import sqrt, cos, acos, sin, pi


def _Rx(deg):
    """ Generate a rotation matrix for an angle about the X-axis. """
    rad = pi*(deg / 180.0)
    c = cos(rad)
    s = sin(rad)
    R = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    return R


def _Ry(deg):
    """ Generate a rotation matrix for an angle about the Y-axis. """
    rad = pi*(deg / 180.0)
    c = cos(rad)
    s = sin(rad)
    R = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    return R


def _vrrotmat2vec(R):
    """ Extract the equivalent rotation about an axis from a rotation matrix.
    """
    # print("-------------- _vrrotmat2vec -------------- ")
    r00 = R[0, 0]
    r01 = R[0, 1]
    r02 = R[0, 2]
    r10 = R[1, 0]
    r11 = R[1, 1]
    r12 = R[1, 2]
    r20 = R[2, 0]
    r21 = R[2, 1]
    r22 = R[2, 2]
    # print("R %s " % str(R))
    # print("r21-r12 %s " % str(r21-r12))
    # print("r02-r20 %s " % str(r02-r20))
    # print("r10-r01 %s " % str(r10-r01))

    angle = acos((r00 + r11 + r22 - 1)/2.0)
    # print("angle %g " % angle)

    # Check for singularities at angle=0.
    if abs(angle) < 0.0001:
        angle = 0.0
        x = 1.0
        y = 0.0
        z = 0.0

    # Check for singularity at angle=180. Calculate axis using largest diagonal
    # term.
    elif abs(angle - pi) < 0.0001:
        epsilon = 0.01
        xx = (r00+1) / 2.0
        yy = (r11+1) / 2.0
        zz = (r22+1) / 2.0
        xy = (r01+r10) / 4.0
        xz = (r02+r20) / 4.0
        yz = (r12+r21) / 4.0

        if (xx > yy) and (xx > zz):
            if xx < epsilon:
                x = 0.0
                y = 0.7071
                z = 0.7071
            else:
                x = sqrt(xx)
                y = xy / x
                z = xz / x
        elif yy > zz:
            if yy < epsilon:
                x = 0.7071
                y = 0.0
                z = 0.7071
            else:
                y = sqrt(yy)
                x = xy / y
                z = yz / y
        else:
            if zz < epsilon:
                x = 0.7071
                y = 0.7071
                z = 0.0
            else:
                z = sqrt(zz)
                x = xz / z
                y = yz / z

    # Calculate axis for non-singular angle.
    else:
        x = (r21 - r12) / sqrt(pow(r21-r12, 2) +
                               pow(r02-r20, 2) + pow(r10-r01, 2))
        y = (r02 - r20) / sqrt(pow(r21-r12, 2) +
                               pow(r02-r20, 2) + pow(r10-r01, 2))
        z = (r10 - r01) / sqrt(pow(r21-r12, 2) +
                               pow(r02-r20, 2) + pow(r10-r01, 2))

    # Check the result.
    check_result = False
    check_result = True
    if check_result:
        c = cos(angle)
        s = sin(angle)
        t = 1.0 - c
        M = np.array([[0, 0, -1], [-1, 0, 0], [0, 1, 0]], dtype=float)
        M[0, 0] = t*x*x + c
        M[0, 1] = t*x*y - z*s
        M[0, 2] = t*x*z + y*s
        M[1, 0] = t*x*y + z*s
        M[1, 1] = t*y*y + c
        M[1, 2] = t*y*z - x*s
        M[2, 0] = t*x*z - y*s
        M[2, 1] = t*y*z + x*s
        M[2, 2] = t*z*z + c

        for i in range(0, 3):
            for j in range(0, 3):
                d = R[i, j] - M[i, j]
                if abs(d) > 0.0001:
                    print(
                        "[_vrrotmat2vec] **** WARNING: R-M not zero: %g " % (
                            R[i, j] - M[i, j])
                    )

    return np.array([x, y, z], dtype=float), angle


def _vrrotvec2mat(axis, theta):
    """ Create a rotation matrix to rotate theta degrees about the axis defined
    by vec. """
    s = np.sin(theta)
    c = np.cos(theta)
    t = 1 - c
    # print("[vrrotvec2mat] theta=%s" % str(theta))

    # normalize the vector
    x, y, z = axis / np.linalg.norm(axis)
    # print("[vrrotvec2mat] x=%s" % str(x))
    # print("[vrrotvec2mat] y=%s" % str(y))
    # print("[vrrotvec2mat] z=%s" % str(z))

    return np.array([[t*x*x + c,   t*x*y - s*z,  t*x*z + s*y],
                     [t*x*y + s*z, t*y*y + c,    t*y*z - s*x],
                     [t*x*z - s*y, t*y*z + s*x,  t*z*z + c]])
