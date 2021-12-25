# -*- coding: utf-8 -*-

def B_Quad8(xi, eta, X):
    import numpy as np
    D = np.matrix(np.zeros([4, 16]))  # Initialize D with zeros
    # Derivatives of shape functions wrt xi & eta
    dNdxi = 1/4.*np.matrix([[eta + 2. * xi * (1 - eta) - eta ** 2,
                             -eta + 2. * xi * (1 - eta) + eta ** 2,
                             eta + eta ** 2 + 2. * xi * (1 + eta),
                             -eta + 2. * xi * (1 + eta) - eta ** 2,
                             4*(-xi * (1 - eta)),
                             2. - 2. * eta ** 2,
                             4*(-xi * (1 + eta)),
                             -2. + 2. * eta ** 2],
                            [xi - xi ** 2 + (2. - 2. * xi) * eta,
                             -xi - xi ** 2 + (2. + 2. * xi) * eta,
                             xi + (2. + 2. * xi) * eta + xi ** 2,
                             -xi + xi ** 2 + (2. - 2. * xi) * eta,
                             -2. + 2. * xi ** 2,
                             -4 * (1 + xi) * eta,
                             2. - 2. * xi ** 2,
                             -4 * (1 - xi) * eta]])

    for i in range(2):
        for j in range(8):
            D[i, 2 * j] = dNdxi[i, j]
            D[i + 2, 2 * j + 1] = dNdxi[i, j]
              # Arrange shape function derivatives into D
    J = dNdxi * X
    detJ = np.linalg.det(J)  # Determinant of Jacobian J
    invJ = np.linalg.inv(J)  # Inverse of Jacobian
    dxidX = np.matrix([[invJ[0, 0], invJ[0, 1], 0, 0],
                       [0, 0, invJ[1, 0], invJ[1, 1]],
                       [invJ[1, 0], invJ[1, 1], invJ[0, 0], invJ[0, 1]]])

    B = dxidX * D  # Shape function derivatives wrt x and y: Eq.(2.39)
    return B, detJ


def Quad8R_Stiffness(X, Cmat, t):
    import numpy as np
    Tangent = np.matrix(np.zeros([16, 16]))  # Initialize tangent with zeros
    Gauss_pos = 1./np.sqrt(3.0)  # Gauss point location

    for jGauss in range(2):    # 2 by 2 Gauss integration loops
        eta = (-1) ** (jGauss + 1) * Gauss_pos  # Natural coordinate eta
        for iGauss in range(2):
            xi = ((-1) ** (iGauss + 1)) * Gauss_pos   # Natural coordinate xi
            [B, detJ] = B_Quad8(xi, eta, X)     # B for Eq.(4.83)
            Tangent = Tangent + B.T * Cmat * B * detJ * t
    return Tangent

def Quad8_Stiffness(X, Cmat, t):
    import numpy as np
    Tangent = np.matrix(np.zeros([16, 16]))  # Initialize tangent with zeros
    W  = [5.0/9.0, 8.0/9.0, 5.0/9.0]
    GP = [-np.sqrt(0.6), 0.0, np.sqrt(0.6)]
    for jGauss in range(3):   # 3 by 3 Gauss integration loops
        eta = GP[jGauss]      # Natural coordinate eta
        for iGauss in range(3):
            xi = GP[iGauss]   # Natural coordinate xi
            [B, detJ] = B_Quad8(xi, eta, X)     # B for Eq.(4.83)
            Tangent = Tangent + B.T * Cmat * B * detJ * t * W[jGauss]*W[iGauss]
    return Tangent

def B_Quad4(xi, eta, X):
    import numpy as np
    G = np.matrix(np.zeros([4, 8]))  # Initialize D with zeros
    # Derivatives of shape functions wrt xi & eta
    dNdxi = 1/4.*np.matrix([[(eta - 1.0),
                             (1.0 - eta),
                             (1.0+eta),
                             (-1.0-eta)],
                            [(xi - 1.0),
                             (-1.0-xi),
                             (1.0+xi),
                             (1.0-xi)]])

    for i in range(2):
        for j in range(4):
            G[i, 2 * j] = dNdxi[i, j]
            G[i + 2, 2 * j + 1] = dNdxi[i, j]
              # Arrange shape function derivatives into D
    J = dNdxi * X
    detJ = np.linalg.det(J)  # Determinant of Jacobian J
    invJ = np.linalg.inv(J)  # Inverse of Jacobian
    A = np.matrix([[invJ[0, 0], invJ[0, 1], 0, 0],
                   [0, 0, invJ[1, 0], invJ[1, 1]],
                   [invJ[1, 0], invJ[1, 1], invJ[0, 0], invJ[0, 1]]])

    B = A * G  # Shape function derivatives wrt x and y: Eq.(2.39)
    return B, detJ

def Quad4_Stiffness(X, Cmat, t):
    import numpy as np
    Tangent = np.matrix(np.zeros([8, 8]))  # Initialize tangent with zeros
    Gauss_pos = 1./np.sqrt(3.0)  # Gauss point location

    for jGauss in range(2):    # 2 by 2 Gauss integration loops
        eta = (-1) ** (jGauss + 1) * Gauss_pos  # Natural coordinate eta
        for iGauss in range(2):
            xi = ((-1) ** (iGauss + 1)) * Gauss_pos   # Natural coordinate xi
            [B, detJ] = B_Quad4(xi, eta, X)     # B for Eq.(4.83)
            Tangent = Tangent + B.T * Cmat * B * detJ * t
    return Tangent