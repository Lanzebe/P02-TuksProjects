# -*- coding: utf-8 -*-
"""
MKM321 Two-Dimensional Linear Finite Element Program
@author: Dominic Kafka
@contributor: Daniel N. Wilke

Adapted by Dominic Kafka from the postgraduate MEE780
and MEE781 finite element codes developed in Matlab.
"""
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import time

# Import read_input_file from module preprocess
from preprocess import read_input_file

# Import stiffness matrix functions for elements
# Each function depends on (xi, eta, X)
# where xi and eta are the local element coordinates betweem -1 and 1
# where X is the global nodal coordinates for each node
from elements import B_Quad8
from elements import Quad8_Stiffness
from elements import Quad8R_Stiffness
from elements import B_Quad4
from elements import Quad4_Stiffness
from elements import Quad4R_Stiffness

# Import a number of functions to allow us to visualize the results
from postprocess import nodal_values
from postprocess import calc_von_mises
from postprocess import calc_tresca
from postprocess import write_output_file
from postprocess import PlotIso
from postprocess import PlotData
from postprocess import PlotGeom
from postprocess import PlotClearAll
from postprocess import SmoothData
from postprocess import NewTri

def launch_fem(case, MagFac=1.0, IntegrationType = 'FI'):
    print(' ')
    print(' ')
    print(' ')
    print('-----------------------------------------------')
    print('MKM321: Two-Dimensional Linear Finite Element Program')
    print('-----------------------------------------------')

    file_out = case + '.out'
    filename = case + '.inp'
    tic = time.time()

    DD = read_input_file(filename) # DD is the data dictionary

    if len(DD['ELEMENT'][0]) <= 5:
        ElType = 'q4'
        NodesPerEl = 4
    elif len(DD['ELEMENT'][0]) > 4:
        ElType = 'q8'
        NodesPerEl = 8

    #Element types 1, 3 and 5 are plane strain: plane = 0
    #Element types 2, 4 and 6 are plain stress: plane = 1
    plane = 1
    if np.mod(DD['ELTYPE'][0][0], 2):
        plane = 0

    DofPerEl = 2 * NodesPerEl
    TotKvec = DofPerEl ** 2

    toc = time.time()
    finish = toc-tic
    print('-----------------------------------------------')
    print('Done processing input file :','{:8.4f}'.format(finish), 'seconds.')
    nodes = np.array(DD['NODE'])[:,0].tolist()
    coor = np.array(DD['NODE'])[:, 1:]
    elnodes = np.array(DD['ELEMENT'],int)
    nnodes = np.shape(DD['NODE'])[0]
    nelem  = np.shape(DD['ELEMENT'])[0]
    ndispl = np.shape(DD['BOUNDARY'])[0]
    ncload = np.shape(DD['CLOAD'])[0]
    displ = sorted(DD['BOUNDARY'])

    # Find prescribed (pdof) and free (fdof) degrees of freedom
    dof = np.ones([nnodes, 2],int)
    Up = []
    for node, dimension, displacement in displ:
        dof[nodes.index(int(node)), int(dimension) - 1] = -1
        Up.append(displacement)
    Up = np.matrix(Up).T
    dof = dof.flatten()
    pdof = np.flatnonzero(dof == -1)
    fdof = np.flatnonzero(dof != -1)

    Allpg = np.zeros(len(dof),int)
    Pos = 0
    Neg = 0
    for i in range(len(dof)):
        if dof[i]==1:
            Pos = Pos + 1
            Allpg[i]=Pos
        else:
            Neg = Neg-1
            Allpg[i]=Neg

    # Construct elasticity tensor D
    # If plane strain
    elas = DD['EMODULUS'][0][0]; pois = DD['POISSON'][0][0]
    if plane == 0:
        e = elas / float(1 - pois ** 2)
        nu = pois / float(1 - pois)
    else:
        e = elas
        nu = pois

    c = e / float(1 - nu ** 2)
    matD = np.asarray([[c,c*nu,0],[c*nu,c,0],[0,0,c*(1-nu)/2]])
    ndnum = range(1, (1 + NodesPerEl))
    [colpos, rowpos] = np.meshgrid(range(DofPerEl), range(DofPerEl))

    colpos = colpos.flatten()
    rowpos = rowpos.flatten()

    F = np.zeros([2 * nnodes, 1])

    stress = np.asarray(np.zeros([nelem, 12]))
    strain = np.asarray(np.zeros([nelem, 12]))

    # Prescribed displacement applied
    U = np.asarray(np.zeros([2*nnodes,1]))
    U[pdof, :] = Up
    tic = time.time()

    #Main loop over elements. Compute k_elem and assemble
    #Initialize global stiffness matrix vectors;
    row_vec = np.array(np.zeros([TotKvec * nelem]),int)
    col_vec = np.array(np.zeros([TotKvec * nelem]),int)
    stiff_vec = np.array(np.zeros([TotKvec * nelem]))


    for i in range(nelem):
        # Find reference coordinates of element nodes

        X = coor[elnodes[i, ndnum] - 1, 0]
        Y = coor[elnodes[i, ndnum] - 1, 1]
        # Get global degree of freedom numbers per element
        pg = np.asarray(np.zeros([DofPerEl], int))
        pg[::2] = np.asarray(2 * elnodes[i, 1:(1 + NodesPerEl)] - 2)
        pg[1::2] = np.asarray(2 * elnodes[i, 1:(1 + NodesPerEl)] - 1)
        # Get current guess for nodal displacements
        XY = np.matrix([X, Y]).T

        if ElType=='q8' and IntegrationType=='FI':
            k_elem = Quad8_Stiffness(XY,matD,DD['THICKNESS'][0][0])
        elif ElType=='q8' and IntegrationType=='RI':
            k_elem = Quad8R_Stiffness(XY,matD,DD['THICKNESS'][0][0])
        elif ElType=='q4' and IntegrationType=='FI':
            k_elem = Quad4_Stiffness(XY,matD,DD['THICKNESS'][0][0])
        elif ElType=='q4' and IntegrationType=='RI':
            k_elem = Quad4R_Stiffness(XY,matD,DD['THICKNESS'][0][0])
        # Assemble k_elem into sparse k_global using vectors
        k = TotKvec * i + np.matrix(range(TotKvec)).T
        k_elem = np.matrix(k_elem.flatten()).T
        El_pg = Allpg[pg].flatten()
        row_vec.flat[k] = El_pg[rowpos]
        col_vec.flat[k] = El_pg[colpos]
        stiff_vec[k] = k_elem[0:TotKvec, 0]

    Kff_pos = np.flatnonzero((row_vec>0)&(col_vec>0))
    Kff_row = row_vec[Kff_pos]-1
    Kff_col = col_vec[Kff_pos]-1
    Kff_stf = stiff_vec[Kff_pos]
    Kff = coo_matrix((Kff_stf, (Kff_row,Kff_col)),
                          shape=(len(fdof),len(fdof)))
    Kff = Kff.tocsr()

    Kfp_pos = np.flatnonzero((row_vec>0)&(col_vec<0))
    Kfp_row = row_vec[Kfp_pos]-1
    Kfp_col = -col_vec[Kfp_pos]-1
    Kfp_stf = stiff_vec[Kfp_pos]
    Kfp = coo_matrix((Kfp_stf, (Kfp_row,Kfp_col)),
                          shape=(len(fdof),len(pdof)))
    Kfp = Kfp.tocsr()

    Kpp_pos = np.flatnonzero((row_vec<0)&(col_vec<0))
    Kpp_row = -row_vec[Kpp_pos]-1
    Kpp_col = -col_vec[Kpp_pos]-1
    Kpp_stf = stiff_vec[Kpp_pos]
    Kpp = coo_matrix((Kpp_stf, (Kpp_row,Kpp_col)),
                          shape=(len(pdof),len(pdof)))
    Kpp = Kpp.tocsr()

    #Initialize global load vector
    F_ext = np.asarray(np.zeros([2 * nnodes, 1]))
    # Add nodal loads to global load vector
    cload = np.array(DD['CLOAD'])
    for i in range(ncload):
        pos = int(cload[i,0]*2+cload[i,1]-3)
        F_ext[pos, 0] = cload[i,2]

    toc = time.time()
    finish = toc - tic
    print('Done assembling K and F    :','{:8.4f}'.format(finish), 'seconds.')

    Pf  = F_ext[fdof]

    tic = time.time()

    Uf = spsolve(Kff, Pf-Kfp*Up)

    finish = time.time() - tic
    print('Done solving system        :','{:8.4f}'.format(finish), 'seconds.')

    tic = time.time()
    # Place Uf into U
    U[fdof,0] = Uf.T

    # Get support reactions
    Fp = np.matrix(Kfp.T.dot(Uf)).T + Kpp.dot(Up)
    toc = time.time()
    finish = time.time() - tic
    print('Done computing reactions   :','{:8.4f}'.format(finish), 'seconds.')

    tic = time.time()
    Gauss_pos = 1./np.sqrt(3.0)  # Gauss point location
    #Post process strain and stress
    for i in range(nelem):
        # Find reference coordinates of element nodes
        X = coor[elnodes[i, ndnum] - 1, 0]
        Y = coor[elnodes[i, ndnum] - 1, 1]
        # Get global degree of freedom numbers per element
        pg = np.matrix(np.zeros([DofPerEl, 1], int))
        pg[::2, 0] = np.matrix(2 * elnodes[i, 1:(1 + NodesPerEl)] - 2).T
        pg[1::2, 0] = np.matrix(2 * elnodes[i, 1:(1 + NodesPerEl)] - 1).T
        # Get current guess for nodal displacements
        XY = np.matrix([X, Y]).T
        Gn =0
        for jGauss in range(2):    # 2 by 2 Gauss integration loops
            eta = (-1) ** (jGauss + 1) * Gauss_pos  # Natural coordinate eta
            for iGauss in range(2):
                xi = ((-1) ** (iGauss + 1)) * Gauss_pos   # Natural coordinate xi
                if ElType in ['q8','q8r']:
                    [B, detJ] = B_Quad8(xi, eta, XY)     # B for Eq.(4.83)
                elif ElType=='q4':
                    [B, detJ] = B_Quad4(xi, eta, XY)     # B for Eq.(4.83)
                Gp_strain = B*U[pg,0]
                Gp_stress = matD*Gp_strain
                strain[i,3*Gn:3*Gn+3]=Gp_strain.T
                stress[i,3*Gn:3*Gn+3]=Gp_stress.T
                Gn = Gn + 1

    toc = time.time()
    finish = toc-tic
    print('Done computing stresses    :','{:8.4f}'.format(finish), 'seconds.')
    tic = time.time()
    #Compute nodal loads, Von Mises and Tresca
    [StressOut, StressNode] = nodal_values(elnodes, stress)
    [StrainOut, StrainNode] = nodal_values(elnodes, strain)

    IsoData = np.zeros([len(StrainNode),4])
    for i in range(0,len(StrainNode)):
        for j in range(0,4):
            S_tensor=[[StressNode[i,1+3*j],StressNode[i,3+3*j],0.0],
                      [StressNode[i,3+3*j],StressNode[i,2+3*j],0.0],
                      [0.0,0.0,0.0]]
            [Eig,Dir]=np.linalg.eig(S_tensor)
            IsoData[i,j] = max(Eig)-min(Eig)

    VonMises = calc_von_mises(StressNode, pois, plane)

    Tresca = calc_tresca(StressNode, pois, plane)

    toc = time.time()
    finish = toc - tic
    print('Done post-processing stress:','{:8.4f}'.format(finish), 'seconds.')

    tic = time.time()
    #Write output to text based output file
    write_output_file(file_out, U, displ, Fp, nodes, elnodes, StrainOut,
        StressOut, tic)

    toc = time.time()
    finish = toc - tic
    print('Done writing output        :','{:8.4f}'.format(finish), 'seconds.')

    S11 = StressNode[:,[1,4,7,10]]
    S22 = StressNode[:,[2,5,8,11]]
    S12 = StressNode[:,[3,6,9,12]]

    S11sm = SmoothData(S11,nnodes,elnodes,nelem)
    S22sm = SmoothData(S22,nnodes,elnodes,nelem)
    S12sm = SmoothData(S12,nnodes,elnodes,nelem)
    VMism = SmoothData(VonMises,nnodes,elnodes,nelem)
    IsoSm = SmoothData(IsoData,nnodes,elnodes,nelem)
    AllD  = [VMism,S11sm,S22sm,S12sm]
    AllD  = np.asarray(AllD)
    [trisnew,Xvec,Yvec] = NewTri(nelem,nnodes,elnodes,coor)

    Scale  = 8.0
    Titles=[r'$\sigma_{VonMises}$',r'$\sigma_{xx}$',r'$\sigma_{yy}$',r'$\sigma_{xy}$']
    PlotClearAll()
    PlotGeom(elnodes,coor,nelem,nnodes,MagFac,U)
    for i in range (0,4):
        Data = AllD[i,:]
        LB = min(Data)
        UB = max(Data)
        if (UB-LB<0.001):
            LB=LB-0.0005
            UB=UB+0.0005
        PlotData(Xvec,Yvec,Data,trisnew,Titles[i],LB,UB)

#        PlotIso(Xvec,Yvec,IsoSm,trisnew,'Isochromatic',Scale)

    return U, Fp, VonMises, S11, S22, S12, elnodes, coor, nelem, nnodes, U, fdof, pdof, sp.sparse.csr_matrix.todense(Kff)
