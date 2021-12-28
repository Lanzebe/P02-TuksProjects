"""

Author of code: Stephan Schmidt

This code is intended for and provided to Mechanical and Aeronautical Engineering Students at the University of Pretoria (UP) as part of the course MOW 323 - Simulation-based Design.

Please contact stephan.schmidt@up.ac.za if you have questions regarding the use of the code.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

__version__ = "2021-11-10"

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

def fem_input_checker(node_coordinates,
                  elem_connectivity,
                  prescribed_displacement,
                  prescribed_forces,
                  element_properties_dict):
    success = True

    if prescribed_forces is None:
        print("*** WARNING: There are no forces applied to the system. Only a modal analysis can be performed. ***")
    if prescribed_displacement is None:
        print("*** WARNING: There are no displacements applied to the system. Only a modal analysis can be performed. ***")


    PRESCRIBED_DISPLACEMENT_TEXT = "prescribed displacement matrix"
    PRESCRIBED_FORCES_TEXT = "prescribed forces matrix"
    def array_checker(data,text_input,success):
        is_not_an_array = isinstance(data,np.ndarray) == False

        # This text is added to allow a None for the prescribed forces and displacements:
        if (text_input==str(PRESCRIBED_FORCES_TEXT)) or (text_input == str(PRESCRIBED_DISPLACEMENT_TEXT)):
            if data is None:
                is_not_an_array = False

        if is_not_an_array:
            print("ERROR: The {} should be a numpy array and not a {}".format(text_input,type(data)))
            success = False
        return success

    for text_input,data in [["nodal coodinate matrix",node_coordinates],
                 ["element connectivity matrix",elem_connectivity],
                 [PRESCRIBED_DISPLACEMENT_TEXT,prescribed_displacement],
                 [PRESCRIBED_FORCES_TEXT,prescribed_forces]]:
        success = array_checker(data,text_input,success)

    try:

        if prescribed_displacement is not None:
            try:
                rows,columns = prescribed_displacement.shape # This is a test.
            except Exception as e:
                print(e)
                raise ValueError("The matrix should be two-dimensional")
            assert len(prescribed_displacement.shape) == 2 # This is a test.

            prescribed_displacement[:,1] # This is a test.
    except Exception as e:
        print("The prescribed displacement should be a N x 3 numpy array.")
        raise ValueError(e)

    try:
        if prescribed_forces is not None:
            try:
                rows,columns = prescribed_forces.shape # This is a test.
            except Exception as e:
                print(e)
                raise ValueError("The matrix should be two-dimensional")
            assert len(prescribed_forces.shape) == 2 # This is a test.
            prescribed_forces[:,1] # This is a test.
    except Exception as e:
        print("The prescribed displacement should be a N x 3 numpy array.")
        raise ValueError(e)

    shape_node = node_coordinates.shape

    nnodes = node_coordinates.shape[0]
    ec = np.sort(np.unique(elem_connectivity.reshape(elem_connectivity.size)))

    if np.sum(np.diff(np.unique(ec))==2) != 0:
        success = False
        print("ERROR: It seems that some nodes are not connected")
        print("Please check the element connectivity matrix.")

    if prescribed_displacement is not None:
        if np.max(prescribed_displacement[:,0]) >= (node_coordinates.shape[0]):
            success = False
            print("ERROR: There are nodes specified in the prescribed displacement that does not exist.")

    if prescribed_forces is not None:
        if np.max(prescribed_forces[:,0]) >= (node_coordinates.shape[0]):
            success = False
            print("ERROR: There are nodes specified in the prescribed force that does not exist.")

    if len(shape_node) != 2:
        success = False
        print("ERROR: Shape of node coordinates:",shape_node)
        print("The nodal coordinates should be a matrix with a shape N x 2.")

    if shape_node[1] != 2:
        success = False
        print("ERROR: Shape of node coordinates:",shape_node)
        print("The nodal coordinates should be a matrix with a shape N x 2.")

    shape_elem_conn = elem_connectivity.shape

    if len(shape_elem_conn) != 2:
        success = False
        print("ERROR: Shape of element connectivity matrix:",shape_elem_conn)
        print("The element connectivity matrix should be a matrix with a shape N x M.")

    if prescribed_displacement is not None:
        shape_disp = prescribed_displacement.shape
        if len(shape_disp) != 2:
            success = False
            print("ERROR: Shape of the prescribed displacement matrix:",shape_disp)
            print("The prescribed displacement matrix should be a matrix with a shape N x 3.")

        if shape_disp[1] != 3:
            success = False
            print("ERROR: Shape of the prescribed displacement matrix:",shape_disp)
            print("The prescribed displacement matrix should be a matrix with a shape N x 3.")

    if prescribed_forces is not None:
        shape_force = prescribed_forces.shape
        if len(shape_force) != 2:
            success = False
            print("ERROR: Shape of the prescribed force matrix:",shape_force)
            print("The prescribed force matrix should be a matrix with a shape N x 3.")

        if shape_force[1] != 3:
            success = False
            print("ERROR: Shape of the prescribed force matrix:",shape_force)
            print("The prescribed force matrix should be a matrix with a shape N x 3.")

    if prescribed_forces is not None:
        unique_force = np.unique(prescribed_forces[:,:2], axis=0)

        if unique_force.shape[0] != prescribed_forces.shape[0]:
            success = False
            print("ERROR: There are duplicates in the prescribed forces matrix")

    else:
        unique_forces = np.array([[]])
        #prescribed_displacement = None

    if prescribed_displacement is not None:
        unique_displ = np.unique(prescribed_displacement[:,:2], axis=0)

        if unique_displ.shape[0] != prescribed_displacement.shape[0]:
            success = False
            print("ERROR: There are duplicates in the prescribed displacement matrix")
    else:
        unique_displ = np.array([[]])
        #prescribed_displacement = None

    #unique_displ = np.unique(prescribed_displacement[:,:2], axis=0)

    if (prescribed_displacement is not None) and (prescribed_forces is not None):
        force_displacement = np.vstack([prescribed_displacement,prescribed_forces])
        unique_forcedisp = np.unique(force_displacement[:,:2], axis=0)
        if unique_forcedisp.shape[0] != force_displacement.shape[0]:
            success = False
            print("ERROR: There are duplicates in the prescribed forces and prescribed displacement. You cannot assign a force and displacement to the same degree-of-freedom")

    element_properties = element_properties_dict.keys()
    size_list = []
    for ikey in element_properties:
        size_of_element = element_properties_dict[ikey].shape[0]
        if size_of_element != shape_elem_conn[0]:
            print("ERROR: The property {} of the element property dictionary is inconsistent with the size:".format(ikey),size_of_element,"should be",shape_elem_conn[0])
            success = False

            success = array_checker(element_properties_dict[ikey],"element matrix property, {},".format(ikey),success)

    if success == False:
        raise ValueError("The input data need to be fixed before the code can continue.")

    return True

class FEM_Truss_2D:

    """
        This is a class that allows you to model a system using truss elements.

        Parameters
        ------------

        node_coordinates: array_like
            This is an N x 2 array with the N nodal coordinates. The first column is the x-coordinate, the second column is the y-coordinate.
            The row number indicates the node number.

            For example,

            [[0,0],[1,2]]

            [0,0] is the coordinate of node 0 (x = 0, y = 0)
            [1,2] is the coordinate of node 1 (x = 1, y = 2)

        elem_connectivity: array_like
            This is an M x 2 array which indicates which nodes are connected to form an element.
            The row number indicates the element number.

            For example,

            [[0,1],[1,2]]

            indicates that element 0 consists of node 0 and node 1; and
            element 1 consists of node 1 and node 2.

        prescribed_displacement: array_like
            This is an S x 2 array that indicates the prescribed displacements of the nodes.

            [[node_number,direction,value],[node_number,direction,value]]

            node_number indicates the node number.

            direction = 0 means it is in the x-direction.
            direction = 1 means it is in the y-direction.

            value indicates the value of the prescribed displacement on the node.

            For example,

            [[0,1,0],[1,0,0]]

            indicates

            - [0,1,0] means node 0 has a prescribed displacement of 0 in the 1 direction (y-direction).
            - [1,0,2] means node 1 has a prescribed displacement of 2 in the 0 direction (x-direction).

        prescribed_forces: array_like
            This is an S x 2 array that indicates the prescribed forces of the nodes.

            [[node_number,direction,value],[node_number,direction,value]]

            direction = 0 means it is in the x-direction.

            direction = 1 means it is in the y-direction.

            value indicates the value of the prescribed forces on the node.

            For example,

            [[0,1,0],[1,0,0]]

            indicates

            - [0,1,0] means node 0 has a prescribed force of 0 in the 1 direction (y-direction).
            - [1,0,2] means node 1 has a prescribed force of 2 in the 0 direction (x-direction).

            If no forces are applied to the structure supply None.

        element_properties_dict: dictionary

            This is a dictionary that indicates the properties of the elements.

            The following keys are used:

            - "E": Young's modulus or Modulus of Elasticity
            - "A": The cross-sectional area of the element.
            - "rho": The density of the material (optional depending on the analysis)

            For example,

            {"E": [1,2], "A": [3,4]}

            specifies the properties of the two elements:

            Element 0 has an E = 1 (Young's modulus) and an A = 3 (Area)
            Element 1 has an E = 2 (Young's modulus) and an A = 4 (Area)

    """

    def __init__(self,node_coordinates,
                      elem_connectivity,
                      prescribed_displacement,
                      prescribed_forces,
                      element_properties_dict):

        sucess = fem_input_checker(node_coordinates,
                          elem_connectivity,
                          prescribed_displacement,
                          prescribed_forces,
                          element_properties_dict)

        assert "E" in element_properties_dict, "The property E is not in the element dictionary"
        #assert "I" in element_properties_dict, "The property I is not in the element dictionary"
        assert "A" in element_properties_dict, "The property A is not in the element dictionary"

        self.prescribed_displacement = prescribed_displacement
        self.prescribed_forces = prescribed_forces
        self.nodecoor       = node_coordinates
        self.elemconn       = elem_connectivity
        self.elemmatprop    = element_properties_dict
        self.nelem          = elem_connectivity.shape[0]
        self.nnodes         = node_coordinates.shape[0]
        self.ndof           = 2


        if prescribed_displacement is not None:
            if len([1 for i in prescribed_displacement[:,1] if i in [0,1]]) != prescribed_displacement.shape[0]:
                raise ValueError("All prescribed displacement directions should be 0 (x) or 1 (y). For example, you cannot supply negative directions.")
            self._supplied_displacements = True
        else:
            self._supplied_displacements = False

        if prescribed_forces is not None:
            self._supplied_loads = True
            if len([1 for i in prescribed_forces[:,1] if i in [0,1]]) != prescribed_forces.shape[0]:
                raise ValueError("All prescribed force directions should be 0 (force in x), 1 (force in y). For example, you cannot supply negative directions.")
        else:
            self._supplied_loads = False

        self._calculate_L_theta()

        self._solved_system_bool = False

    def _transformation_rad(self,theta_rad):
        """
        Rotate an element theta_rad radians counterclockwise.

        returns a transformation matrix.
        """
        C = np.cos(theta_rad)
        S = np.sin(theta_rad)
        T = np.array([[C,S,0,0],[-S,C,0,0],[0,0,C,S],[0,0,-S,C]])

        return T

    def _transformation_deg(self,theta_deg):
        """
        Rotate an element at theta_deg degrees counterclockwise.

        returns a transformation matrix.
        """
        theta_rad = np.deg2rad(theta_deg)

        T = self._transformation_rad(theta_rad)

        return T

    def _stiffness_local(self,A,E,L):
        """
        # The stiffness of an unrotated element.

        returns the local stiffness matrix.
        """
        k = A * E / L * np.array([[1,0,-1,0],[0,0,0,0],[-1,0,1,0],[0,0,0,0]])

        return k

    def _calculate_L_theta(self):
        """
        Calculates the length and angles of different elements in the mesh.

        inputs
        --------

        None

        returns
        --------

        matrix: array_like

            The first column indicates the length of each element. The second column indicates the angle of rotation of each element.
        """

        nodes  = self.nodecoor#  = nodes
        elem   = self.elemconn# = elemconn

        nelem = elem.shape[0]

        nnodes = (nodes.shape[0])

        ndof = int(nnodes*2)
        K = np.zeros([ndof,ndof])
        Ltheta = []
        for i_elem in range(nelem):

            node_num_A = elem[i_elem][0]
            node_num_B = elem[i_elem][1]

            nodeA = nodes[node_num_A]
            nodeB = nodes[node_num_B]

            L = np.linalg.norm(nodeA - nodeB)

            dydx = nodeB - nodeA

            dx = dydx[0]
            dy = dydx[1]

            theta_rad = np.arctan2(dy,dx)

            Ltheta.append([L,theta_rad])

        self.element_L_theta = np.array(Ltheta)

    def _mass_local(self,rho,A,L):
        """
        The mass matrix of an unrotated element.

        returns the mass matrix.
        """
        m = rho * A * L/6 * np.array([[2,0,1,0],[0,0,0,0],[1,0,2,0],[0,0,0,0]])

        return m

    def mass_global(self):
        """
        The global mass matrix of the system under consideration is returned by this function.

        inputs
        ----------
        None

        returns
        ----------
        Mass matrix.
        """

        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        #nelem = self.nelem

        elem_E = self.elemmatprop["E"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)

        elem_L_theta = self.element_L_theta

        nelem = elemconn.shape[0]

        nnodes = (nodecoor.shape[0])

        ndof = int(nnodes*self.ndof)
        M = np.zeros([ndof,ndof])
        for i_elem in range(nelem):

            node_num_A = elemconn[i_elem][0]
            node_num_B = elemconn[i_elem][1]

            nodeA = nodecoor[node_num_A,:]
            nodeB = nodecoor[node_num_B,:]

            L = np.linalg.norm(nodeA - nodeB)

            dydx = nodeB - nodeA

            dx = dydx[0]
            dy = dydx[1]

            theta_rad = np.arctan2(dy,dx)

            idx_conn = np.array([np.hstack([node_num_A*2,node_num_A*2+1,node_num_B*2,node_num_B*2+1])])

            A = elem_A[i_elem]
            E = elem_E[i_elem]
            rho = elem_rho[i_elem]
            L,theta_rad = self.element_L_theta[i_elem,:]

            k = self._mass_local(rho,A,L)

            T = self._transformation_rad(theta_rad)

            k_new = np.dot(T.T,np.dot(k,T))

            M[idx_conn.T,idx_conn] += k_new

        return M

    def stiffness_global(self):
        """
        The global stiffness matrix of the system under consideration.

        inputs
        --------

        None

        returns
        --------

        matrix: array_like
            stiffness matrix associated with the whole system. The global degrees-of-freedom are returned.
        """

        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        #nelem = self.nelem

        elem_E = self.elemmatprop["E"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)

        # elem_A = self.elem_A# = elem_A
        # elem_E = self.elem_E# = elem_E
        # nodes  = self.nodes#  = nodes
        # elem   = self.elemconn# = elemconn

        elem_L_theta = self.element_L_theta

        nelem = elemconn.shape[0]

        nnodes = (nodecoor.shape[0])

        ndof = int(nnodes*self.ndof)
        K = np.zeros([ndof,ndof])
        for i_elem in range(nelem):

            node_num_A = elemconn[i_elem][0]
            node_num_B = elemconn[i_elem][1]

            nodeA = nodecoor[node_num_A,:]
            nodeB = nodecoor[node_num_B,:]

            L = np.linalg.norm(nodeA - nodeB)

            dydx = nodeB - nodeA

            dx = dydx[0]
            dy = dydx[1]

            theta_rad = np.arctan2(dy,dx)

            idx_conn = np.array([np.hstack([node_num_A*2,node_num_A*2+1,node_num_B*2,node_num_B*2+1])])

            A = elem_A[i_elem]
            E = elem_E[i_elem]
            L,theta_rad = self.element_L_theta[i_elem,:]

            k = self._stiffness_local(A,E,L)

            T = self._transformation_rad(theta_rad)

            k_new = np.dot(T.T,np.dot(k,T))

            K[idx_conn.T,idx_conn] += k_new

        return K

    def solve(self):
        """
        Solve the system of equations for the prescribed displacements and loads.

        inputs
        ---------
        None

        returns
        ---------

        solution_dict: dictionary

            This dictionary contains the following keys:
                "displacements": array_like
                    The displacement solution for the system is given in this field.
                    The displacements of the nodes are presented. If we denote the data of the displacements as u, the notation is as follows:
                    u[0] = displacement of node0 in direction 0
                    u[1] = displacement of node0 in direction 1
                    u[2] = displacement of node1 in direction 0
                    u[3] = displacement of node1 in direction 1

                "loads": array_like
                    The loads for the system is given in this field.
                    The loads applied to the nodes are presented here. This will contain applied loads and reaction forces. If we denote the data of the loads as f, the notation is as follows:
                    f[0] = load of node0 in direction 0
                    f[1] = load of node0 in direction 1
                    f[2] = load of node1 in direction 0
                    f[3] = load of node1 in direction 1
        """
        if self._supplied_displacements == False:
            raise ValueError("Cannot solve the system. No displacements are specified.")
        if self._supplied_loads == False:
            raise ValueError("Cannot solve the system. No loads are specified.")

        K = self.stiffness_global()

        ddof = np.arange(K.shape[0])

        pdisp = self.prescribed_displacement

        pforce = self.prescribed_forces
        pdof = pdisp[:,0]*2 + pdisp[:,1]
        fdof = np.array(list(set(ddof) - set(pdof)))

        F = np.zeros(K.shape[0])
        u = np.zeros(K.shape[0])

        F[(pforce[:,0]*2 + pforce[:,1]).astype(int)] = pforce[:,2]

        Kff = K[np.array([fdof]).T,fdof]
        Kfp = K[np.array([fdof]).T,pdof]
        Kpp = K[np.array([pdof]).T,pdof]

        Up = u[pdof]
        x = F[fdof] - np.dot(Kfp,Up)

        Uf = np.linalg.solve(Kff,x)
        Fp = np.dot(Kfp.T,Uf) + np.dot(Kpp,Up)

        u[fdof] = Uf

        F[pdof] = Fp
        #self.calculate_stress()

        self._solved_system_bool = True

        self.solution_dict = {"displacements": u,
                              "loads": F}

        return self.solution_dict

    def post_print_nodal_data(self,nodal_data):
        """
            Print the nodal data using data from the global degrees-of-freedom.

            inputs:
            ---------

            nodal_data: array_like
                The global degrees-of-freedom_data.

            returns:
            ----------
            None
        """
        K = self.stiffness_global()
        nodal_data = nodal_data.squeeze()
        assert (K.shape[0]) == (nodal_data.shape[0]), "The nodal data does not have the expected shape of {}".format(K.shape[0])

        for idx in range(int(len(nodal_data)/2)):
            ux = nodal_data[2*idx]
            uy = nodal_data[2*idx+1]

            print(" "*20,"Node",idx," "*20)
            print("x: {:f}, y: {:f}".format(ux,uy))

    def make_nodal_data_global(self,nodal_data_free,nodal_data_prescribed):
        """
            Convert data associated with the free degrees-of-freedom and prescribed degrees-of-freedom to global degrees-of-freedom.

            inputs:
            ---------
                nodal_data_free: array_like
                    The data associated with the free degrees-of-freedom.

                nodal_data_prescribed: array_like
                    The data associated with the prescribed degrees-of-freedom.

            returns:
            ---------
                nodal_data_global: array_like
                    The data associated with the global degrees-of-freedom.
        """
        K = self.stiffness_global()
        idx_f,idx_p = self._partition_vector(np.arange(K.shape[0]))

        data = np.zeros(K.shape[0])

        if (idx_f.size) != (nodal_data_free.size):
            raise ValueError("The free degree of freedom data should have the shape ({},)".format(idx_f.size))
        if (idx_p.size) != (nodal_data_prescribed.size):
            raise ValueError("The prescribed degree of freedom data should have the shape ({},)".format(idx_p.size))

        data[idx_f] = nodal_data_free
        data[idx_p] = nodal_data_prescribed

        return data

    def post_print_element_data(self,array):

        if isinstance(array,(np.ndarray,list)) == False:
            raise ValueError("The input data should be an array.")
        array = list(array)
        for idx in range(len(array)):
            print(" "*20,"Element",idx," "*20)
            print("Data: {:f}".format(array[idx]/1E6))


    def post_print_solution(self,displacement=None):
        """
        This function prints the displacements, forces and stresses. The units are consistent with the units supplied in the input of the class.

        inputs
        ---------

        None required.

        Optional: displacement of the global degrees-of-freedom.

        returns
        ---------
        None

        """

        if displacement is None:
            if self._solved_system_bool == False:
                raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
            u = self.solution_dict["displacements"]
            F = self.solution_dict["loads"]
        else:
            u = np.array(displacement)
            K = self.stiffness_global()
            if K.shape[0] != u.shape[0]:
                raise ValueError("The displacement input is the wrong shape.")
            F = np.dot(K,u)


        for idx in range(int(len(u)/2)):
            ux = u[2*idx]
            uy = u[2*idx+1]
            Fx = F[2*idx]
            Fy = F[2*idx+1]

            print(" "*20,"Node",idx," "*20)
            print("u_x: {:f} x 10^-3, u_y: {:f} x 10^-3, F_x: {:f}, F_y: {:f}".format(ux*1000,uy*1000,Fx,Fy))

        dict_post = self.post_stress_strain(displacement=u)
        stress = dict_post["stress"]
        strain = dict_post["strain"]

        print("")
        for idx in range(len(stress)):
            print(" "*20,"Element",idx," "*20)
            print("stress: {:f} x 10^6, strain: {:f}".format(stress[idx]/1E6,strain[idx]))

    def make_force_global(self,node,direction,force_magnitude):
        """
        inputs
        -------

        node: integer
            The node where the foce will be applied.

        direction: integer
            The direction in which the force will be applied.

        force_magnitude: integer
            The magnitude of the force.

        Returns
        ---------

        F: array_like
            The global force array that has the same dimensions as the global degrees-of-freedom of the structure

        Ff: array_like
            The force associated with the free degrees-of-freedom. This is used in dynamic analyses.

        Fp: array_like
            The force associated with prescribed degrees-of-freedom. This contains the reaction forces.
        """
        direction = int(direction)
        self._checks_direction(direction)
        self._checks_node_number(node)

        F = np.zeros(self.nnodes*self.ndof)
        # prescribed_forces = self.prescribed_forces
        # nforces = prescribed_forces.shape[0]
        # for iforce in range(nforces):
        node,dof,value = node,direction,force_magnitude#prescribed_forces[iforce,:]
        nodedof = int(dof + 2 * node)
        F[nodedof] = value

        fdof,pdof = self._get_fdof_pdof()
        return F,F[fdof],F[pdof]

    def _checks_direction(self,direction):
        assert direction in [0,1],"The direction must be 0 or 1."
    def _checks_node_number(self,node):
        assert node in np.arange(self.nnodes),"The node should be between 0 and N-1 (where N is the number of nodes)."
    def _checks_element_number(self,element):
        elemconn   = self.elemconn# = elemconn
        nelem = self.elemconn.shape[0]
        assert element in np.arange(nelem),"The element should be between 0 and N-1 (where N is the number of elements)."

    def _get_fdof_pdof(self):
        """
            Returns the free-degree-of-freedoms and the prescribed-degrees-of-freedom
        """

        prescribed_displacement = self.prescribed_displacement
        # ndisp = prescribed_displacement.shape[0]
        if prescribed_displacement is not None:
            ndisp = prescribed_displacement.shape[0]
        else:
            ndisp = 0

        prescribed_dof = []
        for iforce in range(ndisp):

            node,dof,_ = prescribed_displacement[iforce,:]
            nodedof = int(dof + self.ndof * node)
            prescribed_dof.append(nodedof)

        self.prescribed_dof = prescribed_dof
        pdof = np.array(self.prescribed_dof)
        fdof = np.array(list(set(list(np.arange(self.nnodes*self.ndof))) - set(list(pdof))))
        return fdof,pdof

    #-----
    # prescribed_displacement = self.prescribed_displacement
    #
    # if prescribed_displacement is not None:
    #     ndisp = prescribed_displacement.shape[0]
    # else:
    #     ndisp = 0
    # prescribed_dof = []
    # for iforce in range(ndisp):
    #
    #     node,dof,_ = prescribed_displacement[iforce,:]
    #     nodedof = int(dof + self.ndof * node)
    #     prescribed_dof.append(nodedof)
    #
    # self.prescribed_dof = prescribed_dof
    # pdof = np.array(self.prescribed_dof)
    # fdof = np.array(list(set(list(np.arange(self.nnodes*self.ndof))) - set(list(pdof))))
    #
    # #print(fdof),print(pdof)
    #
    # return fdof,pdof


    def _partition_matrix(self,A):
        """

        returns Aff,Afp,App
        """
        fdof,pdof = self._get_fdof_pdof()
        fdof = np.array([fdof])
        pdof = np.array([pdof])


        if len(pdof.squeeze()) > 0:
            App = A[pdof.T,pdof]
            Afp = A[pdof.T,fdof].T
            Aff = A[fdof.T,fdof]
        else:
            return A,np.array([]),np.array([])

        return Aff,Afp,App

    def _partition_vector(self,B):
        """
        returns vector_f,vector_p
        """
        fdof,pdof = self._get_fdof_pdof()

        pdof = np.array(pdof)
        if len(pdof.squeeze()) > 0:
            return B[fdof],B[pdof]
        else:
            return B,np.array([])


        #return B[fdof],B[pdof]

    def stiffness_global_partition(self):
        """
        The governing equations for a linear elastic solid mechanics problem can be written in the following form when using the finite element method:

            [K]{u} = {f}

        We can decompose the stiffness matrix in terms of ff, fp, pp, pf degrees of freedom as follows

        [K] = [[Kff], [Kfp]
               [Kpf], [Kpp]]

        Since the matrix is symmetric, Kfp = Kpf.T

        inputs
        -------

        None

        returns
        --------

        Kff: array_like

        Kfp: array_like

        Kpp: array_like

        """
        K = self.stiffness_global()

        Kff,Kfp,Kpp = self._partition_matrix(K)

        return Kff,Kfp,Kpp

    def mass_global_partition(self):
        """
        We can decompose the mass matrix in terms of ff, fp, pp, pf degrees of freedom as follows

        [M] = [[Mff], [Mfp]
               [Mpf], [Mpp]]

        Since the matrix is symmetric, Mfp = Mpf.T

        inputs
        -------

        None

        returns
        --------

        Mff: array_like

        Mfp: array_like

        Mpp: array_like
        """

        M = self.mass_global()

        Mff,Mfp,Mpp = self._partition_matrix(M)

        return Mff,Mfp,Mpp

    def dynamic_make_matrices(self):
        """
        In a dynamic analysis, the following governing equation is solved for an undamped system:

        Mff x'' + Kff x = Ff(t)

        This function returns Mff and Kff.

        inputs
        ----------

        None

        returns
        -----------

        Mff: array_like
             The mass matrix corresponding to the free degrees-of-freedom.

        Kff: array_like
             The stiffness matrix corresponding to the free degrees-of-freedom.

        """
        Mff,_,_ = self.mass_global_partition()
        Kff,_,_ = self.stiffness_global_partition()

        return Mff,Kff

    def dynamic_make_global_displacement(self,displacement_free):
        """

        This function transforms the displacements associated with the free degrees-of-freedom into the global displacement vector.

        inputs
        ---------

        displacement_free: array_like
            This is the displacement associated with the free-degrees-of-freedom. The order should be consistent with the notation used in this code.

        WARNING: It is assumed that the prescribed displacements are zero.

        returns
        ----------
        global displacement: array_like
            This is the global displacement vector. If we denote the displacements as u, the notation is as follows:
                    u[0] = displacement of node0 in direction 0
                    u[1] = displacement of node0 in direction 1
                    u[2] = displacement of node1 in direction 0
                    u[3] = displacement of node1 in direction 1

        """

        fdof,pdof = self._get_fdof_pdof()

        u = np.zeros(len(fdof) + len(pdof))

        u[fdof] = displacement_free

        return u

    def eigen(self,number_of_eigenvalues=10):
        """
        This performs an eigenvalue analysis.


        inputs:
        -----------

        number_of_eigenvalues: integer

            The number of modes that need to be calculated.

        returns:
        -----------
        results: dictionary with the following keys

            "natural_frequencies_hz": Natural frequencies in Hz.
            "modes": The eigenvectors associated with the free-degrees of freedom.

            Not important:
            "coor_global": The global coordinates.
            "index_f": Indices associated with the f-dofs.
            "index_p": Indices associated with the p-dofs.

        data["natural_frequencies_hz"] = fn #The natural frequencies of the system are in Hz
        data["modes"] = V                #Raw modes
        data["coor_global"] = XY
        data["index_p"] = idx_p
        data["index_f"] = idx_f

        """
        K = self.stiffness_global()
        M = self.mass_global()

        if (self._supplied_displacements == True):
            Kff,Kfp,Kpp = self._partition_matrix(K)
            Mff,Mfp,Mpp = self._partition_matrix(M)
        else:
            print("Using global matrices.")
            Kff = np.array(K)
            Mff = np.array(M)

        idx_f,idx_p = self._partition_vector(np.arange(K.shape[0]))

        import scipy.linalg as scilin
        number_of_eigenvalues = np.min([number_of_eigenvalues,Kff.shape[0]-1])
        eigvals = (0,number_of_eigenvalues)
        try:
            D,V = scilin.eigh(
                Kff,
                b=Mff,
                lower=True,
                eigvals_only=False,
                overwrite_a=False,
                overwrite_b=False,
                turbo=True,
                eigvals=None,
                type=1,
                check_finite=True,
                subset_by_index=eigvals,
                subset_by_value=None,
                driver=None,
            )
        except Exception as e:
            print("*"*100)
            print("The eigenvalue decomposition does not work. This might be due to an older version of scipy.linalg being used.")
            print("The error message that was received is:",e)
            print("Trying a simplified eigenvalue decomposition....")
            print("*"*100)
            D,V = scilin.eigh(
                Kff,
                b=Mff,
                lower=True
            )
            print("It seems to be working.")


        fn = np.sqrt(np.abs(D)) / (2*np.pi)

        data = {}

        y = self.nodecoor[:,1]
        x = self.nodecoor[:,0]

        XY = np.hstack([x,y])

        data["natural_frequencies_hz"] = fn #The natural frequencies of the system are in Hz
        data["modes"] = V                #Raw modes
        data["coor_global"] = XY
        data["index_p"] = idx_p
        data["index_f"] = idx_f
        return data

    def plot_modes(self,magnification=1,number_of_modes=3):
        """

        This function plots the modes.

        inputs:
        ---------
        magnification: integer,float
            Magnifies the displacement of the modes with this factor.

        number_of_modes: integer
            NUmber of modes that need to be calculated.

        returns:
        -----------
        None


        """

        eigen_output = self.eigen(number_of_eigenvalues=number_of_modes+1)

        modes = eigen_output["modes"]
        Xg    = eigen_output["coor_global"]
        idx_p = eigen_output["index_p"]
        idx_f = eigen_output["index_f"]

        #print(modes.shape)
        for eigennumber in range(np.min([modes.shape[1],number_of_modes])):
            V = modes[:,eigennumber]
            #print(V.shape)
            Vg = np.zeros(Xg.shape[0])
            #print(Vg.shape)
            Vg[idx_f] = V
            #print(idx_f,Vg)
            self.plot_system(fignum=None,displacement_vector=Vg,magnification=magnification)
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.title(str(eigen_output["natural_frequencies_hz"][eigennumber]) + " Hz")

    def post_stress_strain(self,displacement=None):
        """
        Post-processing calculation for the stress and strain.

        inputs:
        ---------

        displacements: array_like (optional input)

            This displacement should be the global degrees-of-freedom and should be in the order as expected by the code.

        returns:
        ---------

            results: dictionary

                strain: array_like
                    The strain in the order of the elements.

                    strain[0] = strain of element 0
                    strain[1] = strain of element 1
                    strain[2] = strain of element 2
                    ...

                stress: array_like
                    The stress in the order of the elements.

                    stress[0] = stress of element 0
                    stress[1] = stress of element 1
                    stress[2] = stress of element 2
                    ...

        """
        if displacement is None:
            if self._solved_system_bool == False:
                raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
            u = self.solution_dict["displacements"]
        else:
            u = np.array(displacement)
        elem_E = self.elemmatprop["E"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)

        nodecoor  = self.nodecoor#  = nodes
        elemconn   = self.elemconn# = elemconn
        nelem = elem_E.shape[0]

        strain = np.zeros(nelem)
        stress = np.zeros(nelem)
        for i_elem in range(nelem):

            node_num_A = elemconn[i_elem][0]
            node_num_B = elemconn[i_elem][1]

            nodeA = nodecoor[node_num_A,:]
            nodeB = nodecoor[node_num_B,:]

            uA = u[[node_num_A*2,node_num_A*2+1]]
            uB = u[[node_num_B*2,node_num_B*2+1]]

            L = np.linalg.norm(nodeB - nodeA)
            #-------------------------------------------------------------------
            # Small strain assumption:
            disp = uB - uA #Displacement of the vector in the global space
            unitv = (nodeB - nodeA)/np.linalg.norm(nodeB - nodeA) #Directional vector of element
            disp = np.dot(unitv,disp) #Displacement of element that is related to strain
            #-------------------------------------------------------------------

            strain[i_elem]  = disp/L
            stress[i_elem]  = strain[i_elem] * elem_E[i_elem]

        return {"strain": strain,"stress": stress}

    def post_get_displacement(self,node,direction):
        """
        Returns the displacement at a node.

        inputs:
        ---------

        node: integer
            The node number starting at 0

        direction: integer
            The direction (0 = x, 1 = y).

        returns:
        ----------

        displacement: float
            The displacement of the node in the specified direction.
        """
        self._checks_direction(direction)
        self._checks_node_number(node)

        nnodes = self.nodecoor.shape[0]

        node = int(node)
        if (node >= nnodes) | (node <0 ):
            raise ValueError("The node does not exist")
        direction = int(direction)
        if direction not in [0,1]:
            raise ValueError("The direction does not exist.")

        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        return self.solution_dict["displacements"][int(node*2+direction)]

    def post_get_force(self,node,direction):
        self._checks_direction(direction)
        self._checks_node_number(node)

        """
            Returns the force at a node.

            inputs:
            ---------

            node: integer
                The node number starting at 0

            direction: integer
                The direction (0 = x, 1 = y).

            returns:
            ----------

            force: float
                The force applied to the node in the specified direction.

        """
        nnodes = self.nodecoor.shape[0]
        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")

        node = int(node)
        if (node >= nnodes) | (node <0 ):
            raise ValueError("The node does not exist")
        direction = int(direction)
        if direction not in [0,1]:
            raise ValueError("The direction does not exist.")

        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        return self.solution_dict["loads"][int(node*2+direction)]

    def post_get_stress(self,element_number):
        self._checks_element_number(element_number)
        """
        This function gets the stress of an element.

        input:
        --------

        element_number: integer
            The element number of interest. Numbering starts at 0.

        returns:
        --------

        stress: float
            The stress in the specified element.
        """
        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        dict_ss = self.post_stress_strain()

        S = dict_ss["stress"]

        element_number = int(element_number)

        if element_number >= self.elemconn.shape[0]:
            raise ValueError("The specified element number does not exist.")

        return S[element_number]

    def get_mass(self):
        """
        This function calculates the mass of the structure.


        inputs
        -------
        None

        returns
        --------

        mass: float
            The mass of the structure.

        """
        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        #nelem = self.nelem

        elem_E = self.elemmatprop["E"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)

        if elem_rho is None:
            raise ValueError("The density of the elements should be defined.")

        # elem_A = self.elem_A# = elem_A
        # elem_E = self.elem_E# = elem_E
        # nodes  = self.nodes#  = nodes
        # elem   = self.elemconn# = elemconn

        elem_L_theta = self.element_L_theta

        nelem = elemconn.shape[0]
        #
        # nnodes = (nodecoor.shape[0])
        #
        # ndof = int(nnodes*self.ndof)
        # K = np.zeros([ndof,ndof])

        mass = 0
        for i_elem in range(nelem):
            A = elem_A[i_elem]
            E = elem_E[i_elem]
            rho = elem_rho[i_elem]
            L,theta_rad = self.element_L_theta[i_elem,:]
            mass += rho * A * L

        return mass


    def plot_system(self,fignum=None,
                         magnification=1,
                         displacement_vector=None,
                         show_node_numbers=True,
                         show_displacement=False,
                         print_text=None,
                         scaling_text=1,
                         show_element_numbers=True,
                         scale_cross_area=False,
                         figsize=None):

        """

        The function plots the system.

        inputs: All inputs are optional
        ----------

            fignum: integer
                The figure number. Default (None)

            magnification: float
                The magnification of the displacement (if applicable). Default (1).

            displacement_vector: array_like
                If the displacement vector is specified, it will superimpose it on the plot. It should be the global dof displacements.

            show_node_numbers: bool
                Plot the node numbers (True) or not (False).

            show_element_numbers: bool
                Plot the element numbers (True) or not (False)

            show_displacement:
                Show the displacement (True) or not (False). Note that if the displacement_vector is not specified, it will try to find the solved displacement field when set to true.
                Therefore, the user first needs to solve the FEM to use this functionality.

        returns:
        ------------
        None

        """


        elem_text_data = None

        if (show_displacement == True):
            if (self._solved_system_bool == False) and (displacement_vector is None):
                raise ValueError("The system is not solved yet. Run .solve() before applying this function.")

            if displacement_vector is None:
                u = np.array(self.solution_dict["displacements"])
            else:
                u = np.array(displacement_vector)

        if displacement_vector is not None:
            show_displacement = True
            u = np.array(displacement_vector)

        elem_E = self.elemmatprop["E"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)

        nodes  = self.nodecoor#  = nodes
        elem   = self.elemconn# = elemconn

        if isinstance(fignum,(int,float)):
            plt.figure(int(fignum),figsize=figsize)
        else:
            plt.figure(figsize=figsize)

        if show_node_numbers == True:
            dxy_max = np.abs(np.max(nodes,axis = 0) - np.min(nodes,axis = 0))
            for idx_nodes in range(nodes.shape[0]):
                rnd = np.random.rand(2) * dxy_max * 0.01
                plt.text(nodes[idx_nodes,0]+rnd[0],nodes[idx_nodes,1]+rnd[1],"N" + str(idx_nodes))

            # for idx in range(self.prescribed_displacement.shape[0]):
            #     pd = self.prescribed_displacement[idx,:]
            #     print("Node:",pd[0],"Direction",pd[1],"Displacement",pd[2])

        lw_min = 0.5
        lw_max = 3.0
        elem_A = self.elemmatprop["A"]
        min_A,max_A = np.min(elem_A),np.max(elem_A)
        A = np.array([[max_A,1.0],[min_A,1.0]])
        if np.linalg.det(A) > 1E-16:
            scale_line_A_m,scale_line_A_c = np.linalg.solve(A,np.array([lw_max,lw_min]))
        else:
            scale_line_A_m = 0.0
            scale_line_A_c = 1.0


        # min_A,max_A = np.min(elem_A),np.max(elem_A)
        #
        # lw_min = 0.5
        # lw_max = 3.0
        #
        # scale_line_A_m,scale_line_A_c = np.linalg.solve(np.array([[max_A,1],[min_A,1]]),np.array([lw_max,lw_min]))

        for i_elem in range(elem.shape[0]):

            area_i_elem = elem_A[i_elem]
            line_width_element = 1.0
            if scale_cross_area == True:
                line_width_element = area_i_elem * scale_line_A_m + scale_line_A_c

            node_num_A = elem[i_elem][0]
            node_num_B = elem[i_elem][1]

            nodeA = nodes[node_num_A]
            nodeB = nodes[node_num_B]

            if show_element_numbers == True:
                node_vec = (nodeB - nodeA)*0.5 + nodeA
                plt.text(node_vec[0],node_vec[1],"E" + str(i_elem))

            plt.plot([nodeA[0],nodeB[0]],[nodeA[1],nodeB[1]],'b',linewidth=3)
            if show_displacement==True:
                uA = u[[node_num_A*2,node_num_A*2+1]]*magnification
                uB = u[[node_num_B*2,node_num_B*2+1]]*magnification
                dnodeA = np.array(nodeA) + uA
                dnodeB = np.array(nodeB) + uB
                plt.plot([dnodeA[0],dnodeB[0]],[dnodeA[1],dnodeB[1]],'r--',lw=line_width_element)

class FEM_Beam_2D:

    """
        This is a class that allows you to model a system using beam elements (transverse displacement and moment) with axial loading.

        Parameters
        ------------

        node_coordinates: array_like
            This is an N x 2 array with the N nodal coordinates. The first column is the x-coordinate, the second column is the y-coordinate.
            The row number indicates the node number.

            For example,

            [[0,0],[1,2]]

            [0,0] is the coordinate of node 0 (x = 0, y = 0)
            [1,2] is the coordinate of node 1 (x = 1, y = 2)

        elem_connectivity: array_like
            This is an M x 2 array which indicates which nodes are connected to form an element.
            The row number indicates the element number.

            For example,

            [[0,1],[1,2]]

            indicates that element 0 consists of node 0 and node 1; and
            element 1 consists of node 1 and node 2.

        prescribed_displacement: array_like
            This is an S x 2 array that indicates the prescribed displacements of the nodes.

            [[node_number,direction,value],[node_number,direction,value]]

            node_number indicates the node number.

            direction = 0 means it is in the x-direction.
            direction = 1 means it is in the y-direction.
            direction = 2 means it is a rotation in the counterclockwise-direction.

            value indicates the value of the prescribed displacement on the node.

            For example,

            [[0,1,0],[1,0,0]]

            indicates

            - [0,1,0] means node 0 has a prescribed displacement of 0 in the 1 direction (y-direction).
            - [1,0,2] means node 1 has a prescribed displacement of 2 in the 0 direction (x-direction).
            - [1,2,2] means node 1 has a prescribed rotation of 2 in the 2 direction (counterclockwise-direction).

        prescribed_forces: array_like
            This is an S x 2 array that indicates the prescribed forces of the nodes.

            [[node_number,direction,value],[node_number,direction,value]]

            direction = 0 means it is in the x-direction.
            direction = 1 means it is in the y-direction.
            direction = 2 means it is a rotation in the counterclockwise-direction.

            value indicates the value of the prescribed forces on the node.

            For example,

            [[0,1,0],[1,0,0]]

            indicates

            - [0,1,0] means node 0 has a prescribed force of 0 in the 1 direction (y-direction).
            - [1,0,2] means node 1 has a prescribed force of 2 in the 0 direction (x-direction).

        element_properties_dict: dictionary

            This is a dictionary that indicates the properties of the elements.

            The following keys are used:

            - "E": Young's modulus or Modulus of Elasticity
            - "A": The cross-sectional area of the element.
            - "I": The area moment of inertia of the element.
            - "rho": The density of the material (optional depending on the analysis)

            For example,

            {"E": [1,2], "A": [3,4]}

            specifies the properties of the two elements:

            Element 0 has an E = 1 (Young's modulus) and an A = 3 (Area)
            Element 1 has an E = 2 (Young's modulus) and an A = 4 (Area)

    """
    def __init__(self,
                node_coordinates,
                elem_connectivity,
                prescribed_displacement,
                prescribed_forces,
                element_properties_dict):

        success = fem_input_checker(node_coordinates,
                          elem_connectivity,
                          prescribed_displacement,
                          prescribed_forces,
                          element_properties_dict)




        if prescribed_displacement is not None:
            prescribed_displacement = np.array(prescribed_displacement)
            if len([1 for i in prescribed_displacement[:,1] if i in [0,1,2]]) != prescribed_displacement.shape[0]:
                raise ValueError("All prescribed displacement directions should be 0 (x), 1 (y) and 2 (rotation). For example, you cannot supply negative directions.")
            self._supplied_displacements = True
        else:
            self._supplied_displacements = False

        if prescribed_forces is not None:
            prescribed_forces = np.array(prescribed_forces)
            if len([1 for i in prescribed_forces[:,1] if i in [0,1,2]]) != prescribed_forces.shape[0]:
                raise ValueError("All prescribed force directions should be 0 (force in x), 1 (force in y) and 2 (moment around z). For example, you cannot supply negative directions.")
            self._supplied_loads = True
        else:
            self._supplied_loads = False

        # if len([1 for i in prescribed_displacement[:,1] if i in [0,1,2]]) != prescribed_displacement.shape[0]:
        #     raise ValueError("All prescribed displacement directions should be 0 (x), 1 (y) and 2 (rotation). For example, you cannot supply negative directions.")
        # if len([1 for i in prescribed_forces[:,1] if i in [0,1,2]]) != prescribed_forces.shape[0]:
        #     raise ValueError("All prescribed force directions should be 0 (force in x), 1 (force in y) and 2 (moment around z). For example, you cannot supply negative directions.")

        assert "E" in element_properties_dict, "The property E is not in the element dictionary"
        assert "I" in element_properties_dict, "The property I is not in the element dictionary"
        assert "A" in element_properties_dict, "The property A is not in the element dictionary"

        self._solved_system_bool = False
        self.prescribed_displacement = prescribed_displacement
        self.prescribed_forces = prescribed_forces
        self.nodecoor       = node_coordinates
        self.elemconn       = elem_connectivity
        self.elemmatprop    = element_properties_dict
        self.nelem          = elem_connectivity.shape[0]
        self.nnodes         = node_coordinates.shape[0]
        self.ndof           = 3

    def _checks_direction(self,direction):
        assert direction in [0,1,2],"The direction must be 0, 1 or 2."
    def _checks_node_number(self,node):
        assert node in np.arange(self.nnodes),"The node should be between 0 and N-1 (where N is the number of nodes)."
    def _checks_element_number(self,element):
        nelem = self.nelem
        assert element in np.arange(nelem),"The element should be between 0 and N-1 (where N is the number of elements)."

    # self._checks_direction(direction)
    # self._checks_node_number(node)
    # self._checks_element_number(element)

    def plot_system(self,fignum=None,
                         magnification=1,
                         show_node_numbers=True,
                         show_element_numbers=True,
                         show_displacement=False,
                         displacement_vector=None,
                         scale_cross_area=False,
                         scaling_text=1,
                         npoints=100,
                         figsize=None):
        """
            The function plots the system.

            inputs: All inputs are optional
            ----------

                fignum: integer
                    The figure number. Default (None)

                magnification: float
                    The magnification of the displacement (if applicable). Default (1).

                displacement_vector: array_like
                    If the displacement vector is specified, it will superimpose it on the plot. It should be the global dof displacements.

                show_node_numbers: bool
                    Plot the node numbers (True) or not (False).

                show_element_numbers: bool
                    Plot the element numbers (True) or not (False)

                show_displacement:
                    Show the displacement (True) or not (False). Note that if the displacement_vector is not specified, it will try to find the solved displacement field when set to true.
                    Therefore, the user first needs to solve the FEM to use this functionality.

            returns:
            ------------
            None
        """

        lw_min = 1
        lw_max = 5

    # def plot_system(self,fignum=None,
    #                      magnification=1,
    #                      displacement_vector=None,
    #                      show_node_numbers=True,
    #                      show_displacement=False,
    #                      print_text=None,
    #                      scaling_text=1,
    #                      show_element_numbers=True,
    #                      figsize=None):


        nodecoor = self.nodecoor
        prescribed_displacement = self.prescribed_displacement
        if prescribed_displacement is not None:

            unode  = prescribed_displacement[:,0]
            udof   = prescribed_displacement[:,1]
            uvalue = prescribed_displacement[:,2]

        if (show_displacement == True):
            if (self._solved_system_bool == False) and (displacement_vector is None):
                raise ValueError("The system is not solved yet. Run .solve() before applying this function.")

            if displacement_vector is None:
                u = np.array(self.solution_dict["displacements"])
            else:
                u = np.array(displacement_vector)

        if displacement_vector is not None:
            show_displacement = True
            u = np.array(displacement_vector)

        x = nodecoor[:,0]
        y = nodecoor[:,1]

        elemconn = self.elemconn

        elem_A = self.elemmatprop["A"]
        min_A,max_A = np.min(elem_A),np.max(elem_A)
        A = np.array([[max_A,1.0],[min_A,1.0]])
        if np.linalg.det(A) > 1E-16:
            scale_line_A_m,scale_line_A_c = np.linalg.solve(A,np.array([lw_max,lw_min]))
        else:
            scale_line_A_m = 0.0
            scale_line_A_c = 1.0

        plt.figure(figsize=figsize)
        for idx in range(elemconn.shape[0]):

            area_i_elem = elem_A[idx]
            line_width_element = 1.0
            if scale_cross_area == True:
                line_width_element = area_i_elem * scale_line_A_m + scale_line_A_c

            node_i,node_j = elemconn[idx,:]

            x1,y1 = nodecoor[node_i,:]
            x2,y2 = nodecoor[node_j,:]

            plt.plot([x1,x2],[y1,y2],'bo-',lw=line_width_element)
            if show_element_numbers == True:
                xv = np.mean([x1,x2])
                yv = np.mean([y1,y2])
                plt.text(xv,yv,"E" + str(idx))

        plt.xlabel("$x$")
        plt.ylabel("$y$")
        if show_node_numbers == True:
            rnd = np.random.rand(2)*0
            for idx in range(nodecoor.shape[0]):
                plt.text(nodecoor[idx,0]+rnd[0],nodecoor[idx,1]+rnd[1],"N" + str(idx))

        if (show_displacement == True):
            out = self._calc_displacement_and_bending_moment(u,N_points=npoints)

        # plt.plot(,out["y"][0,:] + out["displ_y"][0,:])
        # plt.plot(out["displ_x"][1,:] + out["x"][1,:],out["y"][0,:] + out["displ_y"][1,:])

            #u = self.solution_dict["displacements"]
            for idx in range(elemconn.shape[0]):

                x = magnification * out["displ_x"][idx,:] + out["x"][idx,:]
                y = magnification * out["displ_y"][idx,:] + out["y"][idx,:]

                # node_i,node_j = elemconn[idx,:]
                #
                # x1,y1 = nodecoor[node_i,:]
                # x2,y2 = nodecoor[node_j,:]
                #
                # u1 = u[node_i*3]
                # u2 = u[node_j*3]
                # v1 = u[node_i*3+1]
                # v2 = u[node_j*3+1]
                #
                # x1 += u1 * magnification
                # x2 += u2 * magnification
                # y1 += v1 * magnification
                # y2 += v2 * magnification

                plt.plot(x,y,'r--')

            #plt.title("Warning: Curvature is not yet included in the plots. \nThe nodal displacements are connected with a straight line.")

    def post_get_bending_moment(self,element_number,npoints=100,displacement=None):
        """
            This function calculates the bending moment in a specific element.

            element_number: integer
                This is an integer that indicates the element number, starting from 0 to N-1, where N is the number of elements.

            npoints: integer
                This specifies the number of points that is used to evaluate the moment of the beam.

            displacement: None or array_like
                This is an optional input. If this is supplied, it should be a displacement array with the same dimension as the degrees-of-freedom as the model.

            returns
            ---------

            data: dictionary
            This is a dictionary with the following fields:
            - "moment": array_like
                The bending moment with npoints length.
            - "x_coor": array_like
                The x-coordinate over which the bending moment is defined.
            - "y_coor": array_like
                The y-coordinate over which the bending moment is defined.
        """
        # self._checks_direction(direction)
        # self._checks_node_number(node)
        self._checks_element_number(element_number)

        if displacement is None:
            if self._solved_system_bool == False:
                raise ValueError("The system must first be solved.")
            u = self.solution_dict["displacements"]
        else:
            u = np.array(displacement)
        dict_post_M = self._calc_displacement_and_bending_moment(u,N_points=npoints)

        assert (element_number < dict_post_M["x"].shape[0]), "The element number exceeds the number of elements."

        return {"moment": dict_post_M["moment"][element_number,:],
                "x_coor": dict_post_M["x"][element_number],
                "y_coor": dict_post_M["y"][element_number]}

    def post_plot_bending_moment(self,npoints=200,function=None,figsize=None,displacement=None):
        """
            This function plots the bending moment in the structure.

            npoints: integer
                This specifies the number of points that is used to evaluate the moment of the beam.

            function: None or a function
                This is a function of the bending moment. This can be used to transform the bending moment using an equation. For example if the entire structure has the same cross-sectional properties the flexure formula can be used to calculate
                the stress due to bending, i.e.

                function = M: M * y/I

            figsize: None or array_like
                This can be used to alter the figure size. It should define the x and y-dimensions as follows: (x,y)

            displacement: None or array_like
                This is an optional input. If this is supplied, it should be a displacement array with the same dimension as the degrees-of-freedom as the model.
        """

        from matplotlib.collections import LineCollection
        from matplotlib.colors import ListedColormap, BoundaryNorm

        if displacement is None:
            if self._solved_system_bool == False:
                raise ValueError("The system must first be solved.")
            u = self.solution_dict["displacements"]
        else:
            u = np.array(displacement)
        dict_post_M = self._calc_displacement_and_bending_moment(u,npoints)


        M = dict_post_M["moment"]
        print("Minimum bending moment: {}".format(M.min()))
        print("Maximum bending moment: {}".format(M.max()))

        label_axis = "Bending moment"
        if function is not None:
            M = function(M)
            label_axis = "Transformed bending moment"

            print("Minimum TRANSFORMED bending moment: {}".format(M.min()))
            print("Maximum TRANSFORMED bending moment: {}".format(M.max()))

        x_coor = dict_post_M["x"]
        y_coor = dict_post_M["y"]

        nnodes = M.shape[0] * (M.shape[1] - 1)
        xy = np.zeros([nnodes,2,2])
        M_mat = []
        counter = 0
        for idx in range(M.shape[0]):
            for idy in range(M.shape[1]-1):
                M_mat.append(M[idx,idy:idy+2].mean())
                xy[counter,0,:] = x_coor[idx,idy],y_coor[idx,idy]
                xy[counter,1,:] = x_coor[idx,idy+1],y_coor[idx,idy+1]
                counter += 1
        M_mat = np.array(M_mat)

        node_coordinates = self.nodecoor
        #
        # moment_plot = np.zeros([elem_connectivity.shape[0],2,2])
        #
        # node_coordinates = np.array(node_coordinates)
        # for ielem in range(elem_connectivity.shape[0]):
        #
        #     idx,idy = elem_connectivity[ielem,:]
        #
        #     moment_plot[ielem,0,:] = node_coordinates[idx,:]
        #     moment_plot[ielem,1,:] = node_coordinates[idy,:]

        #plt.figure(figsize=figsize)
        fig, axs = plt.subplots(1, 1)

        if figsize is not None:
            fig.set_size_inches(figsize[0],figsize[1])
        # Create a continuous norm to map from data points to colors
        norm = plt.Normalize(M_mat.min()-1E-3, M_mat.max()+1E-3)
        lc = LineCollection(xy, cmap='viridis', norm=norm)
        # Set the values used for colormapping
        lc.set_array(M_mat)
        lc.set_linewidth(5)
        line = axs.add_collection(lc)
        fig.colorbar(line, ax=axs,label=label_axis)
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        miny = node_coordinates[:,1].min()
        minx = node_coordinates[:,0].min()
        maxy = node_coordinates[:,1].max()
        maxx = node_coordinates[:,0].max()
        dy = np.abs(maxy - miny) + 1E-6
        dx = np.abs(maxx - minx) + 1E-6
        plt.ylim([miny - dy*0.1,maxy + dy*0.1])
        plt.xlim([minx - dx*0.1,maxx + dx*0.1])


    def _stiffness_local(self,
                        E,
                        I,
                        A,
                        L):
        """
            The local stiffness matrix.
        """
        EI = E*I

        k_bar = E*A/L

        k_local = np.zeros([6,6])
        k_local[0,0] =  k_bar
        k_local[0,3] = -k_bar
        k_local[3,3] =  k_bar
        k_local[3,0] = -k_bar

        k_local[1,[1,2,4,5]] = (EI/L**3)*np.array([12,6*L,-12,6*L])
        k_local[2,[1,2,4,5]] = (EI/L**3)*np.array([6*L,4*L**2,-6*L,2*L**2])
        k_local[4,[1,2,4,5]] = (EI/L**3)*np.array([-12,-6*L,12,-6*L])
        k_local[5,[1,2,4,5]] = (EI/L**3)*np.array([6*L,2*L**2,-6*L,4*L**2])

        return k_local

    def _mass_local(self,density,A,L):
        """
            The local mass matrix.
        """
        rho = float(density)
        k_bar = rho * A * L/6
        k_local = np.zeros([6,6])
        k_local[0,0] =  2 * k_bar
        k_local[0,3] = 1 * k_bar
        k_local[3,3] =  2 * k_bar
        k_local[3,0] = 1 * k_bar

        k_local[1,[1,2,4,5]] = (rho*A*L/420.0)*np.array([156,22*L,54,-13*L])
        k_local[2,[1,2,4,5]] = (rho*A*L/420.0)*np.array([22*L,4*L**2,13*L,-3*L**2])
        k_local[4,[1,2,4,5]] = (rho*A*L/420.0)*np.array([54,13*L,156,-22*L])
        k_local[5,[1,2,4,5]] = (rho*A*L/420.0)*np.array([-13*L,-3*L**2,-22*L,4*L**2])
        return k_local

    def stiffness_global(self):
        """
        The global stiffness matrix of the system under consideration.

        inputs
        --------
            None

        returns
        --------

            matrix: array_like
                stiffness matrix associated with the whole system. The global degrees-of-freedom are returned.
        """
        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        nelem = self.nelem

        elem_E = self.elemmatprop["E"]
        elem_I = self.elemmatprop["I"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)

        K = np.zeros([self.nnodes*self.ndof,self.nnodes*self.ndof])
        dofidx_node = np.arange(self.ndof)
        for ielem in range(nelem):

            nodeidx_i,nodeidx_j = elemconn[ielem,:]

            dofidx_i = dofidx_node + 3 * nodeidx_i
            dofidx_j = dofidx_node + 3 * nodeidx_j

            coor_i = nodecoor[nodeidx_i,:]
            coor_j = nodecoor[nodeidx_j,:]

            vec_ji = coor_j - coor_i

            L = np.linalg.norm(vec_ji)

            E = elem_E[ielem]
            I = elem_I[ielem]
            A = elem_A[ielem]

            angle_rad = np.arctan2(vec_ji[1],vec_ji[0])

            T = self._rotation_matrix(angle_rad)
            k = self._stiffness_local(E,I,A,L)

            k = np.dot(T.T,np.dot(k,T))

            dofidx = (np.array([np.hstack([dofidx_i,dofidx_j])]).T).astype(int)

            K[dofidx,dofidx.T] += k

        return K


    def mass_global(self):
        """
        The global mass matrix of the system under consideration.

        inputs
        --------

        None

        returns
        --------

        matrix: array_like
            The mass matrix associated with the whole system. The global degrees-of-freedom are returned.
        """
        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        nelem = self.nelem
        elem_E = self.elemmatprop["E"]
        elem_I = self.elemmatprop["I"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)
        if elem_rho is None:
            raise ValueError("Density is not supplied in the input dictionary.")

        M = np.zeros([self.nnodes*self.ndof,self.nnodes*self.ndof])
        dofidx_node = np.arange(self.ndof)
        for ielem in range(nelem):

            nodeidx_i,nodeidx_j = elemconn[ielem,:]

            dofidx_i = dofidx_node + 3 * nodeidx_i
            dofidx_j = dofidx_node + 3 * nodeidx_j

            coor_i = nodecoor[nodeidx_i,:]
            coor_j = nodecoor[nodeidx_j,:]

            vec_ji = coor_j - coor_i

            L = np.linalg.norm(vec_ji)

            E   = elem_E[ielem]
            I   = elem_I[ielem]
            A   = elem_A[ielem]
            rho = elem_rho[ielem]

            angle_rad = np.arctan2(vec_ji[1],vec_ji[0])

            T = self._rotation_matrix(angle_rad)
            m = self._mass_local(rho,A,L)

            m = np.dot(T.T,np.dot(m,T))

            dofidx = (np.array([np.hstack([dofidx_i,dofidx_j])]).T).astype(int)

            M[dofidx,dofidx.T] += m

        return M

    def _force_global(self):
        """
            This calculates and returns the global force array
        """
        F = np.zeros(self.nnodes*self.ndof)
        prescribed_forces = self.prescribed_forces
        nforces = prescribed_forces.shape[0]
        for iforce in range(nforces):
            node,dof,value = prescribed_forces[iforce,:]
            nodedof = int(dof + 3 * node)
            F[nodedof] = value
        return F

    def solve(self):
        """
        Solve the system of equations for the prescribed displacements and loads.

        [K] {u} = {f}

        inputs
        ---------
        None

        returns
        ---------

        solution_dict: dictionary

            This dictionary contains the following keys:
                "displacements": array_like
                    The displacement solution for the system is given in this field.
                    The displacements of the nodes are presented. If we denote the data of the displacements as u, the notation is as follows:
                    u[0] = displacement of node0 in direction 0
                    u[1] = displacement of node0 in direction 1
                    u[2] = rotation of node0 in direction 2
                    u[3] = displacement of node1 in direction 0
                    u[4] = displacement of node1 in direction 1
                    u[5] = rotation of node1 in direction 2
                    ....

                "loads": array_like
                    The loads for the system is given in this field.
                    The loads applied to the nodes are presented here. This will contain applied loads and reaction forces. If we denote the data of the loads as f, the notation is as follows:
                    f[0] = force applied to node0 in direction 0
                    f[1] = force applied to node0 in direction 1
                    f[2] = moment applied to node0 in direction 2
                    f[3] = force applied to node1 in direction 0
                    f[4] = force applied to node1 in direction 1
                    f[5] = moment applied to node1 in direction 2
                    ....
        """
        if self._supplied_displacements == False:
            raise ValueError("Cannot solve the system. No displacements are specified.")
        if self._supplied_loads == False:
            raise ValueError("Cannot solve the system. No loads are specified.")

        K = self.stiffness_global()
        F = self._force_global()
        U = self._displacement_global()

        Kff,Kfp,Kpp = self._partition_matrix(K)
        Ff,_ = self._partition_vector(F)
        _,Up = self._partition_vector(U)

        Uf = np.linalg.solve(Kff,Ff - np.dot(Kfp,Up))
        Fp = np.dot(Kfp.T,Uf) + np.dot(Kpp,Up)

        U = self._make_global_vector(Uf,Up)
        F = self._make_global_vector(Ff,Fp)

        self.solution_dict = {"displacements": U,
                              "loads": F}

        self._solved_system_bool = True
        return self.solution_dict

    def post_print_nodal_data(self,nodal_data):
        """
            Print the nodal data using data from the global degrees-of-freedom.

            inputs:
            ---------

            nodal_data: array_like
                The global degrees-of-freedom_data.

            returns:
            ----------
            None
        """
        K = self.stiffness_global()
        nodal_data = nodal_data.squeeze()
        assert (K.shape[0]) == (nodal_data.shape[0]), "The nodal data does not have the expected shape of {}".format(K.shape[0])


        for idx in range(int(len(nodal_data)/3)):
            ux = nodal_data[3*idx]
            uy = nodal_data[3*idx+1]
            angle = nodal_data[3*idx+2]

            print(" "*20,"Node",idx," "*20)
            print("x: {:f}, y: {:f}, angle: {:f}".format(ux,uy,angle))

    # def post_print_nodal_data(self,nodal_data):
    #
    #     for idx in range(int(len(nodal_data)/2)):
    #         ux = nodal_data[2*idx]
    #         uy = nodal_data[2*idx+1]
    #
    #         print(" "*20,"Node",idx," "*20)
    #         print("x: {:f}, y: {:f}".format(ux,uy))

    def make_nodal_data_global(self,nodal_data_free,nodal_data_prescribed):
        """
            Convert data associated with the free degrees-of-freedom and prescribed degrees-of-freedom to global degrees-of-freedom.

            inputs:
            ---------
                nodal_data_free: array_like
                    The data associated with the free degrees-of-freedom.

                nodal_data_prescribed: array_like
                    The data associated with the prescribed degrees-of-freedom.

            returns:
            ---------
                nodal_data_global: array_like
                    The data associated with the global degrees-of-freedom.
        """

        K = self.stiffness_global()
        idx_f,idx_p = self._partition_vector(np.arange(K.shape[0]))

        data = np.zeros(K.shape[0])

        if (idx_f.size) != (nodal_data_free.size):
            raise ValueError("The free degree of freedom data should have the shape ({},)".format(idx_f.size))
        if (idx_p.size) != (nodal_data_prescribed.size):
            raise ValueError("The prescribed degree of freedom data should have the shape ({},)".format(idx_p.size))

        data[idx_f] = nodal_data_free
        data[idx_p] = nodal_data_prescribed

        return data

    def post_print_element_data(self,array):
        if isinstance(array,(np.ndarray,list)) == False:
            raise ValueError("The input data should be an array.")
        array = list(array)
        for idx in range(len(array)):
            print(" "*20,"Element",idx," "*20)
            print("Data: {:f}".format(array[idx]/1E6))

    def post_print_solution(self,displacement=None):
        """
        This function prints the displacements, forces and stresses. The units are consistent with the units supplied in the input of the class.

        inputs
        ---------
            None

        returns
        ---------
            None
        """
        if displacement is None:
            if self._solved_system_bool == False:
                raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
            u = self.solution_dict["displacements"]
            F = self.solution_dict["loads"]
        else:
            u = np.array(displacement)
            K = self.stiffness_global()
            if K.shape[0] != u.shape[0]:
                raise ValueError("The displacement input is the wrong shape.")
            F = np.dot(K,u)

        # if self._solved_system_bool == False:
        #     raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        #
        # u = self.solution_dict["displacements"]
        # F = self.solution_dict["loads"]

        for idx in range(int(len(u)/3)):
            ux = u[3*idx]
            uy = u[3*idx+1]
            angle = u[3*idx+2]
            Fx = F[3*idx]
            Fy = F[3*idx+1]
            Mz = F[3*idx+2]

            print(" "*20,"Node",idx," "*20)
            print("u_x: {:f} x 10^-3, u_y: {:f} x 10^-3, angle: {:f}, F_x: {:f}, F_y: {:f}, M_z: {:f}".format(ux*1000,uy*1000,angle,Fx,Fy,Mz))

        stress = self.post_axial_stress_strain()["stress"]
        #strain = dict_post["strain"]

        print("")
        print("*"*20,"WARNING","*"*20)
        print("*** Bending stresses not calculated ***")
        print("*"*20,"WARNING","*"*20)

        print("")
        for idx in range(len(stress)):
            print(" "*20,"Element",idx," "*20)
            print("Normal stress due to axial load (only): {:f} x 10^6".format(stress[idx]/1E6))

        print("")
        print("*"*20,"WARNING","*"*20)
        print("*** Bending stresses not calculated ***")
        print("*"*20,"WARNING","*"*20)

    def post_get_displacement_element(self,element_number,npoints=20,displacement=None):

        if (displacement is None):
            try:
                u = self.solution_dict["displacements"]
            except Exception as e:
                raise ValueError("The system is not solved yet. Please solve the system.")
        else:
            u = np.array(displacement).squeeze()

            Ks = self.stiffness_global().shape[0]
            assert u.size == Ks,"The displacement size should match the global dof of {}".format(self.stiffness_global().shape[0])

        dict_post_M = self._calc_displacement_and_bending_moment(u,N_points=npoints)

        assert (element_number < dict_post_M["x"].shape[0]), "The element number exceeds the number of elements."

        return {"x_displ": dict_post_M["displ_x"][element_number,:],
                "y_displ": dict_post_M["displ_y"][element_number,:],
                "x_coor": dict_post_M["x"][element_number],
                "y_coor": dict_post_M["y"][element_number]}

    def eigen(self,number_of_eigenvalues=10):
        """
        This performs an eigenvalue analysis.


        inputs:
        -----------

        number_of_eigenvalues: integer

            The number of modes that needs to be calculated.

        returns:
        -----------
        results: dictionary with the following keys

            "natural_frequencies_hz": Natural frequencies in Hz.
            "modes": The eigenvectors associated with the free-degrees of freedom.

            Not important:
            "coor_global": The global coordinates.
            "index_f": Indices associated with the f-dofs.
            "index_p": Indices associated with the p-dofs.

        data["natural_frequencies_hz"] = fn #The natural frequencies of the system are in Hz
        data["modes"] = V                #Raw modes
        data["coor_global"] = XY
        data["index_p"] = idx_p
        data["index_f"] = idx_f

        """
        K = self.stiffness_global()
        M = self.mass_global()

        Kff,Kfp,Kpp = self._partition_matrix(K)
        Mff,Mfp,Mpp = self._partition_matrix(M)

        idx_f,idx_p = self._partition_vector(np.arange(K.shape[0]))

        import scipy.linalg as scilin
        number_of_eigenvalues = np.min([number_of_eigenvalues,Kff.shape[0]-1])
        eigvals = (0,number_of_eigenvalues)

        try:
            D,V = scilin.eigh(
                Kff,
                b=Mff,
                lower=True,
                eigvals_only=False,
                overwrite_a=False,
                overwrite_b=False,
                turbo=True,
                eigvals=None,
                type=1,
                check_finite=True,
                subset_by_index=eigvals,
                subset_by_value=None,
                driver=None,
            )
        except Exception as e:
            print("*"*100)
            print("The eigenvalue decomposition does not work. This might be due to an older version of scipy.linalg being used.")
            print("The error message that was received is:",e)
            print("Trying a simplified eigenvalue decomposition....")
            print("*"*100)
            D,V = scilin.eigh(
                Kff,
                b=Mff,
                lower=True
            )
            print("It seems to be working.")



        fn = np.sqrt(np.abs(D)) / (2*np.pi)

        data = {}

        y = self.nodecoor[:,1]
        x = self.nodecoor[:,0]

        XY = np.hstack([x,y])

        data["natural_frequencies_hz"] = fn #The natural frequencies of the system are in Hz
        data["modes"] = V                #Raw modes
        data["coor_global"] = XY
        data["index_p"] = idx_p
        data["index_f"] = idx_f
        return data



        #
        # x = self.nodecoor[:,0]
        # y = self.nodecoor[:,1]
        #
        # X = np.zeros(self.ndof * len(x))
        # Y = np.zeros(self.ndof * len(x))
        # X[::3]  = x
        # Y[1::3] = y
        #
        # coor = np.zeros(X.shape)
        #
        # coor[::3] = x
        # coor[1::3] = y
        #
        # direction = 2 * np.ones(self.ndof * len(x))
        # direction[0::3] = 0
        # direction[1::3] = 1
        # d,_ = self.partition_vector(direction)
        # x,_ = self.partition_vector(X)
        # y,_ = self.partition_vector(Y)
        # x = x[d == 0]
        # y = y[d == 1]
        #
        # K = self.stiffness_global()
        # M = self.mass_global()
        #
        # Kff,Kfp,Kpp = self.partition_matrix(K)
        # Mff,Mfp,Mpp = self.partition_matrix(M)
        #
        # import scipy.linalg as scilin
        # eigvals = (0,number_of_eigenvalues)
        #
        # D,V = scilin.eigh(
        #     Kff,
        #     b=Mff,
        #     lower=True,
        #     eigvals_only=False,
        #     overwrite_a=False,
        #     overwrite_b=False,
        #     turbo=True,
        #     eigvals=None,
        #     type=1,
        #     check_finite=True,
        #     subset_by_index=eigvals,
        #     subset_by_value=None,
        #     driver=None,
        # )
        #
        # fn = np.abs(D)**0.5 / (2*np.pi)
        #
        # data = {}
        #
        # data["natural_frequencies"] = fn
        # data["modes"] = V
        # data["x"] = x
        # data["y"] = y
        # data["coor"] = coor
        # data["directions"] = d
        # return data

    # def plot_modes(self,nmodes=3,tol = 1E-8):
    #
    #     print("Not working as expected. Needs to be checked.")
    #
    #     raise ValueError("Skip")
    #
    #     data = self.eigen()
    #     fn = data["natural_frequencies"]
    #
    #     print(fn[:10])
    #
    #     coor = data["coor"]
    #     V = data["modes"]
    #     x = data["x"]
    #     y = data["y"]
    #     directions = data["directions"]
    #
    #     print(x)
    #     print(y)
    #
    #     V = V * (V > tol)
    #
    #     fignum = np.zeros((nmodes))
    #     scaling_vector = np.linspace(-1,1,5)
    #     for idx in range(len(scaling_vector)):
    #
    #         scaling = scaling_vector[idx]
    #
    #         # mx = (V[::3,:].T * scaling+x).T
    #         # my = (V[1::3,:].T * scaling+y).T
    #
    #         fdof,pdof = self.get_fdof_pdof()
    #         for modenum in range(nmodes):
    #             if idx == 0:
    #                 plt.figure()
    #                 fignum[modenum] = plt.gcf().number
    #
    #             dvx = V[directions == 0,modenum]
    #             dvy = V[directions == 1,modenum]
    #
    #             x_plot = dvx * scaling + x
    #             y_plot = dvy * scaling + y
    #             coor[fdof[directions == 0]] = x_plot
    #             coor[fdof[directions == 1]] = y_plot
    #
    #             plt.figure(fignum[modenum])
    #             plt.plot(coor[::3],coor[1::3])
    #             plt.title("$f_{n}$" + "$ = {}$ Hz".format(np.round(fn[modenum],4)))

    def _displacement_global(self):

        u = np.zeros(self.nnodes*self.ndof)
        prescribed_displacement = self.prescribed_displacement
        ndisp = prescribed_displacement.shape[0]

        prescribed_dof = []
        for iforce in range(ndisp):

            node,dof,value = prescribed_displacement[iforce,:]

            nodedof = int(dof + 3 * node)

            u[nodedof] = value

            prescribed_dof.append(nodedof)

        self.prescribed_dof = prescribed_dof
        return u

    def _partition_vector(self,B):
        fdof,pdof = self._get_fdof_pdof()

        pdof = np.array(pdof)
        if len(pdof.squeeze()) > 0:
            return B[fdof],B[pdof]
        else:
            return B,np.array([])

    def _get_fdof_pdof(self):
        #
        # pdof = np.array(self.prescribed_dof)
        # fdof = np.array(list(set(list(np.arange(self.nnodes*self.ndof))) - set(list(pdof))))
        #
        # return fdof,pdof
        prescribed_displacement = self.prescribed_displacement

        if prescribed_displacement is not None:
            ndisp = prescribed_displacement.shape[0]
        else:
            ndisp = 0
        prescribed_dof = []
        for iforce in range(ndisp):

            node,dof,_ = prescribed_displacement[iforce,:]
            nodedof = int(dof + self.ndof * node)
            prescribed_dof.append(nodedof)

        self.prescribed_dof = prescribed_dof
        pdof = np.array(self.prescribed_dof)
        fdof = np.array(list(set(list(np.arange(self.nnodes*self.ndof))) - set(list(pdof))))

        #print(fdof),print(pdof)

        return fdof,pdof


    def _partition_matrix(self,A):


        fdof,pdof = self._get_fdof_pdof()
        fdof = np.array([fdof])
        pdof = np.array([pdof])

        if len(pdof.squeeze()) > 0:
            App = A[pdof.T,pdof]
            Afp = A[pdof.T,fdof].T
        else:
            App = np.array([])
            Afp = np.array([])
        Aff = A[fdof.T,fdof]

        return Aff,Afp,App

    def _make_global_vector(self,Uf,Up=None):
        fdof,pdof = self._get_fdof_pdof()
        if Up is None:
            Up = np.zeros(pdof.shape)


        ndof_tot = self.nnodes*self.ndof
        U = np.zeros(ndof_tot)
        U[fdof] = Uf
        U[pdof] = Up
        return U

    def dynamic_make_global_displacement(self,Uf):
        """

        This function transforms the displacements associated with the free degrees-of-freedom into the global displacement vector.

        inputs
        ---------

        displacement_free: array_like
            This is the displacement associated with the free-degrees-of-freedom. The order should be consistent with the notation used in this code.

        WARNING: It is assumed that the prescribed displacements are zero.

        returns
        ----------
        global displacement: array_like
            This is the global displacement vector. If we denote the displacements as u, the notation is as follows:
                    u[0] = displacement of node0 in direction 0
                    u[1] = displacement of node0 in direction 1
                    u[2] = displacement of node0 in direction 2 (rotation)
                    ...

        """
        return self._make_global_vector(Uf)

    def make_force_global(self,node,direction,force_magnitude):
        """
        inputs
        -------

        node: integer
            The node where the foce will be applied.

        direction: integer
            The direction in which the force will be applied.

        force_magnitude: integer
            The magnitude of the force.

        Returns
        ---------

        F: array_like
            The global force array that has the same dimensions as the global degrees-of-freedom of the structure

        Ff: array_like
            The force associated with the free degrees-of-freedom. This is used in dynamic analyses.

        Fp: array_like
            The force associated with prescribed degrees-of-freedom. This contains the reaction forces.
        """
        self._checks_direction(direction)
        self._checks_node_number(node)
        #self._checks_element_number(element)

        F = np.zeros(self.nnodes*self.ndof)
        # prescribed_forces = self.prescribed_forces
        # nforces = prescribed_forces.shape[0]
        # for iforce in range(nforces):
        node,dof,value = node,direction,force_magnitude#prescribed_forces[iforce,:]
        nodedof = int(dof + 3 * node)
        F[nodedof] = value

        fdof,pdof = self._get_fdof_pdof()
        return F,F[fdof],F[pdof]
        # return F


    def post_axial_stress_strain(self,displacement=None):
        """
            This function calculates the axial stress for all elements.

            This does not include the stress due to bending moments.
        """

        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        nelem = self.nelem

        if displacement is None:
            u = np.array(self.solution_dict["displacements"])
        else:
            u = np.array(displacement)

        assert u.shape[0] == self.stiffness_global().shape[0],"The size of the displacement vector is wrong. It needs to have a length of {}. Make sure that you are supplying the global DOF".format(self.stiffness_global().shape[0])

        elem_E = self.elemmatprop["E"]
        elem_I = self.elemmatprop["I"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)
        dofidx_node = np.arange(self.ndof)
        elem_stress = np.zeros(len(elem_E))
        elem_strain = np.zeros(len(elem_E))
        for ielem in range(len(elem_E)):

            nodeidx_i,nodeidx_j = elemconn[ielem,:]

            dofidx_i = dofidx_node + 3 * nodeidx_i
            dofidx_j = dofidx_node + 3 * nodeidx_j

            coor_i = nodecoor[nodeidx_i,:]
            coor_j = nodecoor[nodeidx_j,:]

            vec_ji = coor_j - coor_i

            L = np.linalg.norm(vec_ji)
            E = elem_E[ielem]
            I = elem_I[ielem]

            angle = np.arctan2(vec_ji[1],vec_ji[0])

            uj = np.array([u[3*nodeidx_j],u[3*nodeidx_j+1]])
            ui = np.array([u[3*nodeidx_i],u[3*nodeidx_i+1]])
            #-------------------------------------------------------------------
            # Small strain assumption:
            disp = uj - ui #Displacement of the vector in the global space
            unitv = (coor_j - coor_i)/np.linalg.norm(coor_j - coor_i) #Directional vector of element
            disp = np.dot(unitv,disp) #Displacement of element that is related to strain
            #-------------------------------------------------------------------

            # T = self._rotation_matrix(angle)
            #
            # index = np.hstack([ nodeidx_i * self.ndof + np.array([0,1,2]),nodeidx_j * self.ndof + np.array([0,1,2])])
            # uv = u[index]
            #
            # u_new = np.dot(T,uv)
            #
            elem_strain[ielem] = disp/L
            elem_stress[ielem] = elem_E[ielem] * elem_strain[ielem]

        return {"stress": elem_stress,"strain": elem_strain}

    def plot_modes(self,magnification=1,
                        number_of_modes=3,
                        show_node_numbers=True,
                        show_element_numbers=True,
                        npoints=100):

         # show_node_numbers=True,
         # show_displacement=False,
         # print_text=None,
         # scaling_text=1,
         # show_element_numbers=True,
        """
            This function plots the modes.

            inputs:
            ---------
            magnification: integer,float
                Magnifies the displacement of the modes with this factor.

            number_of_modes: integer
                NUmber of modes that need to be calculated.

            returns:
            -----------
            None
        """

        eigen_output = self.eigen(number_of_eigenvalues=number_of_modes+1)
        #
        # modes = eigen_output["modes"]
        # Xg    = eigen_output["coor_global"]
        # idx_p = eigen_output["index_p"]
        # idx_f = np.array(eigen_output["index_f"])
        #
        # print(np.vstack([idx_f,idx_f/2]))
        #
        # dum_va = np.abs(idx_f+1)/3 # Isolate moments - they will be integers
        # idx = np.where(np.abs(np.round(dum_va) - dum_va) > 1E-6)
        #
        # idx_fn = idx_f[idx]
        # #idx_f = idx_f[np.abs(idx_f - 2) != 0]

        #eigen_output= fembeam.eigen()
        modes = eigen_output["modes"]
        Xg    = eigen_output["coor_global"]
        idx_p = eigen_output["index_p"]
        idx_f = np.array(eigen_output["index_f"])

        #self = fembeam

        # print(XYZ.shape)
        # XYZ[::3] = self.nodecoor[:,0]
        # XYZ[1::3] = self.nodecoor[:,1]

        mode = modes[:,0]


        for eigennumber in range(np.min([modes.shape[1],number_of_modes])):
            V = modes[:,eigennumber]
            U = np.zeros(self.nodecoor.shape[0]*self.ndof)
            dum_va = np.abs(idx_f+1)/3 # Isolate moments - they will be integers
            idx = np.where(np.abs(np.round(dum_va) - dum_va) > 1E-6)
            U[idx_f] = V
            #U = U[idx]

            # Vg = np.zeros(Xg.shape[0]*self.ndof)
            # dum_va = np.abs(np.arange(Vg.shape[0])+1)/3
            # idx = np.where(np.abs(np.round(dum_va) - dum_va) > 1E-6)
            #
            # Vg[idx_f] = V
            # Vg = Vg[idx]

            self.plot_system(fignum=None,displacement_vector=U,magnification=magnification,show_node_numbers=show_node_numbers,show_element_numbers=show_element_numbers,npoints=npoints)
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.title(str(eigen_output["natural_frequencies_hz"][eigennumber]) + " Hz")

    def stiffness_global_partition(self):
        """
        The governing equations for a linear elastic solid mechanics problem can be written in the following form when using the finite element method:

            [K]{u} = {f}

        We can decompose the stiffness matrix in terms of ff, fp, pp, pf degrees of freedom as follows

        [K] = [[Kff], [Kfp]
               [Kpf], [Kpp]]

        Since the matrix is symmetric, Kfp = Kpf.T

        inputs
        -------
            None

        returns
        --------

            Kff: array_like

            Kfp: array_like

            Kpp: array_like

        """

        K = self.stiffness_global()

        Kff,Kfp,Kpp = self._partition_matrix(K)

        return Kff,Kfp,Kpp

    def mass_global_partition(self):
        """
        We can decompose the mass matrix in terms of ff, fp, pp, pf degrees of freedom as follows

        [M] = [[Mff], [Mfp]
               [Mpf], [Mpp]]

        Since the matrix is symmetric, Mfp = Mpf.T

        inputs
        -------

            None

        returns
        --------

            Mff: array_like

            Mfp: array_like

            Mpp: array_like
        """
        M = self.mass_global()

        Mff,Mfp,Mpp = self._partition_matrix(M)

        return Mff,Mfp,Mpp

    def dynamic_make_matrices(self):
        """
        In a dynamic analysis, the following governing equation is solved for an undamped system:

        Mff x'' + Kff x = Ff(t)

        This function returns Mff and Kff.

        inputs
        ----------

        None

        returns
        -----------

        Mff: array_like
             The mass matrix corresponding to the free degrees-of-freedom.

        Kff: array_like
             The stiffness matrix corresponding to the free degrees-of-freedom.
        """

        Mff,_,_ = self.mass_global_partition()
        Kff,_,_ = self.stiffness_global_partition()

        return Mff,Kff

    def _calc_displacement_and_bending_moment(self,displacement,N_points=200):
        """
            This function calculates the bending moment for the entire problem for a given displacement field.
        """

        # F = self.solution_dict["loads"]
        elemconn = self.elemconn
        nodecoor = self.nodecoor
        elemmatprop = self.elemmatprop
        nelem = self.nelem

        u = np.array(displacement)
        x = self.nodecoor[:,0]
        elem_E = self.elemmatprop["E"]
        elem_I = self.elemmatprop["I"]
        elem_A = self.elemmatprop["A"]
        elem_rho = self.elemmatprop.get("density",None)
        dofidx_node = np.arange(self.ndof)
        M = np.zeros([elem_E.shape[0],N_points])
        x_coor = np.zeros(M.shape)
        y_coor = np.zeros(M.shape)
        v      = np.zeros(M.shape)
        U      = np.zeros(M.shape)
        for ielem in range(len(elem_E)):

            nodeidx_i,nodeidx_j = elemconn[ielem,:]

            dofidx_i = dofidx_node + 3 * nodeidx_i
            dofidx_j = dofidx_node + 3 * nodeidx_j

            coor_i = nodecoor[nodeidx_i,:]
            coor_j = nodecoor[nodeidx_j,:]

            vec_ji = coor_j - coor_i

            L = np.linalg.norm(vec_ji)
            E = elem_E[ielem]
            I = elem_I[ielem]

            angle = np.arctan2(vec_ji[1],vec_ji[0])

            T = self._rotation_matrix(angle)
            # T = T[:,[1,2,4,5]]
            # T = T[[1,2,4,5],:]

            fraction_vector = np.linspace(0,1-1E-10,N_points)
            x_coor[ielem,:] = L * fraction_vector * np.cos(angle) + coor_i[0]
            y_coor[ielem,:] = L * fraction_vector * np.sin(angle) + coor_i[1]
            for idx_pos,fraction in enumerate(fraction_vector):

                x = L * fraction

                kappa = np.array([-6/L**2 + 12*x/L**3,
                                    -4/L + 6*x/L**2,
                                    6/L**2 - 12*x/L**3,
                                    -2/L + 6*x/L**2,
                                    ])

                index = np.hstack([ nodeidx_i * self.ndof + np.array([0,1,2]),nodeidx_j * self.ndof + np.array([0,1,2])])

                uv = np.dot(T,u[index])

                M[ielem,idx_pos] = E * I * np.dot(kappa,uv[[1,2,4,5]])

                v1,th1,v2,th2 = uv[[1,2,4,5]]

                U[ielem,idx_pos] = np.dot([1-x/L,x/L],[uv[0],uv[3]])
                v[ielem,idx_pos] = (2 / L**3.0 * (v1 - v2) + 1/L**2.0 * (th1 + th2)) * x**3 + (-3/L**2 * (v1 - v2) - 1/L * (2 * th1 + th2))*x**2 + th1*x + v1


            dumU = np.linalg.solve(T[:3,:3],np.vstack([U[ielem,:],v[ielem,:],v[ielem,:]*0]))

            # print("M.shape",M.shape,x_coor.shape)
            # plt.figure(1001)
            # plt.plot(x_coor[0,:],M[0,:],'b')
            # plt.plot(x_coor[1,:],M[1,:],'r')
            # plt.show()

            U[ielem,:] = dumU[0,:]
            v[ielem,:] = dumU[1,:]

        return {"moment": M,"x": x_coor, "y": y_coor,"displ_y": v,"displ_x": U}

    # def _displacement_interpolator(x):

    def _rotation_matrix(self,angle_rad):

        T = np.zeros([6,6])

        C = np.cos(angle_rad)
        S = np.sin(angle_rad)
        T[0,:] = [C,S,0,0,0,0]
        T[1,:] = [-S,C,0,0,0,0]
        T[2,:] = [0,0,1,0,0,0]
        T[3,:] = [0,0,0,C,S,0]
        T[4,:] = [0,0,0,-S,C,0]
        T[5,:] = [0,0,0,0,0,1]

        return T

    def post_get_displacement(self,node,direction):
        """
        Returns the displacement at a node.

        inputs:
        ---------

        node: integer
            The node number starting at 0

        direction: integer
            The direction (0 = x, 1 = y,2 = rotation).

        returns:
        ----------

        displacement: float
            The displacement of the node in the specified direction.
        """

        self._checks_direction(direction)
        self._checks_node_number(node)
        #self._checks_element_number(element)

        nnodes = self.nodecoor.shape[0]

        node = int(node)
        if (node >= nnodes) | (node <0 ):
            raise ValueError("The node does not exist")
        direction = int(direction)
        if direction not in [0,1,2]:
            raise ValueError("The direction does not exist.")

        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        return self.solution_dict["displacements"][int(node*self.ndof+direction)]

    def post_get_force(self,node,direction):
        """
            Returns the force at a node.

            inputs:
            ---------

            node: integer
                The node number starting at 0

            direction: integer
                The direction (0 = x, 1 = y, 2=moment).

            returns:
            ----------

            force: float
                The force applied to the node in the specified direction.

        """

        self._checks_direction(direction)
        self._checks_node_number(node)
        #self._checks_element_number(element)

        nnodes = self.nodecoor.shape[0]

        node = int(node)
        if (node >= nnodes) | (node <0 ):
            raise ValueError("The node does not exist")
        direction = int(direction)
        if direction not in [0,1,2]:
            raise ValueError("The direction does not exist.")

        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        return self.solution_dict["loads"][int(node*self.ndof+direction)]

    def post_get_stress(self,element_number):
        """
        This function gets the stress of an element.

        input:
        --------

        element_number: integer
            The element number of interest. Numbering starts at 0.

        returns:
        --------

        stress: float
            The stress in the specified element.

        WARNING:
        It does not include bending moments.

        """

        # self._checks_direction(direction)
        # self._checks_node_number(node)
        self._checks_element_number(element_number)

        if self._solved_system_bool == False:
            raise ValueError("The system is not solved yet. Run .solve() before applying this function.")
        dict_ss = self.post_axial_stress_strain()
        print("WARNING: It does not include the stresses due to the bending moments.")
        S = dict_ss["stress"]

        element_number = int(element_number)
        if element_number >= self.elemconn.shape[0]:
            raise ValueError("The specified element number does not exist.")

        return S[element_number]
