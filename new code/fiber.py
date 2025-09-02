import c_fibers_obj_clean
import numpy as np

class Fiber:
    def __init__(self):
        print("initialized")

    def get_bishop_frame():
        # To Do: infer bishop frame from vector soup
        print('Get bishop frame')

def get_bishop(tangent, u_0 = None):
    '''
    Return the coordinates of the blobs.
    '''
    # Get location and orientation
    if u_0 is None:
        t_0 = tangent[0,:]
        u_0 = np.array([1.0,0.0,0.0])
        u_0 = u_0 - np.dot(u_0,t_0)*t_0
        u_0 /= np.linalg.norm(u_0)

    # Compute blobs coordinates
    u = 0*tangent 
    v = 0*tangent 

    u[0,:] = u_0
    v[0,:] = np.cross(tangent[0,:],u_0)

    for k in range((tangent.size // 3)-1):
        t_k = tangent[k,:]
        t_kp = tangent[k+1,:]
        cos_th = np.dot(t_k,t_kp)
        rot_x = np.cross(t_k,t_kp)
        
        u_k = u[k,:]
        
        u[k+1,:] = u_k + np.cross(rot_x,u_k) + (1.0/(1.0 + cos_th)) * np.cross(rot_x,np.cross(rot_x,u_k))
        v[k+1,:] = np.cross(t_kp,u[k+1,:])

    return u,v

if __name__ == "__main__":
    ### Set some physical parameters
    eta = 1.0 # viscosity
    L = 2.0 # fiber length

    ### Make an angled fiber
    T_base = np.array([0.0,1.0,0.0])

    ###############################################################################################################
    # I made Nlinks the primary control variable, but it's pretty straightforward to change things around
    Nlinks = 10
    ### Make array of tangent vectors and bishop frame
    tangents = np.tile(T_base, (Nlinks, 1))
    u,v = get_bishop(tangents)
    # set a,ds based on Nlinks
    ds = L/(1.0*Nlinks)
    a = ds/2.0
    ###############################################################################################################


    ### Make an Nfib x 1 array of fibers at z=0,5a,10a,...
    Nfibs = 10
    X_0 = np.array([0.0,0.0,0.0])

    ### Make big matricees of size (3*Nfibers x Nlinks) for all of the fiber tangent vectors/bishop frame
    ### In this example, all of the fibers have the same tangents, but differnt x_0 to arrange them in a grid
    All_Taus = tangents.T
    All_X = X_0
    All_U = u.T
    All_V = v.T
    for kx in range(1,Nfibs):
        All_Taus = np.vstack((All_Taus, tangents.T))
        X_1 = X_0 + np.array([0.0,0.0,kx*5.0*a])
        All_X = np.r_[All_X, X_1]
        All_U = np.vstack((All_U, u.T))
        All_V = np.vstack((All_V, v.T))

    ### Print the shapes of the fibers
    print(All_Taus)
    print(All_X)
    print(All_U.shape)
    print(All_V.shape)
    ### Print the toal number of fibers
    Nfibs = int(len(All_X)/3)
    print('Fibers: ' + str(Nfibs))


    ### Set the thermal energy and compute the bending stiffness to have 
    ### A prescribed persistance length
    kBT = 0.004142
    #######################################
    k_b = 4.0*L*kBT # bending stiffness
    ######################################


    ### Physical scale for the mobility
    M0 = 1.0/(6.0*np.pi*eta*a) #np.log(L/a)/(8.0*np.pi*eta*a)
    ### Bending timescale
    alpha =((ds**3)/M0)/k_b
    ### dimensionless timestep
    alpha_fact = 1.0 ##################### Change this to change the time step
    ### Physcial time step
    dt = alpha_fact*alpha
    ### Factor for implicit time stepping
    ########################################################################################
    ########################################################################################
    ########################################################################################
    # Use impl_c=0.5 for CN
    # Use impl_c=1.0 for BE
    impl_c = 1.0
    impl = impl_c*dt*k_b/(ds*ds*ds)
    ########################################################################################
    ########################################################################################
    ########################################################################################


    ### Print timestep and size of fiber
    print('dt = ' + str(dt))
    Nlk = np.shape(All_U)[1]
    Nblobs = Nlk+1
    print('Blobs: ' + str(Nblobs))


    ###########################################
    #### Sets the number of neighbor blobs ####
    #### to be used in the PC              ####
    ###########################################
    M_bands = 0 # NOTE: 0 is fastest in most cases
    ###########################################
    ###########################################


    ###########################################
    #### Set kBT to zero to do a           ####
    #### determinsitic simulation          ####
    ###########################################
    ###########################################
    #kBT = 0.0 # NOTE: kBT is used to define bending stiffness ealier which is why we set it to zero here
    ###########################################
    ###########################################

    ###########################################
    #### Set Clamp to 'True' for a fiber   ####
    #### Bound to it's initial position at ####
    #### One of it's ends.                 ####
    #### Set Clamp to 'False' for a fiber  ####
    #### That is free at both ends         ####
    ###########################################
    ###########################################
    Clamp = False
    ###########################################
    ###########################################

    ###########################################
    ########## Periodic Length   ##############
    ###########################################
    ###########################################
    Lperiodic = -1 # Currently not used unless PSE is specified for the mobility
    ###########################################
    ###########################################

    ###########################################
    # Domain int specifies the mobility product 
    # to be used in the simulation. Right now
    # the values are:
    # DomainInt=0: Batched RPY (1 fiber per batch)
    # DomainInt=1: Full RPY
    # DomainInt=2: Batched RPB (1 fiber per batch)
    # DomainInt=3: Full RPB
    DomainInt = 0
    ###########################################

    cf = c_fibers_obj_clean.CManyFibers()
    cf.setParameters(DomainInt, Nfibs, Nblobs, a, ds, dt, k_b, M0, impl_c, kBT, eta, Lperiodic, Clamp, T_base)
    pos = cf.multi_fiber_Pos(All_Taus,All_X)
    print(pos)