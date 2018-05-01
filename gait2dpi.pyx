import numpy as np
cimport numpy as np

cdef extern from "gait2dpi.h":

    cdef enum:
        NDOF = 9
        NSTICK = 10

    ctypedef struct param_struct:
        double TrunkMass
        double TrunkInertia
        double TrunkCMy
        double ThighMass
        double ThighInertia
        double ThighCMy
        double ThighLen
        double ShankMass
        double ShankInertia
        double ShankCMy
        double ShankLen
        double FootMass
        double FootInertia
        double FootCMx
        double FootCMy
        double ContactY
        double ContactHeelX
        double ContactToeX
        double ContactStiff
        double ContactDamp
        double ContactV0
        double ContactFric
        double gravity

    void gait2dpi(param_struct* par,
                   double q[NDOF],
                   double qd[NDOF],
                   double qdd[NDOF],
                   double Vsurface[2],
                   double QQ[NDOF],
                   double dQQdq[NDOF*NDOF],
                   double dQQdqd[NDOF*NDOF],
                   double dQQdqdd[NDOF*NDOF],
                   double GRF[6],
                   double stick[NSTICK*2])


def evaluate_autolev_rhs(np.ndarray[np.double_t, ndim=1, mode='c'] generalized_coordinates,
                         np.ndarray[np.double_t, ndim=1, mode='c'] generalized_speeds,
                         np.ndarray[np.double_t, ndim=1, mode='c'] generalized_accelerations,
                         np.ndarray[np.double_t, ndim=1, mode='c'] velocity_surface,
                         constants):
    """This function takes the current values of the coordinates, speeds,
    and accelerations and returns the required torques, derrivatives i.e. dQQdp, dQQdpd, dQQdpdd, GRF, and stick coordinates.

    [QQ, dQQdp, dQQdpd, dQQdpdd, grf, stick] = f(q, q', q'',vel_surf u)

    Parameters
    ----------
    generalized_coordinates : ndarray of floats, shape(9,)
        Trunk translation (x, y) and the joint angles.
    generalized_speeds : ndarray of floats, shape(9,)
        Trunk translation rate and the joint rates.
    generalized_accelerations : ndarray of floats, shape(9,)
        Trunk translation acceleration and the joint acceleration.
    velocity_surface : ndarray of floats, shape(2,)
        The velocity of surface perturbation, contains left and right
    constants : dictionary
        A dictionary of floats with keys that match the Autolev constants'
        names.

    Returns
    -------
    specified_quantities : ndarray, shape(9,)
        Required trunk horizontal and vertical force and the joint torques.
    dQQ_dq: ndarray, shape(81,)
        The derivative of QQ to q.
    dQQ_dqd: ndarray, shape(81,)
        The derivative of QQ to qd.
    dQQ_dqdd: ndarray, shape(81,)
        The derivative of QQ to qdd.
    ground_reaction_forces : ndarray, shape(6,)
        The right and left ground reaction forces.
    stick_figure_coordinates : ndarray, shape(20,)
        The x and y coordinates of the important points.

    Notes
    -----

    Generalized Coordinates

    q1: x hip translation wrt ground
    q2: y hip translation wrt ground
    q3: trunk z rotation wrt ground
    q4: right thigh z rotation wrt trunk
    q5: right shank z rotation wrt right thigh
    q6: right foot z rotation wrt right shank
    q7: left thigh z rotation wrt trunk
    q8: left shank z rotation wrt left thigh
    q9: left foot z rotation wrt left shank

    Surface velocity
    
    s1: surface velocity under right foot
    s2: surface velocity under left foot

    Specified outputs

    t1: x force applied to trunk mass center
    t2: y force applied to trunk mass center
    t3: torque between ground and trunk
    t4: torque between right thigh and trunk
    t5: torque between right thigh and right shank
    t6: torque between right foot and right shank
    t7: torque between left thigh and trunk
    t8: torque between left thigh and left shank
    t9: torque between left foot and left shank

    GRFs

    grf1: right horizontal
    grf2: right vertical
    grf3: right moment
    grf4: left horizontal
    grf5: left vertical
    grf6: left moment

    Sticks:

    stick1: trunkCM position in x corrodinate
    stick2: trunkCM position in y corrodinate
    stick3: hip position in x corrodinate
    stick4: hip position in y corrodinate
    stick5: Rknee position in x corrodinate
    stick6: Rknee position in y corrodinate
    stick7: Rankle position in x corrodinate
    stick8: Rankle position in y corrodinate
    stick9: Rheel position in x corrodinate
    stick10: Rheel position in y corrodinate
    stick11: Rtoe position in x corrodinate
    stick12: Rtoe position in y corrodinate
    stick13: Lknee position in x corrodinate
    stick14: Lknee position in y corrodinate
    stick15: Lankle position in x corrodinate
    stick16: Lankle position in y corrodinate
    stick17: Lheel position in x corrodinate
    stick18: Lheel position in y corrodinate
    stick19: Ltoe position in x corrodinate
    stick20: Ltoe position in y corrodinate

    """

    cdef param_struct p = param_struct(
        TrunkMass=constants['TrunkMass'],
        TrunkInertia=constants['TrunkInertia'],
        TrunkCMy=constants['TrunkCMy'],
        ThighMass=constants['ThighMass'],
        ThighInertia=constants['ThighInertia'],
        ThighCMy=constants['ThighCMy'],
        ThighLen=constants['ThighLen'],
        ShankMass=constants['ShankMass'],
        ShankInertia=constants['ShankInertia'],
        ShankCMy=constants['ShankCMy'],
        ShankLen=constants['ShankLen'],
        FootMass=constants['FootMass'],
        FootInertia=constants['FootInertia'],
        FootCMx=constants['FootCMx'],
        FootCMy=constants['FootCMy'],
        ContactY=constants['ContactY'],
        ContactHeelX=constants['ContactHeelX'],
        ContactToeX=constants['ContactToeX'],
        ContactStiff=constants['ContactStiff'],
        ContactDamp=constants['ContactDamp'],
        ContactV0=constants['ContactV0'],
        ContactFric=constants['ContactFric'],
	gravity=constants['gravity'])

    # TODO: Should allow the option to pass these in, instead of creating a
    # new array on each call to this function. It would be faster.
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] specified_quantities = np.zeros(9)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] jac_dQQ_dq = np.zeros(81)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] jac_dQQ_dqd = np.zeros(81)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] jac_dQQ_dqdd = np.zeros(81)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] ground_reaction_forces = np.zeros(6)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] stick_figure_coordinates = np.zeros(20)

    gait2dpi(&p,
              <double*> generalized_coordinates.data,
              <double*> generalized_speeds.data,
              <double*> generalized_accelerations.data,
              <double*>  velocity_surface.data,
              <double*> specified_quantities.data,
              <double*> jac_dQQ_dq.data,
              <double*> jac_dQQ_dqd.data,
              <double*> jac_dQQ_dqdd.data,
              <double*> ground_reaction_forces.data,
              <double*> stick_figure_coordinates.data)

    return specified_quantities, jac_dQQ_dq, jac_dQQ_dqd, jac_dQQ_dqdd, ground_reaction_forces, stick_figure_coordinates
