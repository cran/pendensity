D.m                package:pendensity                R Documentation

_C_a_l_c_u_l_a_t_i_n_g _t_h_e _p_e_n_a_l_t_y _m_a_t_r_i_x

_D_e_s_c_r_i_p_t_i_o_n:

     calculating the penalty matrix depending on the number of
     covariates 'p', the order of differences to be penalized 'm', the
     corresponding difference matrix 'L' of order 'm', the covariate
     matrix 'Z', the number of observations 'n' and the number of knots
     'K'.

_U_s_a_g_e:

     D.m(penden.env)

_A_r_g_u_m_e_n_t_s:

penden.env: Containing all information, environment of pendensity()

_D_e_t_a_i_l_s:

     The penalty matrix is calculated as

     D_m=(L^T otimes I_p) (I_{K-m} otimes frac{Z^T Z}{n}) (L otimes
     I_p)

     The needed values are saved in the environment.

_V_a_l_u_e:

     Returning the penalty matrix.

_A_u_t_h_o_r(_s):

     Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>

_R_e_f_e_r_e_n_c_e_s:

     Penalized Density Estimation, Kauermann G. and Schellhase C.
     (2009), to appear.

