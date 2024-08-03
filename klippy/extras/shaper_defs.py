# Definitions of the supported input shapers
#
# Copyright (C) 2020-2023  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import collections, math

SHAPER_VIBRATION_REDUCTION=20.
DEFAULT_DAMPING_RATIO = 0.1

InputShaperCfg = collections.namedtuple(
        'InputShaperCfg', ('name', 'init_func', 'min_freq'))
InputSmootherCfg = collections.namedtuple(
        'InputSmootherCfg', ('name', 'init_func', 'min_freq'))

def get_none_shaper():
    return ([], [])

def get_zv_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)
    A = [1., K]
    T = [0., .5*t_d]
    return (A, T)

def get_zvd_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)
    A = [1., 2.*K, K**2]
    T = [0., .5*t_d, t_d]
    return (A, T)

def get_mzv_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-.75 * damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)

    a1 = 1. - 1. / math.sqrt(2.)
    a2 = (math.sqrt(2.) - 1.) * K
    a3 = a1 * K * K

    A = [a1, a2, a3]
    T = [0., .375*t_d, .75*t_d]
    return (A, T)

def get_ei_shaper(shaper_freq, damping_ratio):
    v_tol = 1. / SHAPER_VIBRATION_REDUCTION # vibration tolerance
    df = math.sqrt(1. - damping_ratio**2)
    t_d = 1. / (shaper_freq * df)
    dr = damping_ratio

    a1 = (0.24968 + 0.24961 * v_tol) + (( 0.80008 + 1.23328 * v_tol) +
                                        ( 0.49599 + 3.17316 * v_tol) * dr) * dr
    a3 = (0.25149 + 0.21474 * v_tol) + ((-0.83249 + 1.41498 * v_tol) +
                                        ( 0.85181 - 4.90094 * v_tol) * dr) * dr
    a2 = 1. - a1 - a3

    t2 = 0.4999 + ((( 0.46159 + 8.57843 * v_tol) * v_tol) +
                   (((4.26169 - 108.644 * v_tol) * v_tol) +
                    ((1.75601 + 336.989 * v_tol) * v_tol) * dr) * dr) * dr

    A = [a1, a2, a3]
    T = [0., t2 * t_d, t_d]
    return (A, T)

def _get_shaper_from_expansion_coeffs(shaper_freq, damping_ratio, t, a):
    tau = 1. / shaper_freq
    T = []
    A = []
    n = len(a)
    k = len(a[0])
    for i in range(n):
        u = t[i][k-1]
        v = a[i][k-1]
        for j in range(k-1):
            u = u * damping_ratio + t[i][k-j-2]
            v = v * damping_ratio + a[i][k-j-2]
        T.append(u * tau)
        A.append(v)
    return (A, T)

def get_2hump_ei_shaper(shaper_freq, damping_ratio):
    t = [[0., 0., 0., 0.],
         [0.49890,  0.16270, -0.54262, 6.16180],
         [0.99748,  0.18382, -1.58270, 8.17120],
         [1.49920, -0.09297, -0.28338, 1.85710]]
    a = [[0.16054,  0.76699,  2.26560, -1.22750],
         [0.33911,  0.45081, -2.58080,  1.73650],
         [0.34089, -0.61533, -0.68765,  0.42261],
         [0.15997, -0.60246,  1.00280, -0.93145]]
    return _get_shaper_from_expansion_coeffs(shaper_freq, damping_ratio, t, a)

def get_3hump_ei_shaper(shaper_freq, damping_ratio):
    t = [[0., 0., 0., 0.],
         [0.49974,  0.23834,  0.44559, 12.4720],
         [0.99849,  0.29808, -2.36460, 23.3990],
         [1.49870,  0.10306, -2.01390, 17.0320],
         [1.99960, -0.28231,  0.61536, 5.40450]]
    a = [[0.11275,  0.76632,  3.29160 -1.44380],
         [0.23698,  0.61164, -2.57850,  4.85220],
         [0.30008, -0.19062, -2.14560,  0.13744],
         [0.23775, -0.73297,  0.46885, -2.08650],
         [0.11244, -0.45439,  0.96382, -1.46000]]
    return _get_shaper_from_expansion_coeffs(shaper_freq, damping_ratio, t, a)

def get_shaper_offset(A, T):
    if not A:
        return 0.
    return sum([a * t for a, t in zip(A, T)]) / sum(A)

def init_smoother(coeffs, smooth_time, normalize_coeffs):
    if not normalize_coeffs:
        return (list(reversed(coeffs)), smooth_time)
    n = len(coeffs)
    int_t0 = 0.
    for i in range(0, n, 2):
        int_t0 += coeffs[n-i-1] / (2**i * (i+1))
    inv_int = 1. / int_t0
    c = [0] * n
    inv_t_sm = inv_t_sm_n = 1. / smooth_time
    for i in range(n-1,-1,-1):
        c[n-i-1] = coeffs[i] * inv_t_sm_n * inv_int
        inv_t_sm_n *= inv_t_sm
    return (c, smooth_time)

def get_none_smoother():
    return ([1.], 0.)

# Input smoother is a polynomial function w(t) which is convolved with the input
# trajectory x(t) as x_sm(T) = integral(x(t) * w(T-t) * dt, t=[T-T_sm...T+T_sm])
# to obtain a modified trajectory which has low residual vibrations in a certain
# frequency range. Similar to discrete input shapers, the amplitude of residual
# vibrations of a step input signal can be calculated as
# *) V_c = integral(exp(-zeta * omega * t) * cos(omega_d * t) *
#                   * w(T_sm/2 - t) * dt, t=[0...T_sm])
# *) V_s = integral(exp(-zeta * omega * t) * sin(omega_d * t) *
#                   * w(T_sm/2 - t) * dt, t=[0...T_sm])
# *) V = sqrt(V_c^2 + V_s^)
#
# The input smoothers defined below were calculated and optimized in the Maxima
# algebra system by calculating V(w, T_sm, omega, zeta) in a closed form for
# several (depending on smoother) omega_i, zeta_i pairs and optimizing T_sm,
# omega_i, and zeta_i in order to obtain a polynomial w(t) with a minimal T_sm
# such that w(t) >= 0 and V(omega) <= 0.05 in a certain frequency range.
#
# Additional constraints that were enforced for w(t) during optimization:
# *) integral(w(t) * dt, t=[-T_sm/2...T_sm/2]) = 1
# *) w(-T_sm/2) = w(T_sm/2) = 0
# The first condition ensures unity norm of w(t) and the second one makes the
# input smoothers usable as extruder smoothers (for the ease of synchronization
# of the extruder with the toolhead motion).

def _get_smoother_from_expansion_coeffs(damping_ratio, t, a):
    if damping_ratio is None:
        damping_ratio = DEFAULT_DAMPING_RATIO
    n = len(a)
    k = len(a[0])
    C = []
    for i in range(n):
        v = a[i][k-1]
        for j in range(k-1):
            v = v * damping_ratio + a[i][k-j-2]
        C.append(v)
    t_s = t[k-1]
    for j in range(k-1):
        t_s = t_s * damping_ratio + t[k-j-2]
    return t_s, C[::-1]

def get_zv_smoother(shaper_freq, damping_ratio=None, normalize_coeffs=True):
    t = [0.66471804,0.0039585139,0.33614293,0.2157273]
    a = [[0.1418285825,-0.4423501425,0.696114325,-1.3609939],
         [0.0,3.5635179,-3.83613985,-6.2783295],
         [-15.50842883,29.06586482,-12.3383868,112.6474131],
         [0.0,-206.8822741,52.8629899,152.492858],
         [452.7401905,-281.3808295,-2.555613,-1175.326975],
         [0.0,770.51281,-150.073722,-509.51816],
         [-1571.90293,688.77989,163.085324,2986.0529]]
    t_s, coeffs = _get_smoother_from_expansion_coeffs(damping_ratio, t, a)
    return init_smoother(coeffs, t_s / shaper_freq, normalize_coeffs)

def get_mzv_smoother(shaper_freq, damping_ratio=None, normalize_coeffs=True):
    t = [0.85165483,0.0183642417,0.39000898,0.23209004]
    a = [[2.239795725,-0.3957072675,-5.57438075,0.07195255],
         [0.0,29.10368225,-3.461668825,-22.41514],
         [-154.4864979,92.80041407,23.0961391,774.3910148],
         [0.0,-1362.168729,-35.59728045,1147.640935],
         [3870.466035,-3292.00459,7849.3652856,-36092.1938],
         [0.0,21476.8835,2071.594648,-21131.4185],
         [-43110.2746,39792.0705,-138379.754,483530.1615],
         [0.0,-134532.955,-15615.6118,139113.6955],
         [214526.885,-191572.8415,791628.475,-2450648.74],
         [0.0,274229.94,32481.346,-286062.83],
         [-378797.98,316954.83,-1455001.66,4177695.0]]
    t_s, coeffs = _get_smoother_from_expansion_coeffs(damping_ratio, t, a)
    return init_smoother(coeffs, t_s / shaper_freq, normalize_coeffs)

def get_ei_smoother(shaper_freq, damping_ratio=None, normalize_coeffs=True):
    t = [0.96753911,0.48592866,-7.2991788,18.2075661]
    a = [[2.4923607,1.1935848,-46.9202175,127.76906],
         [0.0,33.3731105,103.46339,-262.7533525],
         [-166.2319578,-317.4835492,5178.61787,-13957.64074],
         [0.0,-1914.663767,-2867.941885,6546.52041],
         [4043.65056,13078.06184,-149368.2255,386432.888],
         [0.0,30664.02555,11115.66415,4553.999],
         [-43879.103,-162723.6625,1672995.285,-4352601.97],
         [0.0,-184244.5385,56995.0041,-422243.748],
         [213582.6565,773861.825,-7739194.35,2.06293873e7],
         [0.0,360348.71,-248768.99,1264398.56],
         [-371055.41,-1248811.1,1.24707394e7,-3.4165302e7]]
    t_s, coeffs = _get_smoother_from_expansion_coeffs(damping_ratio, t, a)
    return init_smoother(coeffs, t_s / shaper_freq, normalize_coeffs)

def get_2hump_ei_smoother(shaper_freq, damping_ratio_unused=None,
                          normalize_coeffs=True):
    coeffs = [-22525.88434486782,2524.826047114184,10554.22832043971,
               -1051.778387878068,-1475.914693073253,121.2177946817349,
              57.95603221424528,-4.018706414213658,0.8375784787864095]
    return init_smoother(coeffs, 1.14875 / shaper_freq, normalize_coeffs)

def get_si_smoother(shaper_freq, damping_ratio_unused=None,
                    normalize_coeffs=True):
    coeffs = [-6186.76006449789,1206.747198930197,2579.985143622855,
              -476.8554763069169,-295.546608490564,52.69679971161049,
              4.234582468800491,-2.226157642004671,1.267781046297883]
    return init_smoother(coeffs, 1.245 / shaper_freq, normalize_coeffs)

def get_zvd_ei_smoother(shaper_freq, damping_ratio_unused=None,
                        normalize_coeffs=True):
    coeffs = [-18835.07746719777,1914.349309746547,8786.608981369287,
              -807.3061869131075,-1209.429748155012,96.48879052981883,
              43.1595785340444,-3.577268915175282,1.083220648523371]
    return init_smoother(coeffs, 1.475 / shaper_freq, normalize_coeffs)

def get_smoother_offset(C, t_sm, normalized=True):
    int_t0 = int_t1 = 0.
    if normalized:
        for i, c in enumerate(C):
            if i & 1:
                int_t1 += c * t_sm**(i+2) / (2**(i+1) * (i+2))
            else:
                int_t0 += c * t_sm**(i+1) / (2**i * (i+1))
    else:
        for i, c in enumerate(C):
            if i & 1:
                int_t1 += c / (2**(i+1) * (i+2))
            else:
                int_t0 += c / (2**i * (i+1))
        int_t1 *= t_sm
    return int_t1 / int_t0

# min_freq for each shaper is chosen to have projected max_accel ~= 1500
INPUT_SHAPERS = [
    InputShaperCfg('zv', get_zv_shaper, min_freq=21.),
    InputShaperCfg('mzv', get_mzv_shaper, min_freq=23.),
    InputShaperCfg('zvd', get_zvd_shaper, min_freq=29.),
    InputShaperCfg('ei', get_ei_shaper, min_freq=29.),
    InputShaperCfg('2hump_ei', get_2hump_ei_shaper, min_freq=39.),
    InputShaperCfg('3hump_ei', get_3hump_ei_shaper, min_freq=48.),
]

# min_freq for each smoother is chosen to have projected max_accel ~= 1000
INPUT_SMOOTHERS = [
    InputSmootherCfg('smooth_zv', get_zv_smoother, min_freq=18.),
    InputSmootherCfg('smooth_mzv', get_mzv_smoother, min_freq=20.),
    InputSmootherCfg('smooth_ei', get_ei_smoother, min_freq=21.),
    InputSmootherCfg('smooth_2hump_ei', get_2hump_ei_smoother, min_freq=21.5),
    InputSmootherCfg('smooth_zvd_ei', get_zvd_ei_smoother, min_freq=26.),
    InputSmootherCfg('smooth_si', get_si_smoother, min_freq=21.5),
]
