from typing import List, Tuple
import math


def _t_cdf(t: float, df: float) -> float:
    if df <= 0:
        raise ValueError("Degrees of freedom must be positive")
    
    if df == 1:
        return 0.5 + math.atan(t) / math.pi
    
    x = t / math.sqrt(df)
    a = df / 2.0
    b = 0.5
    z = df / (df + t * t)
    
    beta_val = _betai(a, b, z)
    
    if t >= 0:
        return 1.0 - 0.5 * beta_val
    else:
        return 0.5 * beta_val


def _gammln(xx: float) -> float:
    cof = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ]
    
    x = xx - 1.0
    tmp = x + 5.5
    tmp -= (x + 0.5) * math.log(tmp)
    ser = 1.000000000190015
    
    for j in range(6):
        x += 1.0
        ser += cof[j] / x
    
    return -tmp + math.log(2.5066282746310005 * ser)


def _betacf(a: float, b: float, x: float) -> float:
    MAXIT = 200
    EPS = 3.0e-7
    FPMIN = 1.0e-30
    
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        
        if abs(delta - 1.0) < EPS:
            break
    
    return h


def _betai(a: float, b: float, x: float) -> float:
    if x < 0.0 or x > 1.0:
        raise ValueError("x out of range in betai")
    
    if x == 0.0 or x == 1.0:
        return x
    
    bt = math.exp(_gammln(a + b) - _gammln(a) - _gammln(b)
                  + a * math.log(x) + b * math.log(1.0 - x))
    
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    else:
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_test_two_tailed(t_stat: float, df: int) -> float:
    if df <= 0:
        return 1.0
    if t_stat == 0:
        return 1.0
    return 2.0 * (1.0 - _t_cdf(abs(t_stat), float(df)))


def linear_regression(
    x_values: List[float],
    y_values: List[float],
) -> Tuple[float, float, float, float]:
    n = len(x_values)
    if n < 3:
        raise ValueError("At least 3 data points required for linear regression")
    if len(y_values) != n:
        raise ValueError("x and y must have same length")
    
    sum_x = 0.0
    sum_y = 0.0
    sum_xy = 0.0
    sum_x2 = 0.0
    sum_y2 = 0.0
    
    for i in range(n):
        x = x_values[i]
        y = y_values[i]
        sum_x += x
        sum_y += y
        sum_xy += x * y
        sum_x2 += x * x
        sum_y2 += y * y
    
    mean_x = sum_x / n
    mean_y = sum_y / n
    
    ss_xy = sum_xy - n * mean_x * mean_y
    ss_xx = sum_x2 - n * mean_x * mean_x
    ss_yy = sum_y2 - n * mean_y * mean_y
    
    if ss_xx == 0.0:
        return 0.0, mean_y, 0.0, 1.0
    
    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x
    
    y_pred = [slope * x + intercept for x in x_values]
    ss_res = sum((y_values[i] - y_pred[i]) ** 2 for i in range(n))
    ss_tot = ss_yy
    
    if ss_tot == 0.0:
        r_squared = 1.0
        p_value = 1.0
    else:
        r_squared = 1.0 - ss_res / ss_tot
        df = n - 2
        if df <= 0 or ss_xx == 0.0:
            p_value = 1.0
        else:
            mse = ss_res / df
            se_slope = math.sqrt(mse / ss_xx)
            if se_slope == 0.0:
                p_value = 0.0 if slope != 0.0 else 1.0
            else:
                t_stat = slope / se_slope
                p_value = t_test_two_tailed(t_stat, df)
    
    return slope, intercept, r_squared, p_value


def determine_trend(slope: float, p_value: float, p_threshold: float = 0.05) -> str:
    if p_value >= p_threshold:
        return "stable"
    if abs(slope) < 1e-15:
        return "stable"
    if slope > 0:
        return "rising"
    else:
        return "falling"


def mean(values: List[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def median(values: List[float]) -> float:
    if not values:
        return 0.0
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    if n % 2 == 0:
        return (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2.0
    else:
        return sorted_vals[n // 2]
