import math, random


class Statistics:
    def __init__(self):
        self.oracle_call = 0


class ChernoffBoundEstimator:
    def estimate(self, cobs, Nshot, delta):
        diff = math.sqrt(math.log(2.0 / delta) * 12.0 / Nshot)
        c_min = max(-1, cobs - diff)
        c_max = min(1, cobs + diff)
        return c_min, c_max

    def error(self, Nshot, delta):
        return math.sqrt(math.log(2.0 / delta) * 12.0 / Nshot)


class GroverSampler:
    def sample(self, n_shot, m):
        """
        sample cos((2^2m + 2)\theta) where \theta = amplitude/4 by n_shot measurements
        """
        pass

    def get_oracle_calls(self):
        pass


class GroverSimulator(GroverSampler):
    def __init__(self, amplitude):
        self.theta = amplitude / 4
        self.statistics = Statistics()

    def sample(self, n_shot, m):
        n_11 = 0
        for i in range(n_shot):
            self.statistics.oracle_call = self.statistics.oracle_call + m
            v = self._measure(self._p_11(self._angle(m, self.theta)))
            if v == 1:
                n_11 = n_11 + 1
        c2m = 1 - (2.0 * n_11) / n_shot
        return c2m

    def get_oracle_calls(self):
        return self.statistics.oracle_call

    def _p_11(self, angle):
        return (1.0 - math.cos(angle)) / 2

    def _measure(self, one_probability):
        rand = random.random()
        if one_probability > rand:
            return 1
        return 0

    def _angle(cls, m, theta):
        return 2 * (2 * m + 1) * theta


class FasterAmplitudeEstimation:
    def estimate_amplitude(self, sampler, l, delta_c):
        theta, o, j_0 = self._estimate_theta(sampler, l, delta_c)
        return theta * 4, o, j_0

    def _nu(self, theta_min, theta_max, j0):
        return pow(2, j0 + 1) * (theta_min + theta_max) / 2

    def _integer(self, j, rho, theta):
        result = (math.pow(2, j + 1) + 2) * theta - rho + math.pi / 3
        return math.floor(result / (2 * math.pi))

    def _del_rho(self, nu, del_nu, Nsecond, delta):
        bound_estimator = ChernoffBoundEstimator()
        cm_error = bound_estimator.error(Nsecond, delta)
        sin_error = (math.sqrt(2 - 2 * math.cos(del_nu)) + math.fabs(cm_error * math.cos(nu)) + math.fabs(
            cm_error)) / math.sin(nu)
        return max(2 * math.fabs(cm_error) + 2 * math.fabs(sin_error), math.fabs(5 * cm_error))

    def _estimate_theta(self, sampler, l, delta_c):
        first_stage = True
        theta_min = 0
        theta_max = 0.252
        n_first = int(1944 * math.log(2 / delta_c))
        n_second = int(972 * math.log(2 / delta_c))
        bound_estimator = ChernoffBoundEstimator()
        j_0 = l
        nu = 0
        for j in range(1, l + 1):
            if first_stage:
                m = pow(2, j - 1)
                cm = sampler.sample(n_first, m)
                cm_min, cm_max = bound_estimator.estimate(cm, n_first, delta_c)
                theta_min = math.acos(cm_max) / (pow(2, j + 1) + 2)
                theta_max = math.acos(cm_min) / (pow(2, j + 1) + 2)
                if pow(2, j + 1) * theta_max >= 3 * math.pi / 8 and j < l:
                    j_0 = j
                    nu = self._nu(theta_min, theta_max, j_0)
                    first_stage = False
            else:
                m = pow(2, j - 1)
                mprime = pow(2, j - 1) + pow(2, j_0 - 1)
                cm = sampler.sample(n_second, m)
                cmprime = sampler.sample(n_second, mprime)
                sm = (cm * math.cos(nu) - cmprime) / math.sin(nu)
                rho = math.atan2(sm, cm)
                n = self._integer(j, rho, theta_max)

                theta_min = (2 * math.pi * n + rho - math.pi / 3) / (pow(2, j + 1) + 2)
                theta_max = (2 * math.pi * n + rho + math.pi / 3) / (pow(2, j + 1) + 2)
        return (theta_min + theta_max) / 2, sampler.statistics.oracle_call, j_0
