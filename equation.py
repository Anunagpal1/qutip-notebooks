import csv
from math import cos, pi, sin, sqrt

from sympy.physics.wigner import wigner_3j as w3j
from sympy.physics.wigner import wigner_6j

wigner_3j = lambda *args: w3j(*args).n()

# from wigners import clebsch_gordan, wigner_3j

HBAR = 1.0545718e-34
muB = 9.274009994e-24
gJ = -2.00231930436092
theta = pi / 2
# fac = lambda x: factorial(int(x))
B1 = 2e-3


def A(fa, fb, I, J):
    return (
        B1**2
        * gJ**2
        * muB**2
        * (2 * fa + 1)
        * (2 * fb + 1)
        * wigner_6j(J, fa, I, fb, J, 1, prec=20)
        * J
        * (J + 1)
        * (2 * J + 1)
    ) / HBAR**2


def omega(fb, fa, mb, ma, J, I):
    # print(f"ma-mb= {ma-mb}")
    if ma == -mb:
        # print(1)
        # print(f"Wig: {wigner_3j(fb, 1, fa, mb, 0 , ma)} {(fb, 1, fa, mb, 0 , ma)}")
        val = (
            (cos(theta) ** 2) * (wigner_3j(fb, 1, fa, mb, 0, ma) ** 2) * A(fa, fb, I, J)
        )
        # print(f"omega: {val}")
        return val
    elif (mb == ma + 1) or (mb == ma - 1):
        return (
            (0.5)
            * (sin(theta) ** 2)
            * (wigner_3j(fb, 1, fa, -mb, 0, ma) ** 2)
            * A(fa, fb, I, J)
        )
    else:
        return 0


def main(I, J):
    F = I + J
    runs = []  # Fix for float values loop
    for j1 in range(0, F + 1):
        j2 = 1
        for j3 in range(0, F + 1):
            if abs(j1 - j2) <= j3 <= j1 + j2:
                for m1 in range(-j1, j1 + 1):
                    m2 = 0
                    for m3 in range(-j3, j3 + 1):
                        if (m1 + m2 + m3) == 0:
                            ome = omega(j1, j3, m1, m3, J, I)
                            t = 1e-3  # seconds
                            Acap = 0.5  # hyperfine structure constant
                            delta = (
                                0.5 * Acap * (F * (F + 1) - I * (I + 1) - J * (J + 1))
                            )
                            # print(ome)
                            if ome:
                                prob = (ome / (ome + delta)) * sin(
                                    0.5 * t * sqrt(ome + delta)
                                )
                                # print(prob)
                            runs.append(
                                {
                                    "F": F,
                                    "j1": j1,
                                    "j2": j2,
                                    "j3": j3,
                                    "m1": m1,
                                    "m2": m2,
                                    "m3": m3,
                                    "omega^2": ome,
                                }
                            )
    with open("result.csv", "w", newline="") as output_file:
        dict_writer = csv.DictWriter(output_file, runs[0].keys())
        dict_writer.writeheader()
        dict_writer.writerows(runs)


main(2, 2)
