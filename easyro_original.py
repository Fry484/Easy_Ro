def easyro(time_from_fulton, Temp_from_fulton):
    # Time and C are the inputs
    # Fill these lists with the temperatures at known times, and Easy%Ro calculates the expected %Ro for that heating rate
    # Time is in years

    time1 = np.linspace(0.0, (age * 31536000000000.0), num=10)
    C1 = np.linspace(geotherm, (geotherm * (depth / 1000.0) + 20.0), num=10)

    time2 = np.array([])
    for i in time_from_fulton:
        time2 = np.append(time2, (i + time1[-1]))
    C2 = Temp_from_fulton

    time = np.concatenate((time1, time2))
    C = np.concatenate((C1, C2))



    pre = 10000000000000.0  # s^-1
    R = 0.0019872041
    a1 = 2.334733
    a2 = 0.250621
    b1 = 3.330675
    b2 = 1.681534

    weights = (
    0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06, 0.04, 0.04, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04, 0.03, 0.02, 0.02,
    0.01)
    E = (
    34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, 66.0, 68.0, 70.0,
    72.0)
    ER = (17109.47, 18115.90, 19122.34, 20128.78, 21135.22, 22141.66, 23148.10, 24154.54, 25160.98, 26167.42, 27173.86,
          28180.30, 29186.74, 30193.17, 31199.61, 32206.05, 33212.49, 34218.93, 35225.37, 36231.81)

    # 1 my = 31600000000000 seconds

    K = [i + 273 for i in C]

    tarray = np.zeros((len(C), 1))
    TC = np.zeros((len(C), 1))
    TK = np.zeros((len(C), 1))
    H = np.ones((len(C), 1))
    Ro = np.zeros((len(C), 1))
    delI = np.zeros((len(C), 20))
    I = np.zeros((len(C), 20))
    ERT = np.zeros((len(C), 20))
    f = np.zeros((len(C), 20))
    F = np.zeros((len(C), 1))

    for val, i in enumerate(time):
        if val == 0:
            for v, l in enumerate(ER):
                ERT[val, v] = (l / K[val])
            for n in range(20):
                I[val, n] = pre * K[val] * m.exp(-ERT[val, n]) * (
                            1 - (ERT[val, n] ** 2 + a1 * ERT[val, n] + a2) / (ERT[val, n] ** 2 + b1 * ERT[val, n] + b2))
            Ro[val] = m.exp(-1.6 + 3.7 * F[val])
        if val > 0:
            for q in range(len(time) - 1):
                if C[val] - C[(val - 1)] == 0.0:
                    H[val] = 0.00000001
                else:
                    H[val] = (C[val] - C[(val - 1)]) / (time[val] - time[(val - 1)])  # / 31600000000000.0
            for v, l in enumerate(ER):
                ERT[val, v] = l / K[val]
            for n in range(20):
                I[val, n] = pre * K[val] * m.exp(-ERT[val, n]) * (
                            1 - (ERT[val, n] ** 2 + a1 * ERT[val, n] + a2) / (ERT[val, n] ** 2 + b1 * ERT[val, n] + b2))
            for p in range(20):
                delI[val, p] = delI[(val - 1), p] + (I[val, p] - I[(val - 1), p]) / H[val]
            for r in range(20):
                if delI[val, r] < 1e-20:
                    f[val, r] = 0.0
                elif delI[val, r] > 220:
                    f[val, r] = weights[r]
                else:
                    f[val, r] = weights[r] * (1.0 - m.exp(-delI[val, r]))
            F[val] = sum(f[val])
            Ro[val] = m.exp(-1.6 + 3.7 * F[val])
    # return ((Ro[999]),(Ro[3999]),max(C))

    return (C, Ro)
    # plt.subplot(212)
    # plt.plot(time,Ro)
