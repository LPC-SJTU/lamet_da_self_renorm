# %%

import gvar as gv
import numpy as np
import matplotlib.pyplot as plt

fig_width = 6.75  # in inches, 2x as wide as APS column
fig_size_sq = (fig_width * 0.8, fig_width * 0.8)
plt_axes_small = [0.17, 0.15, 0.8, 0.75]
fs_p_l = {"fontsize": 19}
ls_p_l = {"labelsize": 15.5}


color_list = ["orange", "dodgerblue", "blueviolet", "deeppink", "indigo", "rosybrown", "greenyellow", "cyan", "fuchsia", "royalblue", "red", "green", "orange", "dodgerblue", "blueviolet", "deeppink", "indigo", "rosybrown", "greenyellow", "cyan", "fuchsia", "royalblue", "red", "green"]

x_label = r"$x$"
phi_pi_x_label = r"$\phi_{\pi}(x)$"
phi_k_x_label = r"$\phi_{K}(x)$"


def sum_rule(meson, x, a1, a2, a4):
    def C1(x, a):
        return 2 * a * x

    def C2(x, a):
        return 1 / 2 * (2 * x * (2 + a - 1) * C1(x, a) - (2 + 2 * a - 2))

    def C3(x, a):
        return 1 / 3 * (2 * x * (3 + a - 1) * C2(x, a) - (3 + 2 * a - 2) * C1(x, a))

    def C4(x, a):
        return 1 / 4 * (2 * x * (4 + a - 1) * C3(x, a) - (4 + 2 * a - 2) * C2(x, a))

    if meson == "pion":
        return (
            6
            * x
            * (1 - x)
            * (1 + C2(2 * x - 1, 3 / 2) * a2 + C4(2 * x - 1, 3 / 2) * a4)
        )

    if meson == "kaon":
        return (
            6
            * x
            * (1 - x)
            * (1 + C1(2 * x - 1, 3 / 2) * a1 + C2(2 * x - 1, 3 / 2) * a2)
        )


def DSE(x):
    return 18.2 * x * (1 - x) * (1 - 2.33 * np.sqrt(x * (1 - x)) + 1.79 * x * (1 - x))


def DSE_kaon():
    pix_x = [22, 66, 163, 235, 307, 402, 429, 547, 669, 829, 906, 960, 1033, 982, 1112, 1185, 1266, 1306, 1373, 1417, 1454, 1486, 1509, 1525]
    pix_y = [75, 204, 433, 557, 658, 751, 770, 840, 880, 906, 907, 905, 891, 904, 864, 821, 745, 692, 571, 469, 364, 249, 161, 80]

    x1 = 1537  # coor of x=1
    y1 = 685  # coor of y=1

    x = []
    y = []

    for val in pix_x:
        x.append((x1 - val) / x1)

    for val in pix_y:
        y.append(val / y1)

    x = np.array(x)
    y = np.array(y)

    return x, y


def lcda_large_pz_plot(meson, x_ls, mom_n_lic_da, large_mom_lic_da):
    ### replace all [:] with [202:301] for plot with tails ###

    ### sys error of large mom limit
    mom_sys_ls = []
    for idx in range(len(x_ls)):
        mom_sys = abs(
            large_mom_lic_da[idx].mean - mom_n_lic_da[idx].mean
        )  # system error
        mom_sys_ls.append(mom_sys)

    ### sys error of extrapolation
    if meson == "pion":
        ext_y = gv.load("dump/p_dif_ext_y")
    elif meson == "kaon":
        ext_y = gv.load("dump/k_dif_ext_y")

    ext_sys_ls = []
    for idx in range(len(x_ls)):
        ext_sys = abs(large_mom_lic_da[idx].mean - ext_y[idx])
        ext_sys_ls.append(ext_sys)

    ### sys error of continuum limit
    if meson == "pion":
        con_y = gv.load("dump/p_a06_y")
    elif meson == "kaon":
        con_y = gv.load("dump/k_a06_y")

    con_sys_ls = []
    for idx in range(len(x_ls)):
        con_sys = abs(large_mom_lic_da[idx].mean - con_y[idx])
        con_sys_ls.append(con_sys)

    ### sys error of mu
    if meson == "pion":
        mu_y = gv.load("dump/p_mu_y")
    elif meson == "kaon":
        mu_y = gv.load("dump/k_mu_y")

    mu_sys_ls = []
    for idx in range(len(x_ls)):
        mu_sys = abs(large_mom_lic_da[idx].mean - mu_y[idx])
        mu_sys_ls.append(mu_sys)

    y1 = np.array(
        [
            (
                large_mom_lic_da[id].mean
                + np.sqrt(
                    large_mom_lic_da[id].sdev ** 2
                    + mom_sys_ls[id] ** 2
                    + ext_sys_ls[id] ** 2
                    + con_sys_ls[id] ** 2
                    + mu_sys_ls[id] ** 2
                )
            )
            for id in range(len(large_mom_lic_da))
        ]
    )
    y2 = np.array(
        [
            (
                large_mom_lic_da[id].mean
                - np.sqrt(
                    large_mom_lic_da[id].sdev ** 2
                    + mom_sys_ls[id] ** 2
                    + ext_sys_ls[id] ** 2
                    + con_sys_ls[id] ** 2
                    + mu_sys_ls[id] ** 2
                )
            )
            for id in range(len(large_mom_lic_da))
        ]
    )

    x_ls = np.hstack((x_ls, np.array([1])))
    y1 = np.hstack((y1, np.array([0])))
    y2 = np.hstack((y2, np.array([0])))

    fig = plt.figure(figsize=fig_size_sq)
    ax = plt.axes(plt_axes_small)

    a1 = gv.gvar(-0.06, 0.03)  # sum rule
    a2 = gv.gvar(0.25, 0.15)  # same for pion and kaon
    a4 = gv.gvar(-0.015, 0.025)

    ax.fill_between(
        x_ls[:],
        [
            sum_rule(meson, x, a1, a2, a4).mean + sum_rule(meson, x, a1, a2, a4).sdev
            for x in x_ls
        ][:],
        [
            sum_rule(meson, x, a1, a2, a4).mean - sum_rule(meson, x, a1, a2, a4).sdev
            for x in x_ls
        ][:],
        color=color_list[1],
        label="Sum rule",
        alpha=0.4,
    )

    if meson == "pion":
        a2 = gv.gvar(0.101, 0.024)
        ope = [sum_rule(meson, x, 0, a2, 0) for x in x_ls]
    elif meson == "kaon":
        a1 = gv.gvar(-0.0533, 0.0034)
        a2 = gv.gvar(0.090, 0.019)
        ope = [sum_rule(meson, x, a1, a2, 0) for x in x_ls]

    ax.fill_between(
        x_ls[:],
        [val.mean + val.sdev for val in ope][:],
        [val.mean - val.sdev for val in ope][:],
        color=color_list[2],
        label="OPE",
        alpha=0.6,
    )

    if meson == "pion":
        ax.plot(x_ls[:], DSE(x_ls)[:], color="blue", label="DSE", linestyle="dashed")

    elif meson == "kaon":
        dse_x, dse_y = DSE_kaon()
        ax.plot(dse_x, dse_y, color="blue", label="DSE", linestyle="dashed")

    ax.fill_between(x_ls, y1, y2, color=color_list[0], alpha=0.5)

    ax.plot(
        x_ls,
        (y1 + y2) / 2,
        color=color_list[0],
        label="This work",
        linewidth=2,
        linestyle="dotted",
    )

    ax.plot(
        x_ls[:],
        [6 * x * (1 - x) for x in x_ls][:],
        color="red",
        linestyle="dashdot",
        label="Asymptotic",
    )  # only plot between 0 and 1

    ax.fill_between(
        np.linspace(-0.5, 0.1, 500),
        np.ones(500) * -1,
        np.ones(500) * 2,
        color="grey",
        alpha=0.2,
    )
    ax.fill_between(
        np.linspace(0.9, 1.5, 500),
        np.ones(500) * -1,
        np.ones(500) * 2,
        color="grey",
        alpha=0.2,
    )

    ## grey v band to cover fit region

    ax.axvline(0.5, color="green", linestyle="--")

    ax.axhline(0, color="k", linestyle="--")
    ax.set_xlabel(x_label, **fs_p_l)
    if meson == "pion":
        ax.set_ylabel(phi_pi_x_label, **fs_p_l)
    elif meson == "kaon":
        ax.set_ylabel(phi_k_x_label, **fs_p_l)
    ax.set_ylim([-0.19, 1.7])
    ax.set_xlim([0, 1])
    ax.legend(loc="lower center")
    ax.tick_params(direction="in", **ls_p_l)
    plt.savefig("lcda_Pz_to_infty_of_" + meson + ".pdf", transparent=True)
    plt.show()
